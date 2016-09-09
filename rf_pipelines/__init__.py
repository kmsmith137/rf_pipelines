"""
rf_pipelines: plugin-based radio astronomy pipelines.

There are two fundamental classes, wi_streams ("weighted intensity stream")
and wi_transforms ("weighted intensity transform").  

A wi_stream incrementally generates chunks of intensity data which is channelized,
meaning that the data lives in a (freq channel, time sample) 2D array, and weighted, 
meaning that there is a parallel array of floating-point weights.

A wi_transform operates on weighted intensity data, modifying it in place.  For
example, a detrender or an RFI removal stage could be a wi_transform.  There are
also pseudo-transforms which process the data without actually modifying it, for
example plotters or dedispersers.  Transforms run in sequence: each transform's
output is the input to the next transform.

The best way to see how to use the library is to look at the examples in the
examples/ directory:

  example1_toy/    illustrates making noise+FRB sim, writing transforms in python
  example2_gbncc/  bare-bones gbncc example
  example3_chime/  bare-bones chime example

Transforms can be written either in C++ or Python.  A goal for the library is
complete interoperability: it should be possible to write streams and transforms
in either C++ or Python and mix freely.  Currently this is implemented for transforms,
but streams can only be written in C++ (I hope to improve this soon).  If all streams
and transforms are written in C++, the pipeline can be run from python with no
overhead, allowing python to be the configuration language for optimized C++ code.

If you want to write a new transform in python, the API that you'll need to implement
is non-obvious, since it is implicitly defined by the C++ part of the library.  All
documentation is in the docstring for class 'py_wi_transform' below!

To write a transform (or stream) in C++, see comments in 'rf_pipelines.hpp'.  One loose end:
exporting C++ classes to python is currently kind of a mess, see comments in rf_pipelines_c.cpp.  
I hope to improve this soon!  In the meantime, feel free to email me if it's too cryptic.

Everyone should feel free to commit new transforms to the rf_pipelines github repo.
Over time we can build up a library of building-block transforms which can be recycled
for rapid prototyping, e.g. if we want to experiment with different RFI removal schemes
by chaining together different RFI removers.

Here is a list of all streams and transforms currently available.  For documatation, see
the individual docstrings.

Streams:
   chime_stream_from_acqdir()
   chime_stream_from_filename()
   chime_stream_from_filename_list()
   chime_network_stream()
   gaussian_noise_stream()
   psrfits_stream()

Transforms:

   badchannel_mask()          masks list of frequency ranges specified in external file
   bonsai_dedisperser()       runs data through bonsai dedisperser
   chime_file_writer()        write stream to a single file in CHIME hdf5 format
   chime_packetizer()         send stream over network, to a chime_network_stream running on another machine
   clipper_transform()        masks data based on intensity values
   frb_injector_transform()   simulates an FRB (currently S/N calculation only works for toy noise models)
   kurtosis_filter()          masks data based on kurtosis
   plotter_transform()        makes waterfall plots at a specified place in the pipeline, very useful for debugging
   RC_detrender()             exponential detrender, with bidirectional feature intended to remove "step-like" features
   simple_detreneder()        really boneheaded detrending algorithm (better detrending is available in python, but it's slow!)
   std_dev_filter()           masks data based on variance
   thermal_noise_weight()     applies optimal weighting assuming flat gains and variance proportional to intensity
"""


import sys
import numpy as np

# Not sure why, but this import has to be in the toplevel module,
# or png writing will sometimes segfault!
try:
    import PIL.Image
except:
    print >>sys.stderr, "rf_pipelines: import PIL.Image failed; many things will work but plotting will fail"


# The 'wi_stream' and 'wi_transform' base classes are subclassed to define streams and transforms.
# These classes are written in C++ and exported to Python via 'rf_pipelines_c.cpp'
#
# If you want to write a new transform in python, the API that you'll need to implement
# is non-obvious, since it is implicitly defined by the C++ part of the library.  All
# documentation is in the docstring for class 'py_wi_transform' below!

from .rf_pipelines_c import wi_stream, wi_transform, wi_run_state



####################################################################################################
#
# Python interface for defining new transforms.


class py_wi_transform(wi_transform):
    """
    The API for implementing a new wi_transform in Python is non-obvious, since it's implicitly
    defined by the C++ part of the library.  All Python documentation is in this docstring!

    To write a transform in Python, I recommend subclassing the python class 'py_wi_transform' 
    instead of the C++ class 'wi_transform'.  An object of type 'py_wi_transform' can be inserted
    into a pipeline using the usual run() syntax:

       # Here, 's' is an object of type wi_stream, and each t_i is an object of type wi_transform
       s.run([t1,t2,...,tN], outdir='.', noisy=True, clobber=True)

    The 'py_wi_transform' constructor may want to initialize self.name, a "transform name" string
    which ends up in rf_pipelines.json.

    Your subclass of 'py_wi_transform' should override the following methods:
    

    set_stream(self, stream)

        This is called once, at the beginning of a pipeline run.  The 'stream' argument is an
        object of class wi_stream, and the following members will be defined:
           stream.nfreq         number of frequency channels in stream
           stream.freq_lo_MHz   lower limit of frequency band
           stream.freq_hi_MHz   upper limit of frequency band
           stream.dt_sample     sample length in seconds
           stream.nt_maxwrite   internal chunk size of stream, this probably won't be useful.

        The job of set_stream() is to initialize the following four members.
           self.nfreq           number of frequency channels
           self.nt_chunk        chunk size for process_chunk(), see below
           self.nt_prepad       prepad size for process_chunk(), see below
           self.nt_postpad      postpad size for process_chunk(), see below

        These four members can be initialized either in the transform constructor or in set_stream(),
        but must be defined before set_stream() returns.


    start_substream(self, isubstream, t0)

        Each stream can divide its output into one or more "substreams" which are bracketed with
        start_substream() and end_substream() calls.  Currently I don't use the "substreaming" feature
        much: all streams just define a single substream, and not all transforms support multiple substreams.

        However, I anticipate it being a useful feature when we implement real-time network streams,
        since we'll want a way to finalize state when the correlator goes down (or repoints) and
        restart when it comes back.

        The 'isubstream' arg is 0 for the first substream, 1 for the second substream, etc.
        The 't0' arg is the initial time of the substream in seconds, relative to an arbitrary stream-defined origin.


    process_chunk(self, t0, t1, intensity, weights, pp_intensity, pp_weights)
    
        This routine is called to deliver chunks of data to the transform.  Each transform defines
        three buffer sizes (see above): nt_chunk, nt_prepad, and nt_postpad.

        Each call to process_chunk() is responsible for processing one "chunk" of data with 2D shape
        (self.nfreq, self.nt_chunk).  The 'intensity' and 'weights' arguments are floating-point arrays
        with this shape.
        
        Important: the frequency axis of the arrays is ordered from highest frequency to lowest!
        This is the same ordering used by 'bonsai', and in both GBNCC and CHIME data, so this ordering
        seemed most convenient.

        The 't0' and 't1' args are timestamps in seconds at the beginning and end of the chunk.
        (To be precise, t0 is the start time of the first sample in the chunk, and t1 is the end
        time of the last sample in the chunk.)  This information is mostly redundant since the
        initial time of the substream is passed in start_substream(), and the time sample length
        is available in set_stream().  However, I'm anticipating that for long-running streams it 
        may be useful to allow for a small amount of timestamp "drift" via the t0/t1 args.

        Some transforms will need to do read-only inspection of data outside the chunk.  For example,
        a detrending transform may need to look at values of the data a little bit before and after
        the chunk, in order to detrend data in the chunk.  

        A transform which needs to "see into the future" can set self.nt_postpad to a nonzero value.  
        In this case, 'intensity' and 'weights' will have shape (self.nfreq, self.nt_chunk + self.nt_postpad).
        Important: the transform is only allowed to modify the first 'self.nt_chunk' samples!

        A transform which needs to "see into the past" can set self.nt_prepad to a nonzero value.
        In this case, the prepadding data is passed as separate arrays, via the 'pp_intensity' and
        'pp_weights' args which have shapes (self.nfreq, self.nt_prepad).  This is different from
        the postpadded case, where the extra data is passed by extending the chunk array.  If
        self.nt_prepad is zero, then the 'pp_intensity' and 'pp_weights' arguments are None.

        Each transform can initialize its nt_chunk, nt_prepad, and nt_postpad independently of
        the other transforms, and the rf_pipelines library will handle the necessary buffering 
        and rechunking.

        Transforms should make sure to handle the case where there are many zeroes in the 'weights'
        array, including the extreme case where the weights are all zeros.  One situation where
        many zeros arise is at the end of a stream, where the total stream length may not be
        a multiple of nt_chunk, and so the rf_pipelines library will append zero-weight data.


    end_substream(): counterpart to start_substream() above.
    """

    def __init__(self, name=None):
        """Base class constructor is just a reminder to set self.name, and supplies a default value."""
        self.name = name if (name is not None) else self.__class__.__name__

    def set_stream(self, stream):
        pass

    def start_substream(self, isubstream, t0):
        pass

    def process_chunk(self, t0, t1, intensity, weights, pp_intensity, pp_weights):
        pass

    def end_substream(self):
        pass

    def __str__(self):
        return self.name if (self.name is not None) else self.__class__.__name__

    # Note that py_wi_transform inherits the following methods from the C++ base class 
    # 'wi_transform'.  (See method docstrings for more info.)
    #
    #   add_plot_group(name, nt_per_pix, ny) -> integer group_id
    #   add_plot(basename, it0, nt, nx, ny, group_id=0) -> string filename
    #   add_file(basename) -> string filename




####################################################################################################
#
# Python interface for defining new transforms.


class py_wi_stream(wi_stream):
    """
    The API for implementing a new wi_stream in Python is non-obvious, since it's implicitly
    defined by the C++ part of the library.  All Python documentation is in this docstring!

    The wi_stream class defines a method

        run(self, transform_list, outdir='.', noisy=True, clobber=True, return_json=False)

    which is called to run a pipeline.  Here, 

       - 'transform_list' is a list (or generator) of objects of type wi_transform (including
         its subclass py_wi_transform).

       - 'outdir' is the rf_pipelines output directory, where the rf_pipelines json file will
         be written, in addition to other transform-specific output files such as plots. 

         If 'outdir' is None or an empty string, then the json file will not be written,
         and any transform which tries to write an output file (such as a plotter_transform)
         will throw an exception.
    
       - If 'clobber' is False, then an exception will be thrown if the pipeline tries to
         overwrite an old rf_pipelines.json file.

       - If 'return_json' is True, then the return value from run() will be the rf_pipelines
         json output (i.e. same data which is written to rf_pipelines.json).
         
         A kludge: eventually, the run() return value will be a json object, but for now it returns 
         the string representation, which can be converted to a json object by calling json.loads().

    The wi_stream class defines the following members:
         nfreq            number of frequency channels
         freq_lo_MHz      lowest frequency in band (e.g. 400.0 for CHIME)
         freq_hi_MHz      highest frequency in band (e.g. 800.0 for CHIME)
         dt_sample        length of a sample in seconds
         nt_maxwrite      block size of stream (defined as max number of time samples per call to setup_write())

    These can be initialized either in the stream constructor, or in the optional stream_start() method.
    In some cases, it's convenient to defer initialization until stream_start(), since they may not be
    known until the stream starts running.  An example is a network stream, which may not get the value
    of dt_sample (say) until the first packet is received.

    Note: don't set nt_maxwrite to an excessively large value, since there is an internal
    buffer of approximate size (24 bytes) * nfreq * nt_maxwrite.

    To write a new stream in Python, I recommend subclassing the python class 'py_wi_stream' 
    instead of the C++ class 'wi_stream'.  Your subclass of 'py_wi_transform' should initialize
    the above fields { nfreq, ..., nt_maxwrite } in its constructor, and it should define the
    following method:

        stream_body(self, run_state)

    This method will consist of a loop which moves blocks of data from some source (a file or 
    network connection) into the rf_pipelines ring buffer.  The 'run_state' argument is an
    object of class 'wi_run_state'.  This is a low-level C++ class which defines the following
    methods:

        run_state.start_substream(t0)
        run_state.write(intensity_arr, weight_arr, t0=None)
        run_state.end_substream()

    The run_state.write() method is called by stream_body() to copy data into the rf_pipelines
    ring buffer.  Its 'intensity_arr' and 'weight_arr' arguments are 2D arrays of shape (nfreq, nt)
    where 'nt' is a number of samples satisfying 0 < nt <= stream.nt_maxwrite.  The 't0' arg
    is the start time of the block in seconds.  Note that t0 is optional since it can be inferred from 
    the substream start time, the value of wi_stream::dt_sample, and the number of samples written 
    so far.  However, it may be useful to specify t0 occasionally in order to keep track of slow 
    timestamp drifts over time.  For example in the chimefrb pipeline, the intensity samples 
    always correspond to a fixed number of FPGA counts, and the FPGA clock drifts on long timescales.

    The run_state.start_substream() and run_state.end_substream() methods can be used to divide the
    stream into multiple substreams.  The downstream transforms should reset state between substreams.
    The 't0' arg to start_substream() is the initial time of the substream in seconds, relative to an 
    arbitrary stream-defined origin.

    At the moment the "multiple-substream" feature isn't very well-supported, so it's probably best for all 
    streams to represent their data as a single substream. However, I anticipate it being a useful feature 
    when we implement real-time network streams, since we'll want a way to finalize state when the correlator 
    goes down (or repoints) and restart when it comes back.

    Sumarizing, wi_stream.stream_body() should look something like this.

       def stream_body(self, run_state):
           for ...:                          # outer loop over substreams
               run_state.start_substream()
               for ...:
                   run_state.write()         # inner loop over data blocks
               run_state.end_substream()    

    For an example, check out the 'frb_olympics' github repo, and see 'rerunnable_gaussian_noise_stream'
    in frb_olympics/frb_olympics.py.  (Unfortunately, there's no example of a python stream in the rf_pipelines
    repo itself, so I have to recommend an example in a different repo!)

    Comment: the python stream API has an extra copy relative to the C++ API, so python streams may be a little
    slower than C++ streams, but it's hard to do anything about this!
    """

    def __init__(self, nfreq, freq_lo_MHz, freq_hi_MHz, dt_sample, nt_maxwrite):
        # It's not necessary for the subclass constructor to call this base class constructor,
        # but the five members below must be initialized either in the constructor or stream_start().
        self.nfreq = nfreq
        self.freq_lo_MHz = freq_lo_MHz
        self.freq_hi_MHz = freq_hi_MHz
        self.dt_sample = dt_sample
        self.nt_maxwrite = nt_maxwrite


    def stream_start(self):
        """
        This optional method can be defined by the stream class, if it's convenient to defer
        initializations until the stream starts running.  An example is a network stream, which 
        may not get the value of dt_sample (say) until the first packet is received.
        """
        pass


    def stream_body(self, run_state):
        raise RuntimeError('Subclass of py_wi_stream must override stream_body().')


####################################################################################################
#
# Library of built-in streams and transforms


# Streams (all implemented in C++)

from .streams.chime_streams import \
    chime_stream_from_filename, \
    chime_stream_from_filename_list, \
    chime_stream_from_acqdir, \
    chime_network_stream

from .streams.psrfits_stream import psrfits_stream
from .streams.gaussian_noise_stream import gaussian_noise_stream

# Transforms (some implemented in C++, others in python)

from .transforms.chime_packetizer import chime_packetizer
from .transforms.chime_transforms import chime_file_writer
from .transforms.plotter_transform import plotter_transform
from .transforms.simple_detrender import simple_detrender
from .transforms.bonsai_dedisperser import bonsai_dedisperser
from .transforms.frb_injector_transform import frb_injector_transform
from .transforms.badchannel_mask import badchannel_mask
from .transforms.clipper_transform import clipper_transform
from .transforms.legendre_detrender import legendre_detrender
from .transforms.kurtosis_filter import kurtosis_filter
from .transforms.std_dev_filter import std_dev_filter
from .transforms.thermal_noise_weight import thermal_noise_weight
from .transforms.RC_detrender import RC_detrender

# Helper routines for implementing new transforms in python.

from .utils import write_png, wi_downsample, upsample, tile_arr

# Grouping code to handle output of bonsai_dedisperser
from grouper import group_bonsai_output
