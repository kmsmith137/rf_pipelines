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

   chime_stream_from_filename(filename, nt_chunk=0)
   chime_stream_from_filename_list(filename_list, nt_chunk=0)
   chime_stream_from_acqdir(dirname, nt_chunk=0)
   psrfits_stream(fits_filename)
   gaussian_noise_stream(nfreq, nt_tot, freq_lo_MHz, freq_hi_MHz, dt_sample, sample_rms=1.0, nt_chunk=0)

Transforms:

   simple_detrender(nt_chunk)
   plotter_transform(img_prefix, img_nfreq, img_nt, downsample_nt=1, nt_chunk=0)
   frb_injector_transform(snr, undispersed_arrival_time, dm, intrinsic_width=0.0, sm=0.0, spectral_index=0.0, sample_rms=1.0, nt_chunk=1024)
   bonsai_dedisperser(config_hdf5_filename, output_hdf5_filename, ibeam=0)
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

from .rf_pipelines_c import wi_stream, wi_transform



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

       s.run([t1,t2,...,tN])    # where s is an object of type wi_stream
                                # and each t_i is an object of type wi_transform

    Your subclass of 'py_wi_transform' should override the following methods:
    

    set_stream(self, stream)

        This is called once, at the beginning of a pipeline run.  The 'stream' argument is an
        object of class wi_stream, and the following members will be defined:
           stream.nfreq         number of frequency channels in stream
           stream.freq_lo_MHz   lower limit of frequency band
           stream.freq_hi_MHz   upper limit of frequency band
           stream.dt_sample     sample length in seconds
           stream.nt_maxwrite   internal chunk size of stream, this probably won't be useful.

        The job of set_stream() is to initialize the following four members:
           self.nfreq           number of frequency channels
           self.nt_chunk        chunk size for process_chunk(), see below
           self.nt_prepad       prepad size for process_chunk(), see below
           self.nt_postpad      postpad size for process_chunk(), see below

        These four members can be initialized either in the transform constructor or in set_stream(),
        but must be defined before set_stream() returns.


    start_substream(self, isubstream, t0)

        Each stream can divide its output into one or more "substreams" which are delineated with
        start_substream() and end_substream() calls.  Currently I don't use this feature much: all
        streams just define a single substream, and not all transforms support multiple substreams.

        However, I anticipate it being a useful feature when we implement real-time network streams,
        since we'll want a way to finalize state when the correlator goes down (or repoints) and
        restart when it comes back.

        The 'isubstream' arg is 0 for the first substream, 1 for the second substream, etc.
        The 't0' arg is the initial time of the substream in seconds, relative to an arbitrary origin.


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

    def set_stream(self, stream):
        pass

    def start_substream(self, isubstream, t0):
        pass

    def process_chunk(self, t0, t1, intensity, weights, pp_intensity, pp_weights):
        pass

    def end_substream(self):
        pass


####################################################################################################
#
# Library of built-in streams and transforms


# Streams (all implemented in C++, there is currently no interface for writing streams in python)

from .streams.chime_streams import \
    chime_stream_from_filename, \
    chime_stream_from_filename_list, \
    chime_stream_from_acqdir

from .streams.psrfits_stream import psrfits_stream
from .streams.gaussian_noise_stream import gaussian_noise_stream

# Transforms (some implemented in C++, others in python)

from .transforms.chime_transforms import chime_file_writer
from .transforms.plotter_transform import plotter_transform
from .transforms.simple_detrender import simple_detrender
from .transforms.bonsai_dedisperser import bonsai_dedisperser
from .transforms.frb_injector_transform import frb_injector_transform
from .transforms.badchannel_mask import badchannel_mask

# Helper routines for implementing new transforms in python.

from .utils import write_png, wi_downsample
