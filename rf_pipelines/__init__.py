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

The API for implementing a new transform in Python is non-obvious, since it's implicitly
defined by the C++ part of the library.  All Python documentation is in the docstring
for class 'py_wi_transform' below.

Everyone should feel free to commit new transforms to the rf_pipelines github repo.
Over time we can build up a library of building-block transforms which can be recycled
for rapid prototyping, e.g. if we want to experiment with different RFI removal schemes
by chaining together different RFI removers.

To write a transform (or stream) in C++, see comments in 'rf_pipelines.hpp'.  One loose end:
exporting C++ classes to python is currently kind of a mess, see comments in rf_pipelines_c.cpp.  
I hope to improve this soon!  In the meantime, feel free to email me if it's too cryptic.
"""


# These are the base classes for (written in C++ and exported to Python via the
# 'rf_pipelines_c' extension module).  
#
# If you want to write a new transform in python, I recommend subclassing the python
# class 'py_wi_transform' (see below), rather than subclassing the C++ class 'wi_transform'.
#
# For the API that you'll need to implement in order to write a new transform, see the
# docstring for class 'py_wi_transform'.

from .rf_pipelines_c import wi_stream, wi_transform

# Library of transforms written in python.
from .transforms.plotter_transform import plotter_transform
from .transforms.frb_injector_transform import frb_injector_transform

# Helper routines for writing transforms in python.
from .utils import write_png, wi_downsample


####################################################################################################
#
# Stream library.  All functions below return an object of class 'wi_stream'.  They are implemented
# by wrapping a function written in C++, and interfaced with python in rf_pipelines_c.cpp.
#
# FIXME currently, there is no way to write a stream in Python, but this will be fixed soon!


def chime_stream_from_filename(filename, nt_chunk=0):
    """
    Returns a weighted intensity stream (wi_stream) from a single CHIME hdf5 file.

    The 'filename' arg should be an hdf5 file containing CHIME intensity data.

    The 'nt_chunk' arg is the chunk size used internally when moving data from hdf5 file
    into the rf_pipelines buffer.  If unspecified or zero, it will default to a reasonable value.

    Note: a quick way to inspect a CHIME hdf5 file is using the 'ch-show-intensity-file' program,
    in the ch_frb_io github repo.
    """

    return rf_pipelines_c.make_chime_stream_from_filename(filename, nt_chunk)


def chime_stream_from_filename_list(filename_list, nt_chunk=0):
    """
    Returns a weighted intensity stream (wi_stream) from a sequence of CHIME hdf5 files.

    The 'filename_list' arg should be a list (or python generator) of hdf5 filenames.

    The 'nt_chunk' arg is the chunk size used internally when moving data from hdf5 file
    into the rf_pipelines buffer.  If unspecified or zero, it will default to a reasonable value.

    Note: a quick way to inspect a CHIME hdf5 file is using the 'ch-show-intensity-file' program,
    in the ch_frb_io github repo.
    """

    return rf_pipelines_c.make_chime_stream_from_filename_list(filename_list, nt_chunk)


def chime_stream_from_acqdir(dirname, nt_chunk=0):
    """
    Returns a weighted intensity stream (wi_stream) from an acquisition directory containing CHIME hdf5 files.
    The directory is scanned for filenames of the form NNNNNNNN.h5, where N=[0,9].
    
    The 'nt_chunk' arg is the chunk size used internally when moving data from hdf5 file
    into the rf_pipelines buffer.  If unspecified or zero, it will default to a reasonable value.

    Note: a quick way to inspect a CHIME hdf5 file is using the 'ch-show-intensity-file' program,
    in the ch_frb_io github repo.
    """

    return rf_pipelines_c.make_chime_stream_from_acqdir(dirname, nt_chunk)


def psrfits_stream(filename):
    """Returns a weighted intensity stream (wi_stream) from a single PSRFITS source file."""

    return rf_pipelines_c.make_psrfits_stream(filename)


def gaussian_noise_stream(nfreq, nt_tot, freq_lo_MHz, freq_hi_MHz, dt_sample, sample_rms=1.0, nt_chunk=0):
    """
    Returns a weighted intensity stream (wi_stream) which simulates Gaussian random noise for each frequency channel and time sample.
    
    The 'nt_tot' arg is the total length of the stream, in samples.
    The 'dt_sample' arg is the length of a sample in seconds.
    The 'sample_rms' arg is the Gaussian RMS of a single (freq_channel, time_sample) pair.

    The 'nt_chunk' arg is the chunk size used internally when moving data into the rf_pipelines buffer.
    If unspecified or zero, it will default to a reasonable value.
    """

    return rf_pipelines_c.make_gaussian_noise_stream(nfreq, nt_tot, freq_lo_MHz, freq_hi_MHz, dt_sample, sample_rms, nt_chunk)


####################################################################################################
#
# Transform library


def simple_detrender(nt_chunk):
    """
    Returns a transform object (wi_transform) which detrends the data.
    Simplest possible detrender: just divides the data into chunks and subtracts the mean in each chunk.
    """

    return rf_pipelines_c.make_simple_detrender(nt_chunk)


def bonsai_dedisperser(config_hdf5_filename, output_hdf5_filename, ibeam=0):
    """
    Returns a "transform" (object of class wi_transform) which doesn't actually modify the data,
    it just runs the bonsai dedisperser.  The output is a stream of coarse-grained triggers
    which are written to an output hdf5 file.  The dedisperser must be initialized from a
    config hdf5 file produced with the program 'bonsai-mkweight' in the bonsai github repo.

    Note that the program 'bonsai-plot-triggers.py' in the bonsai github repo may be useful
    for quick visual inspection of the bonsai output.

    The 'ibeam' argument determines the assignment of threads to cores and can probably
    be zero except in special situations.

    FIXME 1: Currently the only trigger "processing" which can be done is writing the triggers
    to an hdf5 file for later analysis.  It would be better if we could use bonsai's python
    interface, which allows the dedisperser process_triggers() callback to be written in
    python.  Right now, it's not possible to use bonsai's python interface with rf_pipelines!

    FIXME 2: Currently the dedisperser must be initialized from a config hdf5 file (rather than
    the simpler config text file) since we use analytic weights to normalize the triggers.
    Since the analytic weights are only correct for unit-variance noise, the trigger normalization
    will be wrong for a real experiment, and the triggers won't be meaningfully normalized to
    "sigmas".  All of this is just a placeholder until Monte Carlo trigger variance estimation
    is implemented in bonsai.
    """
    
    return rf_pipelines_c.make_bonsai_dedisperser(config_hdf5_filename, output_hdf5_filename, ibeam)


####################################################################################################
#
# Python interface for defining new transforms.


class py_wi_transform(wi_transform):
    """
    The API for implementing a new wi_transform in Python is non-obvious, since it's implicitly
    defined by the C++ part of the library.  All Python documentation is in this docstring!

    Define a subclass of 'py_wi_transform' which overrides the following methods:
    

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


    end_substream()

       Counterpart to start_substream()
    """

    def set_stream(self, stream):
        pass

    def start_substream(self, isubstream, t0):
        pass

    def process_chunk(self, t0, t1, intensity, weights, pp_intensity, pp_weights):
        pass

    def end_substream(self):
        pass
