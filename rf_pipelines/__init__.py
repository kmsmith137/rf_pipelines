"""High-level docstring here!!  For more examples, see the examples/ dir"""

# Low-level streams and transform classes.
# See below for python interface.
from .rf_pipelines_c import wi_stream, wi_transform

# Helper routines for writing transforms in python.
from .utils import write_png, wi_downsample

# Library of transforms written in python.
from .transforms.plotter_transform import plotter_transform
from .transforms.frb_injector_transform import frb_injector_transform


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
