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
    Returns a weighted intensity stream from a single CHIME hdf5 file.

    The 'filename' arg should be an hdf5 file containing CHIME intensity data.

    The 'nt_chunk' arg is the chunk size used internally when moving data from hdf5 file
    into the rf_pipelines buffer.  If unspecified or zero, it will default to a reasonable value.

    Note: a quick way to inspect a CHIME hdf5 file is using the 'ch-show-intensity-file' program,
    in the ch_frb_io github repo.
    """

    return rf_pipelines_c.make_chime_stream_from_filename(filename, nt_chunk)


def chime_stream_from_filename_list(filename_list, nt_chunk=0):
    """
    Returns a weighted intensity stream from a sequence of CHIME hdf5 files.

    The 'filename_list' arg should be a list (or python generator) of hdf5 filenames.

    The 'nt_chunk' arg is the chunk size used internally when moving data from hdf5 file
    into the rf_pipelines buffer.  If unspecified or zero, it will default to a reasonable value.

    Note: a quick way to inspect a CHIME hdf5 file is using the 'ch-show-intensity-file' program,
    in the ch_frb_io github repo.
    """

    return rf_pipelines_c.make_chime_stream_from_filename_list(filename_list, nt_chunk)


def chime_stream_from_acqdir(dirname, nt_chunk=0):
    """
    Returns a weighted intensity stream from an acquisition directory containing CHIME hdf5 files.
    The directory is scanned for filenames of the form NNNNNNNN.h5, where N=[0,9].
    
    The 'nt_chunk' arg is the chunk size used internally when moving data from hdf5 file
    into the rf_pipelines buffer.  If unspecified or zero, it will default to a reasonable value.

    Note: a quick way to inspect a CHIME hdf5 file is using the 'ch-show-intensity-file' program,
    in the ch_frb_io github repo.
    """

    return rf_pipelines_c.make_chime_stream_from_acqdir(dirname, nt_chunk)

