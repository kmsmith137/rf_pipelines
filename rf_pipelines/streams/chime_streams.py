"""
CHIME streams.

Everything here is a thin wrapper around a C++ implementation, written in chime_streams.cpp
and exported to Python via rf_pipelines_c.cpp.
"""

from rf_pipelines import rf_pipelines_c
from os import listdir
from h5py import File
from sys import stderr


def chime_stream_from_filename(filename, nt_chunk=0, noise_source_align=0):
    """
    Returns a weighted intensity stream (wi_stream) from a single CHIME hdf5 file.

    The 'filename' arg should be an hdf5 file containing CHIME intensity data.

    The 'nt_chunk' arg is the chunk size used internally when moving data from hdf5 file
    into the rf_pipelines buffer.  If unspecified or zero, it will default to a reasonable value.

    If 'noise_source_align' is nonzero, then it should be equal to the DETRENDER chunk size 
    (not the chime_file_stream nt_chunk).  In this case, the stream will align the noise source 
    edges with the detrender chunks, by discarding initial data if necessary.

    Note: a quick way to inspect a CHIME hdf5 file is using the 'ch-show-intensity-file' and 'ch-plot-intensity-file'
    programs, in the ch_frb_io github repo.
    """

    return rf_pipelines_c.make_chime_stream_from_filename(filename, nt_chunk, noise_source_align)


def chime_stream_from_filename_list(filename_list, nt_chunk=0, noise_source_align=0):
    """
    Returns a weighted intensity stream (wi_stream) from a sequence of CHIME hdf5 files.

    The 'filename_list' arg should be a list (or python generator) of hdf5 filenames.

    The 'nt_chunk' arg is the chunk size used internally when moving data from hdf5 file
    into the rf_pipelines buffer.  If unspecified or zero, it will default to a reasonable value.

    If 'noise_source_align' is nonzero, then it should be equal to the DETRENDER chunk size 
    (not the chime_file_stream nt_chunk).  In this case, the stream will align the noise source 
    edges with the detrender chunks, by discarding initial data if necessary.

    Note: a quick way to inspect a CHIME hdf5 file is using the 'ch-show-intensity-file' program,
    in the ch_frb_io github repo.
    """

    return rf_pipelines_c.make_chime_stream_from_filename_list(filename_list, nt_chunk, noise_source_align)


def chime_stream_from_acqdir(dirname, nt_chunk=0, noise_source_align=0):
    """
    Returns a weighted intensity stream (wi_stream) from an acquisition directory containing CHIME hdf5 files.
    The directory is scanned for filenames of the form NNNNNNNN.h5, where N=[0,9].
    
    This routine should be used with caution, e.g. if pointed to /data/pathfinder/16-06-21 on chimer
    it will try to analyze 43GB of pathfinder data as a single stream!
    
    The 'nt_chunk' arg is the chunk size used internally when moving data from hdf5 file
    into the rf_pipelines buffer.  If unspecified or zero, it will default to a reasonable value.

    If 'noise_source_align' is nonzero, then it should be equal to the DETRENDER chunk size 
    (not the chime_file_stream nt_chunk).  In this case, the stream will align the noise source 
    edges with the detrender chunks, by discarding initial data if necessary.

    Note: a quick way to inspect a CHIME hdf5 file is using the 'ch-show-intensity-file' program,
    in the ch_frb_io github repo.
    """

    return rf_pipelines_c.make_chime_stream_from_acqdir(dirname, nt_chunk, noise_source_align)


def chime_network_stream(udp_port=0, beam_id=0):
    """
    CHIME network stream.  Receives UDP packets in "CHIME L0-L1 format".
    
    This interface is less general than the low-level C++ interface in ch_frb_io: 
    only one beam can be received, and not all boolean options are supported.
    
    If the 'udp_port' argument is zero, then the default chimefrb port will be used.
    """

    return rf_pipelines_c.make_chime_network_stream(udp_port, beam_id)


def get_start_time(file):
    f = File(file, 'r')
    timestamp_array = f['index_map']['time'][:]
    return timestamp_array[0]


def chime_stream_from_times(dirname, t0, t1, nt_chunk=0, noise_source_align=0):
    files = listdir(dirname)

    # Make sure all the files are actually data files
    filter(lambda x: x[-3:] == '.h5', files)

    # The files _should_ be ordered, but check just in case
    for i in range(len(files) - 1):
        assert int(files[i+1][:-3]) > int(files[i][:-3])

    # Search for the first file that has a timestamp that is greater than t0. Take the previous one.
    # Then, search for the first file with a timestamp greater than t1 and take that file.
    start, end = -1, -1

    for i in range(len(files)):
        if get_start_time(files[i]) > t0:
            start = i - i
            break

    for i in range(start, len(files)):
        if get_start_time(files[i]) > t1:
            end = i + 1   # i + 1 since indexing excludes final index
            break

    if start == -1 or end == -1:
        print >> stderr, 'The indices specified were not found.'
        return

    return chime_stream_from_filename_list(files[start:end], nt_chunk=ntchunk, noise_source_align=noise_source_align)
