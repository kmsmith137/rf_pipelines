"""
CHIME streams.

Everything here is a thin wrapper around a C++ implementation, written in chime_streams.cpp
and exported to Python via rf_pipelines_c.cpp.
"""

from rf_pipelines import rf_pipelines_c

# For chime_stream_from_times
from h5py import File
from os import listdir
from os.path import join 
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


def get_times(file):
    """Returns [start_time, end_time] for an .h5 file."""
    print 'Examining', file
    f = File(file, 'r')
    timestamp_array = f['index_map']['time'][:]
    return timestamp_array[0], timestamp_array[-1]


def chime_stream_from_times(dirname, t0, t1, nt_chunk=0, noise_source_align=0):
    """Calls chime_stream_from_filename_list for a specific time range from a directory.
    Helpful for re-running small subsections of acquisitions as seen on the web viewer
    (which displays the starting timestamp of each plot)."""

    # Get, filter, and sort the files
    files = listdir(dirname)
    files.sort()
    filter(lambda x: x[-3:] == '.h5', files)
    files_with_paths = [join(dirname, file) for file in files]

    assert t0 < t1, 'First time index must be less than second time index.'
    assert t0 < get_times(files_with_paths[-1])[1], \
        'First time index must be less than final time index for acquisition'
    assert t1 > get_times(files_with_paths[0])[0], \
        'Second time index must be greater than the first time index for the acquisition'

    # Search for start index
    start, end = -1, -1
    for i in range(len(files_with_paths) - 1):
        if t0 < get_times(files_with_paths[i + 1])[0]:
            start = i
            break
    if start == -1:
        # Handles case in which start is in the final data file
        start = len(files_with_paths) - 1

    # Search for end index
    for i in range(start, len(files_with_paths)):
        if t1 < get_times(files_with_paths[i])[1]:
            end = i + 1
            break
    if end == -1:
        # Handles case in which t1 exceeds the range of the files (or is equal to final timestamp)
        end = len(files_with_paths)

    print 'Indexing the files in', dirname, 'for file indices', start, 'to', end - 1, 'inclusive.'

    return chime_stream_from_filename_list(files_with_paths[start:end], nt_chunk=nt_chunk, noise_source_align=noise_source_align)
