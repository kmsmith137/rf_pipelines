from rf_pipelines.rf_pipelines_c import chime_stream_from_filename_list

# For chime_stream_from_times
from h5py import File
from os import listdir
from os.path import join 
from sys import stderr


def get_times(file, noisy=True):
    """Returns [start_time, end_time] for an .h5 file."""

    if noisy:
        print 'Examining', file

    f = File(file, 'r')
    timestamp_array = f['index_map']['time'][:]
    return timestamp_array[0], timestamp_array[-1]


def chime_stream_from_times(dirname, t0, t1, nt_chunk=0, noise_source_align=0, noisy=True):
    """Calls chime_stream_from_filename_list for a specific time range from a directory.
    Helpful for re-running small subsections of acquisitions as seen on the web viewer
    (which displays the starting timestamp of each plot)."""

    # Get, filter, and sort the files
    files = listdir(dirname)
    files.sort()
    filter(lambda x: x[-3:] == '.h5', files)
    files_with_paths = [join(dirname, file) for file in files]

    assert t0 < t1, 'First time index must be less than second time index.'
    
    # FIXME this assertion assumes a complete acquistion, meaning that 
    # the very last file has been written completely to disk. This is not
    # true when files are being copied over. In this case, the last file 
    # is still inaccessible, hence h5py.File fails.
    #assert t0 < get_times(files_with_paths[-1],noisy)[1], \
    #    'First time index must be less than final time index for acquisition'

    assert t1 > get_times(files_with_paths[0],noisy)[0], \
        'Second time index must be greater than the first time index for the acquisition'

    # Search for start index
    start, end = -1, -1
    for i in range(len(files_with_paths) - 1):
        if t0 < get_times(files_with_paths[i+1],noisy)[0]:
            start = i
            break
    if start == -1:
        # Handles case in which start is in the final data file
        start = len(files_with_paths) - 1

    # Search for end index
    for i in range(start, len(files_with_paths)):
        if t1 < get_times(files_with_paths[i],noisy)[1]:
            end = i + 1
            break
    if end == -1:
        # Handles case in which t1 exceeds the range of the files (or is equal to final timestamp)
        end = len(files_with_paths)

    print 'Indexing the files in', dirname, 'for file indices', start, 'to', end - 1, 'inclusive.'

    return chime_stream_from_filename_list(files_with_paths[start:end], nt_chunk=nt_chunk, noise_source_align=noise_source_align)
