"""
Currently all that's here is the 'chime_file_writer'.  This is a thin wrapper around a C++
implementation.  See chime_file_writer.cpp, and python linkage in rf_pipelines_c.cpp.
"""

from rf_pipelines import rf_pipelines_c


def chime_file_writer(filename, clobber=False, bitshuffle=2, nt_chunk=0):
    """
    This is a pseudo-transform which doesn't actually modify the data, it just writes it to a file in
    CHIME hdf5 format.  (For now, the entire stream is written to a single file, I'll generalize later
    to break the stream into multiple files.)

    If 'clobber' is false, and the target file already exists, an exception will be thrown rather than clobbering the old file.
    If 'nt_chunk' is set to zero, a default chunk size will be chosen.

    The meaning of the 'bitshuffle' arg is:
       0 = no compression
       1 = try to compress, but if plugin fails then just write uncompressed data instead
       2 = try to compress, but if plugin fails then print a warning and write uncompressed data instead
       3 = compression mandatory

    Note: a quick way to inspect a CHIME hdf5 file is using the 'ch-show-intensity-file' and 'ch-plot-intensity-file'
    programs, in the ch_frb_io github repo.
    """

    return rf_pipelines_c.make_chime_file_writer(filename, clobber, bitshuffle, nt_chunk)

