"""
Simple detrender.  This is a thin wrapper around a C++ implementation.
See simple_detrender.cpp, and python linkage in rf_pipelines_c.cpp.
"""

from rf_pipelines import rf_pipelines_c


def simple_detrender(nt_chunk):
    """
    Returns a transform object (wi_transform) which detrends the data.
    Simplest possible detrender: just divides the data into chunks and subtracts the mean in each chunk.
    """

    return rf_pipelines_c.make_simple_detrender(nt_chunk)
