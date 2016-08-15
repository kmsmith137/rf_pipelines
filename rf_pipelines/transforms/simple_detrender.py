"""
Simple detrender.  This is a thin wrapper around a C++ implementation.
See simple_detrender.cpp, and python linkage in rf_pipelines_c.cpp.
"""

from rf_pipelines import rf_pipelines_c


def simple_detrender(nt_detrend):
    """
    Returns a transform object (wi_transform) which detrends the data.

    This the simplest possible detrending algorithm.  We really need something better here! 
    It just divides the data into chunk, and subtracts the time-average of the data for every 
    (chunk, frequency_channel) pair.

    The 'nt_detrend' constructor arg determiens the chunk size (in number of samples).
    """

    return rf_pipelines_c.make_simple_detrender(nt_detrend)
