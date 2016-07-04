"""
PSRFITS stream.  This is a thin wrapper around a C++ implmentation.
See psrfits_stream.cpp, and python linkage in rf_pipelines_c.cpp.
"""

from rf_pipelines import rf_pipelines_c


def psrfits_stream(filename):
    """Returns a weighted intensity stream (wi_stream) from a single PSRFITS source file."""

    return rf_pipelines_c.make_psrfits_stream(filename)
