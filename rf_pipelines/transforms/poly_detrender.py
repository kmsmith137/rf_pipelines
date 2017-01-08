"""
Poly detrender.  This is a thin wrapper around a C++ implementation.
See polynomial_detrenders.cpp, and python linkage in rf_pipelines_c.cpp.
"""

from rf_pipelines import rf_pipelines_c

def poly_detrender(nt_detrend, axis=1, deg=0, epsilon=1.0e-2):
    """
    Returns a transform object (wi_transform) which detrends the intensity
    along time or frequency.

    Constructor syntax:

      t = poly_detrender(nt_detrend, axis=1, deg=0, epsilon=1.0e-2)
      
      'nt_detrend=1024' is the chunk size (in number of samples).

      'axis=0' is the axis convention:
        0: along freq; constant time.
        1: along time; constant freq.

      'deg=0' is the degree of fit.

      'epsilon=1.0e-2' is the threshold value for a poorly conditioned fit.

    """ 
    assert (deg >= 0 and type(deg) == int), "degree must be an integer >= 0"
    assert epsilon > 0., "choose a threshold value > 0 or leave it as default (1.0e-2)"

    if axis == 1:
        return rf_pipelines_c.make_polynomial_detrender_time_axis(nt_detrend, deg, epsilon)
    elif axis == 0:
        return rf_pipelines_c.make_polynomial_detrender_freq_axis(nt_detrend, deg, epsilon)
    else:
        raise RuntimeError("axis must be 0 (along freq) or 1 (along time).")
