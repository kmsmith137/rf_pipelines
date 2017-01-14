"""
Polynomial detrender.  This is a thin wrapper around a C++ implementation.
See polynomial_detrenders.cpp, and python linkage in rf_pipelines_c.cpp.
"""

from rf_pipelines import rf_pipelines_c

def polynomial_detrender_cpp(nt_chunk, axis=1, polydeg=0, epsilon=1.0e-2):
    """
    Returns a transform object (wi_transform) which detrends the intensity
    along time or frequency.

    Constructor syntax:

      t = polynomial_detrender_cpp(nt_chunk, axis=1, polydeg=0, epsilon=1.0e-2)
      
      'nt_chunk=1024' is the chunk size (in number of samples).

      'axis=1' is the axis convention:
        0: along freq; constant time.
        1: along time; constant freq.

      'polydeg=0' is the degree of fit.

      'epsilon=1.0e-2' is the threshold value for a poorly conditioned fit.
    """ 
    
    assert (polydeg >= 0 and type(polydeg) == int), "poly degree must be an integer >= 0"
    assert epsilon > 0., "choose a threshold value > 0 or leave it as default (1.0e-2)"

    return rf_pipelines_c.make_polynomial_detrender(nt_chunk, axis, polydeg, epsilon)
