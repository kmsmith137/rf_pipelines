"""
Intensity clipper.  This is a thin wrapper around a C++ implementation.
See FIXME intensity_clippers.cpp, and python linkage in rf_pipelines_c.cpp.
"""

from rf_pipelines import rf_pipelines_c

def intensity_clipper(axis=None, Df=1, Dt=1, nt_chunk=1024, sigma=3, niter=1, iter_sigma=3):

    """
    Returns a transform object (wi_transform) which uses the intensity
    array for clipping extreme values. Results are applied to the weights 
    array (i.e., weights[clipped] = 0.).

    Constructor syntax:

      t = intensity_clipper(axis=None, Df=1, Dt=1, nt_chunk=1024, sigma=3, niter=1, iter_sigma=3)
      
      'axis=None' is the axis convention:
        None: planar; freq and time. 
        0: along freq; constant time.
        1: along time; constant freq.
 
      'Df=1' and 'Dt=1' are the frequency and time downsampling factors, respectively.

      'nt_chunk=1024' is the chunk size (in number of samples).

      'sigma=3' is the threshold (in sigmas from the mean) for clipping. 
        Note: The weights are used when calculating both the mean and rms intensity.

      'niter=1' is the number of internal iterations before the final clipping.

      'iter_sigma=3' is the threshold value for clipping in internal iterations.
        Note: This value need not be the same as 'sigma'

    """
    assert (axis == None or axis == 0 or axis == 1), "axis must be None (planar; freq and time), 0 (along freq; constant time), or 1 (along time; constant freq)."
    assert (sigma >= 1. and iter_sigma >= 1.), "threshold values must be >= 1."
    assert TODO
    assert nt_chunk > 0
    
    if axis == 1:
        return TODO
    elif axis == 0:
        return TODO
    else:
        return rf_pipelines_c.make_intensity_clipper2d(Df, Dt, nt_chunk, sigma, niter, iter_sigma)
