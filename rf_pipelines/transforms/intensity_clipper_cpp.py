"""
Intensity clipper.  This is a thin wrapper around a C++ implementation.
See intensity_clippers.cpp, and python linkage in rf_pipelines_c.cpp.
"""

import numpy as np
from rf_pipelines import rf_pipelines_c

def intensity_clipper_cpp(nt_chunk=1024, axis=None, sigma=3, niter=1, iter_sigma=3, Df=1, Dt=1, two_pass=False):
    """
    Returns a transform object (wi_transform) which uses the intensity
    array for clipping extreme values. Results are applied to the weights 
    array (i.e., weights[clipped] = 0.).

    Constructor syntax:

      t = intensity_clipper_cpp(nt_chunk=1024, axis=None, sigma=3, niter=1, iter_sigma=3, Df=1, Dt=1)
      
      'nt_chunk=1024' is the chunk size (in number of samples).

      'axis=None' is the axis convention:
        None: planar; freq and time. 
        0: along freq; constant time.
        1: along time; constant freq.
 
      'sigma=3' is the threshold (in sigmas from the mean) for clipping. 
       Note: The weights are used when calculating both the mean and rms intensity.

      'niter=1' is the number of internal iterations before the final clipping.

      'iter_sigma=3' is the threshold value for clipping in internal iterations.
       Note: This value need not be the same as 'sigma'

      'Df=1' and 'Dt=1' are the frequency and time downsampling factors, respectively.

      If two_pass=True, then a more numerically stable but somewhat slower algorithm will be used.
      (In CHIME, this should only be needed in a specific case: when analyzing incoherent-beam data, 
       with axis=1, and before the first detrender in the transform chain.
    """

    assert ((np.log2(Df) % 2) in (0., 1.)) and ((np.log2(Dt) % 2) in (0., 1.)), "Downsampling factors must be powers of 2"
    
    return rf_pipelines_c.make_intensity_clipper(nt_chunk, axis, sigma, niter, iter_sigma, Df, Dt, two_pass)
