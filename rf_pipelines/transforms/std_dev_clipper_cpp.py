"""
standard deviation clipper. This is a thin wrapper around a C++ implementation.
See std_dev_clippers.cpp, and python linkage in rf_pipelines_c.cpp.
"""

from rf_pipelines import rf_pipelines_c

def std_dev_clipper_cpp(nt_chunk=1024, axis=None, sigma=3, Df=1, Dt=1):
    """
    Returns a transform object (wi_transform) which clips an intensity array 
    by masking rows/columns whose standard deviation is an outlier. Results 
    are applied to the weights array (i.e., weights[clipped] = 0.).

    Constructor syntax:

      t = std_dev_clipper_cpp(nt_chunk=1024, axis=None, sigma=3, Df=1, Dt=1)
      
      'nt_chunk=1024' is the chunk size (in number of samples).

      'axis=None' is the axis convention:
        None: planar; freq and time. 
        0: along freq; constant time.
        1: along time; constant freq.
 
      'sigma=3' is the threshold (in sigmas from the mean) for clipping. 

      'Df=1' and 'Dt=1' are the frequency and time downsampling factors, respectively.
    """

    return rf_pipelines_c.make_std_dev_clipper(nt_chunk, axis, sigma, Df, Dt)
