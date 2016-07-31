import numpy as np
import rf_pipelines
import numpy.polynomial.legendre as npl

class legendre_detrender(rf_pipelines.py_wi_transform):
    """
   This transform removes a degree-d weighted-fit legendre 
   polynomial from the intensity along a specified axis. 

   Currently based on the minimization of sum( w_i * (y_i-f(x_i)) )^2,
   where w_i, y_i, and f(x_i) are the weights, intensity and model 
   values of sample i, respectively.

    Constructor syntax:

      t = legendre_detrender(deg=0, axis=0, nt_chunk=1024)
      
      'deg=0' is the degree of fit.
      
      'axis=0' is the axis convention:
        0: along freq; constant time
        1: along time; constant freq

      'nt_chunk=1024' is the buffer size.
    """

    def __init__(self, deg=0, axis=0, nt_chunk=1024):
        
        self.deg = deg
        self.axis = axis
        self.nt_chunk = nt_chunk
        self.nt_prepad = 0
        self.nt_postpad = 0

        assert (self.deg >= 0 and type(self.deg) == int), \
            'degree must be an integer >= 0'
        assert (self.axis == 0 or self.axis == 1), \
            'axis must be 0 (along freq; constant time) or 1 (along time; constant freq).'

    def set_stream(self, stream):
        
        self.nfreq = stream.nfreq

    def process_chunk(self, t0, t1, intensity, weights, pp_intensity, pp_weights):
