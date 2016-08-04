import numpy as np
import rf_pipelines

class legendre_detrender(rf_pipelines.py_wi_transform):
    """
   This transform removes a degree-d weighted-fit legendre 
   polynomial from the intensity along a specified axis. 

   Currently based on the minimization of sum( w_i * (y_i-f(x_i)) )^2,
   where w_i, y_i, and f(x_i) are the weights, intensity and model 
   values of sample i, respectively.

   Names are based on "chime_zerodm_notes"

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
        
        # The following statement prepares the code for
        # looping over the unselected axis. self.N is 
        # the number of elements along the selected axis.
        if self.axis == 0:
            (self.N, self.loop) = (self.nfreq, self.nt_chunk)
        else:
            (self.N, self.loop) = (self.nt_chunk, self.nfreq)
        
        # Here we initialize a coefficients array (self.P; 
        # evaluated over self.N) of legendre polynomials 
        # from degree 0 to self.deg.
        self.x = np.linspace(-1, 1, self.N)
        self.P = np.zeros([self.deg+1, self.N])
        for d in xrange(self.deg+1):
            self.P[d,:] = np.polynomial.Legendre(np.eye(self.deg+1)[d,:])(self.x) 
    
    def process_chunk(self, t0, t1, intensity, weights, pp_intensity, pp_weights):
        
        # Checking whether the coefficients array matches
        # (in dimension) with the weights and intensity 
        # arrays along the selected axis.
        assert np.shape(weights)[self.axis] == np.shape(intensity)[self.axis] ==\
                np.shape(self.P)[1]
        
        # Looping over the unselected axis, we subtract
        # the best fit (i.e., output of self.leg_fit()) 
        # from the intensity along the selected axis.
        for n in xrange(self.loop):
            if self.axis == 0:
                intensity[:,n] -= self.leg_fit(weights[:,n], intensity[:,n])
            else:
                intensity[n,:] -= self.leg_fit(weights[n,:], intensity[n,:])

    def leg_fit(self, w, i):     
        #%%%%% assert: input; skip if w=0, etc.
        #%%%%% case: chunks > nt_chunk
        #%%%%% 0.5*nt_chunk; impl self.wstdv; self.step.
        if np.sum(w) == 0.: #%%%%% 0.5?
            return 0.
        else:
            M = np.dot(w * self.P, self.P.T)
            assert np.shape(M) == (self.deg+1, self.deg+1)
            v = np.sum(w * i * self.P, axis=1)
            #%%%%% A.3 >>> CONDITIONAL <<<
            c = np.dot(np.linalg.inv(M), v)
            assert np.size(c) == self.deg+1
            return np.polynomial.Legendre(c)(self.x)
