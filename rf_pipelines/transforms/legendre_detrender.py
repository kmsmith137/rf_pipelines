import numpy as np
import rf_pipelines

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
        
        # <<<<<<< Assuming self.axis = 0 >>>>>>>
        # --------------------------------------
        # A.1 precomputing P(x_i)_alpha; TODO to be moved to set_stream?
        # TODO have to create an inner func taking nfreq/nt_chunk as input
        x = np.linspace(-1, 1, self.nfreq) # TODO could use .tile to combine with p 
        coef = np.eye(self.deg+1)
        p = np.zeros([self.deg+1, self.nfreq])

        for d in xrange(self.deg+1):
            p[d,:] = np.polynomial.Legendre(coef[d,:])(x)
        
        # --------------------------------------
        # A.2 looping over time samples
        assert np.shape(weights)[self.axis] == np.shape(p)[1] # use self.axis convention
        # TODO change iterator `t' to sample/ sth else)
        for t in xrange(self.nt_chunk):
            
            # TODO weights[:,t] vs. weights[n,:] depending on self.axis
            M = np.dot(weights[:,t] * p, p.T) # (d+1) by (d+1)
            assert np.shape(M) == (self.deg+1, self.deg+1) # TODO check if necessary
            
            v = np.sum(weights[:,t] * intensity[:,t] * p, axis=1) # (d+1) by 1
            assert np.shape(v) == (self.deg+1,) # TODO check if necessary

            # --------------------------------------
            # A.3 conditional statements
            # TODO if poorly conditioned, then skip, else continue..
            
            # --------------------------------------
            # A.4
            c = np.dot(np.linalg.inv(M), v)
            assert np.size(c) == self.deg+1 # TODO check if necessary
            fit = np.polynomial.Legendre(c)(x)
            
            # --------------------------------------
            # A.5 TODO: what to do with fit? e.g., TODO tailor to `self.axis'
            intensity[:,t] -= fit # TODO combine with A.4 to save memory
