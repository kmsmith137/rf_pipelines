import numpy as np
import rf_pipelines

class legendre_detrender(rf_pipelines.py_wi_transform):
    """
   This transform removes a degree-d weighted-fit legendre 
   polynomial from the intensity along a specified axis. 

   + Currently based on the minimization of sum( w_i * (y_i-f(x_i)) )^2,
   where w_i, y_i, and f(x_i) are the weights, intensity and model 
   values of sample i, respectively.
   + Names are consistent with "chime_zerodm_notes".

    Constructor syntax:

      t = legendre_detrender(deg=0, axis=0, nt_chunk=1024, test=False)
      
      'deg=0' is the degree of fit.
      
      'axis=0' is the axis convention:
        0: along freq; constant time
        1: along time; constant freq

      'nt_chunk=1024' is the buffer size.

      'test=False' enables a test mode.
    """

    def __init__(self, deg=0, axis=0, nt_chunk=1024, test=False):
        
        self.deg = deg
        self.axis = axis
        self.nt_chunk = nt_chunk
        self.nt_prepad = 0
        self.nt_postpad = 0
        self.test = test

        assert (self.deg >= 0 and type(self.deg) == int), \
            "degree must be an integer >= 0"
        assert (self.axis == 0 or self.axis == 1), \
            "axis must be 0 (along freq; constant time) or 1 (along time; constant freq)."

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
        
        if self.test:
            weights, intensity = self._legendre_detrender__test(weights, intensity)

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

        if self.test:
            mean = np.sum(weights*intensity)/np.sum(weights)
            rms = np.sqrt(np.sum(weights*(intensity-mean)**2)/np.sum(weights))
            print (self.deg, mean, rms), " = (deg, mean, rms)"

    def leg_fit(self, w, i):
        """This method computes the coefficients of the 
        best-fit leg polynomial for the array 'i' weighted by 'w'.
        The output is an array of the evaluated leg polynomial
        over the same domain as 'i'.
        """
        # Input should be 1d.
        assert w.ndim == i.ndim == 1
        assert w.size == i.size
        
        # Ill-conditioned Matrix M:
        # For now, let's just ignore totally-masked arrays. 
        # In the future, we shall implement a robust 
        # algorithm for catching poorly-conditioned matrices 
        # (see "chime_zerodm_notes").
        if np.sum(w) == 0.:
            return 0.
        else:
            # The following is explained in the notes.
            M = np.dot(w * self.P, self.P.T)
            assert np.shape(M) == (self.deg+1, self.deg+1)
            v = np.sum(w * i * self.P, axis=1)
            # This is a good place to call an inner 
            # method for catching an ill-conditioned M.
            c = np.dot(np.linalg.inv(M), v)
            assert np.size(c) == self.deg+1
            return np.dot(c, self.P)
    
    def __test(self, weights, intensity, mask_level=0.3, test_mode=1):
        """This private method replaces the weights and intensity
        arrays with a new set of pre-defined arrays (see below)"""

        # Let's create a weights array using a gaussian noise distribution.
        weights[:] = np.random.normal(0, 1, weights.size).reshape(weights.shape)
        # Mask out weights less than 'mask_level'
        indx, indy = np.where(weights < mask_level)
        weights[indx,indy] = 0.
        print indx.size / float(weights.size) * 100, " % masked"

        if test_mode == 0:
            # Fitting a noise level with randomly masked 
            # elements should result in a consistent chain 
            # of outputs with different degrees of fit.
            # Test the code by running the transform 
            # (while self.test=True) with a few different 
            # values for 'self.deg'.
            intensity[:] = np.random.normal(0, 1, intensity.size).reshape(intensity.shape)
            
        if test_mode == 1:
            # Test by fitting a well-defined poly tiled in 2d.
            # The poly's coefficients are chosen randomly 
            # (centered at 0 with stdv=10)
            c = np.random.normal(0, 10, self.deg+1)
            intensity[:] = rf_pipelines.tile_arr(\
                    np.dot(c, self.P), self.axis,\
                    self.nfreq, self.nt_chunk)

        return weights, intensity
