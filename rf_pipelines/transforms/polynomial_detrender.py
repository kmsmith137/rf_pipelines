import numpy as np

import rf_pipelines
from rf_pipelines import rf_pipelines_c


def polynomial_detrender(nt_chunk=1024, deg=0, axis=0, cpp=True, epsilon=0.01, test=False):
    """
    This transform removes a degree-d weighted-fit legendre 
    polynomial from the intensity along a specified axis. 

    + Currently based on the minimization of sum( w_i * (y_i-f(x_i)) )^2,
    where w_i, y_i, and f(x_i) are the weights, intensity and model 
    values of sample i, respectively.
    + Names are consistent with "chime_zerodm_notes".

    Constructor syntax:

      t = polynomial_detrender(nt_chunk=1024, deg=0, axis=0, cpp=True, epsilon=0.01, test=False)
      
      'nt_chunk=1024' is the buffer size.

      'deg=0' is the degree of fit.
      
      'axis=0' is the axis convention:
        0: along freq; constant time.
        1: along time; constant freq.

      'cpp=True' will use fast C++ transforms
      'cpp=False' will use reference python transforms

      'epsilon' controls the threshold for a polynomial fit being poorly conditioned
         (only meaningful if cpp=True)

      'test=False' enables a test mode (only meaningful if cpp=False)
    """

    if cpp:
        return rf_pipelines_c.make_polynomial_detrender(nt_chunk, axis, deg, epsilon)
    else:
        return polynomial_detrender_python(nt_chunk, deg, axis, test)


class polynomial_detrender_python(rf_pipelines.py_wi_transform):
    def __init__(self, nt_chunk=1024, deg=0, axis=0, test=False):
        
        assert (deg >= 0 and type(deg) == int), "degree must be an integer >= 0"
        assert axis in (0, 1), "axis must be 0 (along freq; constant time) or 1 (along time; constant freq)."
        
        self.nt_chunk = nt_chunk
        self.deg = deg
        self.axis = axis
        self.name = 'polynomial_detrender_python(nt_chunk=%d, deg=%d, axis=%d)' % (nt_chunk, deg, axis)
        self.nt_prepad = 0
        self.nt_postpad = 0
        self.test = test

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
        
        # The test mode replaces the weights and intensity with
        # simulated chunks (see __test below).
        if self.test:
            weights, intensity = self._polynomial_detrender_python__test(weights, intensity)

        # Checking whether the coefficients array matches
        # (in dimension) with the weights and intensity 
        # arrays along the selected axis.
        assert np.shape(weights)[self.axis] == np.shape(intensity)[self.axis] == self.N
        
        # Looping over the unselected axis, we subtract
        # the best fit (i.e., output of self.leg_fit()) 
        # from the intensity along the selected axis.
        for n in xrange(self.loop):
            if self.axis == 0:
                intensity[:,n] -= self.leg_fit(weights[:,n], intensity[:,n])
            else:
                intensity[n,:] -= self.leg_fit(weights[n,:], intensity[n,:])

        if self.test:
            # Print the degree of fit, the weighted mean 
            # and stdv computed over the entirety of each 
            # chunk. This shows the goodness of our detrender.
            mean = np.sum(weights*intensity)/np.sum(weights)
            stdv = np.sqrt(np.sum(weights*(intensity-mean)**2)/np.sum(weights))
            print (self.deg, mean, stdv), "= (deg, mean, stdv)"
            print "#####################################################"

    def leg_fit(self, w, i):
        """
        This method computes the coefficients of the best-fit leg polynomial for the array 'i' weighted by 'w'.
        The output is an array of the evaluated leg polynomial over the same domain as 'i'.
        """

        assert w.shape == i.shape == (self.N,)
        
        # Ill-conditioned Matrix M:
        if np.sum(w) < 20:
            return np.zeros(self.N)

        else:
            # The following section is explained in "chime_zerodm_notes".
            M = np.dot(w * self.P, self.P.T)
            assert np.shape(M) == (self.deg+1, self.deg+1)
            v = np.sum(w * i * self.P, axis=1)
            # This is a good place to call an inner 
            # method for catching an ill-conditioned M.
            c = np.dot(np.linalg.inv(M), v)
            assert np.size(c) == self.deg+1
            
            # Print the computed coeff of the best fit.
            # Warning: a very long print! 
            # But this is proven to be a very useful probe 
            # for catching an inconsistent output in test(_mode=1).
            if self.test:
                print c, "= Computed Coeff"

            # Use self.P to evaluate the poly (defined by c)
            # over the fit domain.
            return np.dot(c, self.P)

    def __test(self, weights, intensity, mask_level=0.9, test_mode=1):
        """This private method replaces the weights and intensity arrays with a new set of simulated arrays.
        
        + Search for the following words in the printed output on your screen: "Masked", "deg", "mean", "stdv", "Actual", "Computed".
        + Set 'mask_level' to a float between -1 and +1 (e.g., mask_level=0 masks out ~50% of the weights).
        + A visual test: In the main loop above, change '-=' to '=' so that the output of leg_fit() gets stored as the data.
        """

        # Let's create a weights array using a gaussian dist.
        weights[:] = np.random.normal(0, 1, weights.size).reshape(weights.shape)
        # Mask out weights less than 'mask_level'
        indy, indx = np.where(weights < mask_level)
        weights[indy,indx] = 0.

        print indx.size / float(weights.size) * 100, " % Masked"

        if test_mode == 0:
            # Fitting a noise level with randomly masked 
            # elements should result in a consistent chain 
            # of outputs. Test the code by running the transform 
            # (while self.test=True) with a few different 
            # values for 'self.deg'.
            intensity[:] = np.random.normal(0, 1, intensity.size).reshape(intensity.shape)
            
        if test_mode == 1:
            # Test by fitting a well-defined poly tiled in 2D.
            # The poly's coefficients are chosen randomly 
            # (centered at 0 with stdv=10)
            rc = np.random.normal(0, 10, self.deg+1)
            
            # Tricky! we have to tile along the unselected
            # axis because our polynomial is already in 
            # the direction of the selected axis.
            if self.axis == 0:
                tile_axis = 1
            else:
                tile_axis = 0

            # Print the actual coefficients, which change 
            # for each chunk -- the random generator gets 
            # called in the loop over the unselected axis.
            print rc, "= Actual Coeff"
            
            # Replace all intensity values by the 2D-tiled
            # simulated chunk.
            intensity[:] = rf_pipelines.tile_arr(\
                    np.dot(rc, self.P), tile_axis,\
                    self.nfreq, self.nt_chunk)
        
        # Simulated weights and intensity arrays.
        return weights, intensity
