import rf_pipelines
import numpy as np
import time

class variance_estimator(rf_pipelines.py_wi_transform):
    """
    This is a pseudo-transform (meaning that it does not actually modify its input). 
    
    This estimates the variance of the intensity/weights arrays passed to it. First, 
    it computes the variance across v1_chunk number of timesamples. Then, it takes the 
    median of v2_chunk number of v1 estimates. 

    Thus, the final variance array produced is of dimensions
        (nfreq, total number of time samples / v1_chunk / v1_chunk)

    The variance array is written out as a .npy file in the directory of the test script.
    """

    def __init__(self, v1_chunk=128, v2_chunk=64, nt_chunk=1024):
        name = "variance_estimator('v1_chunk=%d, v2_chunk=%d, nt_chunk=%d')" % (v1_chunk, v2_chunk, nt_chunk)

        # Call base class constructor
        rf_pipelines.py_wi_transform.__init__(self, name)

        assert nt_chunk % v1_chunk == 0, \
            'For now, nt_chunk(=%d) must be a multiple of v1_chunk(=%d)' % (nt_chunk, v1_chunk)
                self.v1_chunk = v1_chunk
        
        self.v2_chunk = v2_chunk
        self.nt_chunk = nt_chunk
        self.nt_prepad = 0  
        self.nt_postpad = 0 
        

    def set_stream(self, s):
        # As explained in the rf_pipelines.py_wi_transform docstring, this function is
        # called once per pipeline run, when the stream (the 's' argument) has been specified.
        self.nfreq = s.nfreq


    def start_substream(self, isubstream, t0):
        # Called once per substream (a stream can be split into multiple substreams).
        self.v1 = np.zeros((self.nfreq, self.v1_chunk))
        self.iv1 = np.zeros((self.nfreq), dtype=np.int32) # keeps track of which positon we are adding v1 to 
        self.v2 = []


    def process_chunk(self, t0, t1, intensity, weights, pp_intensity, pp_weights):
        # This is the main computational routine defining the transform, which is called
        # once per incoming "block" of data.  For documentation on the interface see
        # rf_pipelines.py_wi_transform docstring.
        for frequency in xrange(self.nfreq):
            for i in xrange(0, len(intensity), self.v1_chunk):
                # Compute v1 for each frequency
                a = self._v1(intensity[frequency, i : i+self.v1_chunk], weights[frequency, i : i+self.v1_chunk])
                self.v1[frequency, self.iv1[frequency]] = a
                self.iv1[frequency] += 1
                if self.iv1[frequency] == self.v2_chunk:
                    # Empty v1[frequency] and compute median for v2 (nan if too many v1 values are nan)
                    if np.isnan(self.v1[frequency]).sum() > self.v2_chunk * 0.75:
                        self.v2.append(np.nan)
                    else:
                        self.v2.append(np.median(self.v1[frequency]))
                    self.v1[frequency] = []
                    self.iv1[frequency] = 0


    def end_substream(self):
        # Reshape self.v2
        out = np.array(self.v2).reshape((self.nfreq, -1))
        print out.shape
        print out

        # Interpolate nans 
        for frequency in xrange(self.nfreq):
            a = out[frequency]
            nans, x = self._nan_helper(a)
            a[nans]= np.interp(x(nans), x(~nans), a[~nans])
        
        # Write data to script directory
        np.save('var-%s.npy' % (time.strftime('%y-%m-%d-%X')), out)


    def _v1(self, i, w):
        """Calculate weighted variance"""
        w_sum = w.sum()
        if np.count_nonzero(w == 0) > self.v1_chunk * 0.75:
            return np.nan
        average = np.average(i, weights=w)
        variance = np.average((i-average)**2, weights=w) 
        return variance
        

    def _nan_helper(self, a):
        """Input 1D array"""
        return np.isnan(a), lambda z: z.nonzero()[0]
