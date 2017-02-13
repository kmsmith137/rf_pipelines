import sys
import rf_pipelines
import numpy as np

# HEY FORGETFUL MAYA: make sure you add interpolation at some point! 


np.set_printoptions(threshold=np.nan)


class variance_estimator(rf_pipelines.py_wi_transform):
    """
    This is a pseudo-transform (meaning that it does not actually modify its input). 

    Current strategy: calculate variance of v1_chunk number of samples return the median of 
    v2_chunk number of v1 values. 
    
    Maybe introduce a frequency downsample parameter later on?
    """

    def __init__(self, v1_chunk=64, v2_chunk = 32, nt_chunk=1024):
        name = "variance_estimator('v1_chunk=%d, nt_chunk=%d')" % (v1_chunk, nt_chunk)

        # Call base class constructor
        rf_pipelines.py_wi_transform.__init__(self, name)

        assert nt_chunk % v1_chunk == 0, \
            'For now, nt_chunk(=%d) must be a multiple of v1_chunk(=%d)' % (nt_chunk, v1_chunk)
        
        self.v1_chunk = v1_chunk
        self.v2_chunk = v2_chunk
        self.nt_chunk = nt_chunk
        self.nt_prepad = 0   # might want to do something with this later
        self.nt_postpad = 0  # but for now, assume nt_chunk is all we are looking at (must be multiple of v1_samples
        
 
    def set_stream(self, s):
        # As explained in the rf_pipelines.py_wi_transform docstring, this function is
        # called once per pipeline run, when the stream (the 's' argument) has been specified.
        self.nfreq = s.nfreq

        # if s.nfreq % self.img_nfreq != 0:
        #     raise RuntimeError("variance_estimator: current implementation requires 'img_nfreq' to be a divisor of stream nfreq (idk why... relic from plotter transform)")


    def start_substream(self, isubstream, t0):
        # Called once per substream (a stream can be split into multiple substreams).
        self.v1 = np.zeros((self.nfreq, self.v1_chunk))
        self.iv1 = np.zeros((self.nfreq), dtype=np.int32) # keeps track of which positon we are adding v1 to 
        self.v2 = []

    def process_chunk(self, t0, t1, intensity, weights, pp_intensity, pp_weights):
        # This is the main computational routine defining the transform, which is called
        # once per incoming "block" of data.  For documentation on the interface see
        # rf_pipelines.py_wi_transform docstring.

        # For testing - we are only viewing 4 frequencies
        for frequency in xrange(self.nfreq-1000-20):
            if frequency == 0:
                print intensity[0]
                print weights[0]
            for i in xrange(0, len(intensity), self.v1_chunk):
                a = self._v1(intensity[frequency, i : i+self.v1_chunk], weights[frequency, i : i+self.v1_chunk])
                self.v1[frequency, self.iv1[frequency]] = a
                self.iv1[frequency] += 1
                if self.iv1[frequency] == self.v2_chunk:
                    if np.isnan(self.v1[frequency]).sum() > self.v2_chunk * 0.75:
                        print "** NAN **"
                        self.v2.append(np.nan)
                    else:
                        print 'med for', frequency, np.median(self.v1[frequency])
                        self.v2.append(np.median(self.v1[frequency]))
                    self.v1[frequency] = []
                    self.iv1[frequency] = 0

    def end_substream(self):
        # Reshape self.v2
        out = np.array(self.v2).reshape((self.nfreq-1000-20, -1))
        print '-' * 80
        print out.shape
        print out
        
        np.save('foo.npy', out)


    def _v1(self, i, w):
        w_sum = w.sum()
        if np.count_nonzero(w == 0) > self.v1_chunk * 0.75:
            print "*** NAN ***"
            return np.nan
        average = np.average(i, weights=w)
        variance = np.average((i-average)**2, weights=w)  # Fast and numerically precise
        return variance
        
        
