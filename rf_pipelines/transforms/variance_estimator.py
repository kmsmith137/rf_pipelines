import sys
import rf_pipelines
import numpy as np
np.set_printoptions(threshold=np.nan)
# import matplotlib.pyplot as plt

# First implementation - plot v_1 over time interval on long run and plot! 
# Good way to test chunk assumption! 


class variance_estimator(rf_pipelines.py_wi_transform):
    """
    This is a pseudo-transform (meaning that it does not actually modify its input). 

    For now, compute v_1 with short timesamples and plot!

    Maybe introduce a frequency downsample parameter later on?
    """

    def __init__(self, v1_chunk=64, nt_chunk=1024):
        name = "variance_estimator('v1_chunk=%d, nt_chunk=%d')" % (v1_chunk, nt_chunk)

        # Call base class constructor
        rf_pipelines.py_wi_transform.__init__(self, name)

        assert nt_chunk % v1_chunk == 0, \
            'For now, nt_chunk(=%d) must be a multiple of v1_chunk(=%d)' % (nt_chunk, v1_chunk)
        
        self.v1_chunk = v1_chunk
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
        self.v1 = []


    def process_chunk(self, t0, t1, intensity, weights, pp_intensity, pp_weights):
        # This is the main computational routine defining the transform, which is called
        # once per incoming "block" of data.  For documentation on the interface see
        # rf_pipelines.py_wi_transform docstring.

        for frequency in xrange(self.nfreq-1000-20):
            for i in xrange(0, len(intensity), self.v1_chunk):
                a = self._v1(intensity[frequency, i : i+self.v1_chunk], weights[frequency, i : i+self.v1_chunk])
                self.v1.append(a)
                print '*' * 10
                if a > 10:
                    print a
                    print intensity[frequency, i : i+self.v1_chunk]
                    print weights[frequency, i : i+self.v1_chunk]



    def end_substream(self):
        # Reshape self.v1
        out = np.array(self.v1).reshape((self.nfreq-1000-20, -1))
        print '-' * 80
        print out.shape
        print out
        
        np.save('foo.npy', out)


    def _v1(self, i, w):
        w_sum = w.sum()
        if w_sum == 0: 
            return None
        average = np.average(i, weights=w)
        variance = np.average((i-average)**2, weights=w)  # Fast and numerically precise
        return variance
        
        
