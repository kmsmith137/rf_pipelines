import numpy as np
import rf_pipelines

class noise_source_detrender(rf_pipelines.py_wi_transform):
    """
   Removes 0th-degree switched-noise-source signal by highpassing 
   weighted intensity at constant time values. 
   The switch turns on/off every 2**23 * 2.56 * 10**(-6) seconds.

    Constructor syntax:

      t = noise_source_detrender(nt_chunk=1024)

      'nt_chunk=1024' is the buffer size.
    """

    def __init__(self, nt_chunk=1024):
        
        self.nt_chunk = nt_chunk
        self.nt_prepad = 0
        self.nt_postpad = 0
        
    def set_stream(self, stream):
        
        self.nfreq = stream.nfreq
 
    def process_chunk(self, t0, t1, intensity, weights, pp_intensity, pp_weights):

        # Let's compute the numerator and denominator of
        # the weighted mean in the direction of frequencies,
        # i.e. constant in time.
        num = np.sum(intensity*weights, axis=0)
        den = np.sum(weights, axis=0)
        
        # Looping over the sum of weights, we subtract the 
        # well-defined (i.e., the sum is greater than zero)
        # weighted mean from the intensity (along the freq axis).
        t = 0
        for k in den:
            if k > 0.:
                intensity[:,t] -= (num[t]/k)
        t += 1
