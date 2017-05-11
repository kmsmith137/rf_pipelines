import rf_pipelines
import rf_pipelines.rf_pipelines_c

import sys
import numpy as np


class py_online_mask_filler(rf_pipelines.py_wi_transform):
    """For documentation of the mask_filler, see the 'rf_pipelines.online_mask_filler' docstring."""
    def __init__(self, v1_chunk=32, var_weight=2e-3, var_clamp_add=3.3e-3, var_clamp_mult=3.3e-3, w_clamp=3.3e-3, w_cutoff=0.5, nt_chunk=1024):
        name = 'online_mask_filler(v1_chunk=%d, var_weight=%.4f, var_clamp_add=%.4f, var_clamp_mult=%.4f, w_clamp=%.4f, w_cutoff=%.2f)' \
               % (v1_chunk, var_weight, var_clamp_add, var_clamp_mult, w_clamp, w_cutoff)

        # Call base class constructor
        rf_pipelines.py_wi_transform.__init__(self, name)

        self.v1_chunk = v1_chunk
        self.nt_chunk = nt_chunk
        self.var_weight = var_weight
        self.var_clamp_add = var_clamp_add
        self.var_clamp_mult = var_clamp_mult
        self.w_clamp = w_clamp
        self.w_cutoff = w_cutoff
        self.nt_prepad = 0
        self.nt_postpad = 0

        assert v1_chunk > 0
        assert nt_chunk > 0
        assert var_weight > 0
        assert var_clamp_add >= 0
        assert var_clamp_mult >= 0
        assert w_clamp > 0
        assert w_cutoff > 0
        assert nt_chunk % v1_chunk == 0


    def set_stream(self, s):
        # As explained in the rf_pipelines.py_wi_transform docstring, this function is
        # called once per pipeline run, when the stream (the 's' argument) has been specified.
        self.nfreq = s.nfreq


    def start_substream(self, isubstream, t0):
        # Called once per substream (a stream can be split into multiple substreams).
        self.running_var = np.zeros((self.nfreq))  # holds the current variance estimate for each frequency
        self.v1_tmp = np.zeros((self.nfreq))  # holds the temporary v1 estimates while they are being updated
        self.running_weights = np.zeros((self.nfreq))
        self.var_init = np.ones((self.nfreq), dtype=bool)

    def process_chunk(self, t0, t1, intensity, weights, pp_intensity, pp_weights):
        # Loop over intensity/weights in chunks of size v1_chunk
        for ichunk in xrange(0, self.nt_chunk, self.v1_chunk):
            for frequency in xrange(self.nfreq):
                # Calculate the v1 for each frequency
                self.v1_tmp[frequency] =  self._v1(intensity[frequency, ichunk:ichunk+self.v1_chunk], weights[frequency, ichunk:ichunk+self.v1_chunk])
                # Check whether the channel has been initialized. If not and a v1 estimate was successfully produced, use that
                if not self.var_init[frequency] and self.v1_tmp[frequency] != 0:
                    self.var_init[frequency] = True
                    self.running_var = self.v1_tmp[frequency]

            # Once v1s have been calculated for each frequency, update the weights and running variance
            non_zero_v1 = np.logical_and(self.v1_tmp != 0, self.var_init)
            zero_v1 = np.logical_and(np.logical_not(non_zero_v1), self.var_init)

            # For nonzero (successful) v1s, increase the weights (if possible) and update the running variance
            self.running_weights[non_zero_v1] = np.minimum(2.0, self.running_weights[non_zero_v1] + self.w_clamp)
            self.v1_tmp[non_zero_v1] = np.minimum(self.v1_tmp[non_zero_v1], self.running_var[non_zero_v1] + self.var_clamp_add + self.running_var[non_zero_v1] * self.var_clamp_mult)
            self.v1_tmp[non_zero_v1] = np.maximum(self.v1_tmp[non_zero_v1], self.running_var[non_zero_v1] - self.var_clamp_add - self.running_var[non_zero_v1] * self.var_clamp_mult)
            self.running_var[non_zero_v1]  = (1-self.var_weight) * self.running_var[non_zero_v1] + self.var_weight * self.v1_tmp[non_zero_v1]

            # For unsuccessful v1s, decrease the weights (if possible) and do not modify the running variance 
            self.running_weights[zero_v1] = np.maximum(0, self.running_weights[zero_v1] - self.w_clamp)
            
            # Mask fill!
            intensity_valid = (weights[:, ichunk:ichunk+self.v1_chunk] > self.w_cutoff)
            rand_intensity = np.random.standard_normal(size=intensity[:, ichunk:ichunk+self.v1_chunk].shape)
            for (ifreq,v) in enumerate(self.running_var):
                if v > 0.0:
                    rand_intensity[ifreq, :] *= v**0.5
            intensity[:, ichunk:ichunk+self.v1_chunk] = np.where(intensity_valid, intensity[:, ichunk:ichunk+self.v1_chunk], rand_intensity)
            weights[:, ichunk:ichunk+self.v1_chunk] = np.repeat(self.running_weights, self.v1_chunk).reshape(self.nfreq, self.v1_chunk)


    def end_substream(self):
        pass


    def _v1(self, i, w):
        """Calculate weighted variance"""
        if np.count_nonzero(w) < self.v1_chunk * 0.25:
            return 0
        return np.average(i**2, weights=w) 


####################################################################################################


# The externally-visible online_mask_filler() is now a function which can return
# either a C++ transform object or a python transform object, depending on the
# value of its 'cpp' argument.

def online_mask_filler(v1_chunk=32, var_weight=2e-3, var_clamp_add=3.3e-3, var_clamp_mult=3.3e-3, w_clamp=3.3e-3, w_cutoff=0.5, nt_chunk=1024, cpp=False):
    """
    Variance estimates are calculated independently for each frequency. A variance value is computed for each v1_chunk 
    samples. Then, an "exponential average" variance is computed, in which the new variance estimate for a frequency is 
    given a weight of var_weight and the previous variance estimate is given a weight of 1-var_weight. These are summed 
    to arrive at the updated variance estimate. Note that the variance estimate for a particular frequency cannot vary 
    by more than var_clamp in a particular chunk to account for outliers. 
    
    If a new variance estimate is successful, the corresponding weight will be increased (if possible) by w_clamp. If too
    much of the channel is masked to produce a variance estimate, the corresponding weight will be decreased (if possible)
    by w_clamp (stored in a temporary weights array).

    Finally, the intensity values are modified. If the corresponding weight is below w_cutoff, the intensity value will
    be filled with gaussian random noise with a variance given by the calculated variance estimate. 

    tldr: Estimates the variance of each frequency independently on ~10s timescales and uses it to uniform the 
    weights/intensity arrays by filling in gaussian random noise with the calculated variance for low weight values. 

    Constructor arguments
    ---------------------
    v1_chunk - Number of samples used to generate a variance estimate (the variance is computed across the sample)
    var_weight - The weighting given to each new variance estimate in the exponential average
    var_clamp_add - The amount the variance can vary by in each v1_chunk per frequency
    var_clamp_mult - The multiplicative amoun the variance can vary by in each v1_chunk per frequency
    w_clamp - The amount the weight can vary by in each v1_chunk per frequency
    w_cutoff - Weight threshold below which the corresponding intensity will be replaced with gaussian random noise
    cpp - If True, then the fast "production" C++ implementation of the transform will be used.
          If False, then the python reference implementation will be used.
          Currently, the C++ implementation is a placeholder that doesn't do anything!
    """

    if cpp:
        # Return C++ transform.
        print >>sys.stderr, "Warning: the C++ online_mask_filler is currently a placeholder that doesn't do anything!"
        return rf_pipelines.rf_pipelines_c.make_online_mask_filler(v1_chunk, v2_chunk, w_cutoff, nt_chunk)

    # Return instance of the py_online_mask_filler class above.
    return py_online_mask_filler(v1_chunk, var_weight, var_clamp_add, var_clamp_mult, w_clamp, w_cutoff, nt_chunk)
