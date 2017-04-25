import rf_pipelines
import numpy as np

class online_mask_filler(rf_pipelines.py_wi_transform):
    """
    An online implementation of the mask_filler and variance_estimator that does not required a pre-generated file of variance 
    estimates! 

    Variance estimates are calculated independently for each frequency, with a final variance estimates being outputted for 
    each frequency on a timescale of ~10s (with the default values). First, variance estimates are made for each frequency 
    using v1_chunk number of samples. Then, once v2_chunk of those variance estimates have accumulated for each frequency, 
    a median is outputted and used as the final variance estimate. If no variance could be estimated because too much of 
    a channel is masked, a variance of 0 will be used (for now). The final estimate is used until enough samples are 
    accumulated for another final variance estimate to be produced. Hence, the variance estimate lags about 10s behind the 
    current variance. 
    
    The weights array is continually examined; if a weight is below w_cutoff, the weight is maximized and the corresponding 
    intensity is replaced with gaussian random noise with variance specified by the variance estimation part of the transform. 
    If if the weight is greater than or equal to w_cutoff, it is maximized and its intensity is left alone. If a variance 
    estiamte could not be generated for a particular frequency and time, the corresponding weights are set to 0 if they are 
    below w_cutoff. 

    tldr: Estimates the variance of each frequency independently on ~10s timescales and uses it to uniform the 
    weights/intensity arrays by filling in gaussian random noise with the calculated variance for low weight values. 

    Constructor arguments
    ---------------------
    v1_chunk - Number of samples used to generate a "v1 estimate" (the variance is computed across the sample)
    v2_chunk - Number of "v1 estimates" used to output a final variance (the median of the "v1 estimates")
    w_cutoff - Weight threshold below which the corresponding intensity will be replaced with gaussian random noise

    """

    def __init__(self, v1_chunk=32, var_weight=2e-3, var_clamp=3.3e-3, w_clamp=3.3e-3, w_cutoff=1.5, nt_chunk=1024):
        name = 'online_mask_filler(v1_chunk=%d, var_weight=%f, var_clamp=%f, w_clamp=%f, w_cutoff=%f)' % (v1_chunk, var_weight, var_clamp, w_clamp, w_cutoff)

        # Call base class constructor
        rf_pipelines.py_wi_transform.__init__(self, name)

        self.v1_chunk = v1_chunk
        self.nt_chunk = nt_chunk
        self.var_weight = var_weight
        self.var_clamp = var_clamp
        self.w_clamp = w_clamp
        self.w_cutoff = w_cutoff
        self.nt_prepad = 0
        self.nt_postpad = 0

        assert v1_chunk > 0
        assert nt_chunk > 0
        assert var_weight > 0
        assert var_clamp > 0
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
        self.tmp = np.zeros((self.nfreq))  # holds the temporary v1 estimates while they are being updated
        self.running_weights = np.zeros((self.nfreq))

    def process_chunk(self, t0, t1, intensity, weights, pp_intensity, pp_weights):
        # Loop over intensity/weights in chunks of size v1_chunk
        for ichunk in xrange(0, self.nt_chunk, self.v1_chunk):
            for frequency in xrange(self.nfreq):
                # Calculate the v1 for each frequency
                self.tmp[frequency] =  self._v1(intensity[frequency, ichunk:ichunk+self.v1_chunk], weights[frequency, ichunk:ichunk+self.v1_chunk])

            # Once v1s have been calculated for each frequency, update the weights and running variance
            non_zero_v1 = (self.tmp != 0)
            zero_v1 = np.logical_not(non_zero_v1)

            # For nonzero (successful) v1s, increase the weights (if possible) and update the running variance
            self.running_weights[non_zero_v1] = np.minimum(2.0, self.running_weights[non_zero_v1] + self.w_clamp)
            self.tmp[non_zero_v1] = np.minimum(self.tmp[non_zero_v1], self.running_var[non_zero_v1] + self.var_clamp)
            self.tmp[non_zero_v1] = np.maximum(self.tmp[non_zero_v1], self.running_var[non_zero_v1] + self.var_clamp)
            self.running_var[non_zero_v1]  = (1-self.var_weight) * self.running_var[non_zero_v1] + self.var_weight * self.tmp[non_zero_v1]

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
