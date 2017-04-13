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

    def __init__(self, v1_chunk=32, v2_chunk=192, w_cutoff=1.8, nt_chunk=1024):
        name = 'online_mask_filler(v1_chunk=%d, v2_chunk=%d, w_cutoff=%d, nt_chunk=%d)' % (v1_chunk, v2_chunk, w_cutoff, nt_chunk)
        rf_pipelines.py_wi_transform.__init__(self, name)
        self.nt_chunk = nt_chunk
        self.nt_prepad = 0
        self.nt_postpad = 0

    def set_stream(self, s):
        self.nfreq = s.nfreq

    def start_substream(self, isubstream):
        self.isubstream = isubstream
        
    def process_chunk(self, t0, t1, intensity, weights, pp_intensity, pp_weights):
        pass

    def end_substream(self):
        pass
