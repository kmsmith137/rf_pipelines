import numpy as np
import rf_pipelines

class thermal_noise_weight(rf_pipelines.py_wi_transform):
    """ 
    This transform is a rewrite of Kiyo's ch_L1Mock code (preprocess.py)
    For thermal noise dominated data, this should yield meaningful S/Ns
    Drawback - I don't think we can use this in tandem with detrending.
    """

    def __init__(self,nt_chunk=512):
        rf_pipelines.py_wi_transform.__init__(self)
        self.nt_chunk = nt_chunk

    def set_stream(self,stream):
        self.nfreq   = stream.nfreq
        self.delta_f = (stream.freq_hi_MHz-stream.freq_lo_MHz)*1e6/self.nfreq
        self.delta_t = stream.dt_sample
    
    def process_chunk(self, t0, t1, intensity, weights, pp_intensity, pp_weights):
        num = np.sum(intensity*weights,-1)
        den = np.sum(weights,-1)
        bad_freq = np.logical_or(den < 0.001 * np.mean(den), num == 0) 
        den[bad_freq] = 1
        time_mean = num / den
        time_mean[bad_freq] = 0
        std_thermal = time_mean/np.sqrt(2*self.delta_t*self.delta_f)
        std_thermal[bad_freq] = 1
        weights /= std_thermal[:,None]
        weights[bad_freq] = 0
