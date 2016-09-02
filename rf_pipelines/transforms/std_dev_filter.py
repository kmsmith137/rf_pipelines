import numpy as np
import rf_pipelines

class std_dev_filter(rf_pipelines.py_wi_transform):
    """
    Masks freqs that have std.dev. deviating from the others by some sigma.
    Warning! There's no detrending to remove the bandpass signature right now,
    so this transform requires some more work.
    
    Constructor syntax:

        t = std_dev_filter(thr=6.,nt_chunk=1024)

        'thr=6.' is the sigma value to clip (computed wrt the array of all channel
            standard deviations). Note: clipping based on absolute value of deviations

        'nt_chunk=1024' is the buffer size
    """

    def __init__(self,thr=6.,nt_chunk=1024):
        self.thr = thr
        self.nt_chunk  = nt_chunk

    def set_stream(self,stream):
        self.nfreq = stream.nfreq

    def process_chunk(self, t0, t1, intensity, weights, pp_intensity, pp_weights):
        sd = np.ma.masked_where(weights==0.,intensity,copy=False).std(1)
        weights[abs(sd-sd.mean()) > self.thr*sd.std()] = 0.
