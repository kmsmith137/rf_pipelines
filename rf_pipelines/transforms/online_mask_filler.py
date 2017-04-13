import rf_pipelines
import numpy as np

class online_mask_filler(rf_pipelines.py_wi_transform):
    """
    Some docstring
    """
    def __init__(self, v1_chunk, v2_chunk, w_cutoff, nt_chunk):
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
