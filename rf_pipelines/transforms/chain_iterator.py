import numpy as np
import rf_pipelines

class chain_iterator(rf_pipelines.py_wi_transform):
    """
    """
    def __init__(self, nt_chunk):
        self.nt_chunk = nt_chunk
    def set_stream(self, stream):
        self.nfreq = stream.nfreq
    def process_chunk(self, t0, t1, intensity, weights, pp_intensity, pp_weights):
        pass
