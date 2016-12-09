import numpy as np
import rf_pipelines

class chain_iterator(rf_pipelines.py_wi_transform):
    """
    chain clippers to enable global constraints
    """
    def __init__(self, nt_chunk, clipper_nt, nc):
        self.nt_chunk = nt_chunk
        self.clipper_nt = clipper_nt
        self.nc = nc

    def set_stream(self, stream):
        self.nfreq = stream.nfreq
    
    def process_chunk(self, t0, t1, intensity, weights, pp_intensity, pp_weights):
        
        # clipper chain
        rf.pipelines.clip_fx(intensity, weights, thr=3, dsample_nfreq=512, dsample_nt=clipper_nt/16)
        rf_pipelines.clip_fx(intensity, weights, thr=3, axis=0, dsample_nt=clipper_nt)
        rf_pipelines.clip_fx(intensity, weights, thr=3, axis=1, dsample_nt=clipper_nt)
        rf_pipelines.filter_stdv(intensity, weights, thr=3, axis=1, dsample_nt=clipper_nt/16)
