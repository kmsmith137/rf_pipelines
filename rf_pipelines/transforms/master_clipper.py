import numpy as np
import rf_pipelines

class master_clipper(rf_pipelines.py_wi_transform):
    """
    """
    def __init__(self, nt_chunk, dsample_nt):
        self.nt_chunk = nt_chunk
        self.dsample_nt = dsample_nt

        self.name = 'master_clipper(nt_chunk=%d, dsample_nt=%d)' % (nt_chunk, dsample_nt)

    def set_stream(self, stream):
        self.nfreq = stream.nfreq
    
    def process_chunk(self, t0, t1, intensity, weights, pp_intensity, pp_weights):
        
        # global constraints....
        (mean, rms) = rf_pipelines.weighted_mean_and_rms(intensity, weights, n_internal, thr)
        if rms > 0.002:
            pass
        else:
            # clipper chain
            rf_pipelines.clip_fx(intensity, weights, thr=3, dsample_nfreq=512, dsample_nt=self.dsample_nt/16)
            rf_pipelines.clip_fx(intensity, weights, thr=3, axis=0, dsample_nt=self.dsample_nt)
            rf_pipelines.clip_fx(intensity, weights, thr=3, axis=1, dsample_nt=self.dsample_nt)
            rf_pipelines.filter_stdv(intensity, weights, thr=3, axis=1, dsample_nt=self.dsample_nt/16)
