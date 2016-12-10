import numpy as np
import rf_pipelines

class master_clipper(rf_pipelines.py_wi_transform):
    """
    """
    def __init__(self, nt_chunk=1024, dsample_nt=None, rms_cut=0., max_niter=1):
        
        assert nt_chunk > 0
        assert rms_cut >= 0., "rms threshold must be >= 0."
        assert max_niter >= 1
        assert (dsample_nt is None or dsample_nt > 0), "Invalid downsampling number along the time axis!"

        self.nt_chunk = nt_chunk
        self.dsample_nt = dsample_nt
        self.rms_cut = rms_cut
        self.max_niter = max_niter
        self.name = 'master_clipper(nt_chunk=%d, dsample_nt=%d, rms_cut=%f, max_niter=%d)' % (nt_chunk, dsample_nt, rms_cut, max_niter)

    def set_stream(self, stream):
        
        self.nfreq = stream.nfreq
        
        if self.dsample_nt is None:
            self.dsample_nt = self.nt_chunk
        if self.nt_chunk % self.dsample_nt != 0:
            raise RuntimeError("clip_fx: current implementation requires 'dsample_nt' to be a divisor of 'nt_chunk'.")

    def process_chunk(self, t0, t1, intensity, weights, pp_intensity, pp_weights):
        
        for ix in xrange(self.max_niter):
            (mean, rms) = rf_pipelines.weighted_mean_and_rms(intensity, weights, niter=6, sigma_clip=3)

            if rms > self.rms_cut:
                weights[:] = 0.
            
            else:
                print "frac_unmasked=" # TODO
                rf_pipelines.clip_fx(intensity, weights, thr=3, dsample_nfreq=512, dsample_nt=self.dsample_nt/16)
                rf_pipelines.clip_fx(intensity, weights, thr=3, axis=0, dsample_nt=self.dsample_nt)
                rf_pipelines.clip_fx(intensity, weights, thr=3, axis=1, dsample_nt=self.dsample_nt)
                rf_pipelines.filter_stdv(intensity, weights, thr=3, axis=1, dsample_nt=self.dsample_nt/16)
