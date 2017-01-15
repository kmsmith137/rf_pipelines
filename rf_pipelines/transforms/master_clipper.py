# TODO pass args (arr.dsample & pars), import .c, print proper
import numpy as np
import rf_pipelines

class master_clipper(rf_pipelines.py_wi_transform):
    """
    This transform utilizes helper functions to perform 
    an iterated clipping. The iteraion is constrained (hence
    optimized) by user-defined parameters.

    Constructor syntax:

      t = master_clipper(nt_chunk=1024, rms_cut=0., max_niter=1, cpp=False):
      
      'nt_chunk=1024' is the buffer size (in number of samples).
      
      'rms_cut=0.' is the rms threshold for the entire chunk.
       If the chunk rms is above this threshold, then all weights 
       are set to zero and further iterations are terminated.
    
      'max_niter=1' is the maximum number of iterations for each chunk.

      'cpp=False' is the flag which enables fast C++ algorithms.
    """

    def __init__(self, nt_chunk=1024, rms_cut=0., max_niter=1, cpp=False):
        
        assert nt_chunk > 0
        assert rms_cut >= 0., "rms threshold must be >= 0."
        assert max_niter >= 1

        self.nt_chunk = nt_chunk
        self.rms_cut = rms_cut
        self.max_niter = max_niter
        self.name = 'master_clipper(nt_chunk=%d, rms_cut=%f, max_niter=%d)' % (nt_chunk, rms_cut, max_niter)

    def set_stream(self, stream):
        
        self.nfreq = stream.nfreq
        
    def process_chunk(self, t0, t1, intensity, weights, pp_intensity, pp_weights):
        
        for ix in xrange(self.max_niter):
            
            (mean, rms) = rf_pipelines.weighted_mean_and_rms(intensity, weights, niter=6, sigma_clip=3) # TODO synch with cpp
            print "(mean, rms)=", (mean, rms)
            
            if rms > self.rms_cut:
                weights[:] = 0.
                print "master_clipper: weights are set to zero"
                break
            else:
                unmasked_percentage = np.count_nonzero(weights) / float(weights.size) * 100.
                print unmasked_percentage, "% not masked." # TODO break if constant
                if not cpp:
                    rf_pipelines.clip_fx(intensity, weights, thr=3)
                    rf_pipelines.clip_fx(intensity, weights, thr=3, axis=0)
                    rf_pipelines.clip_fx(intensity, weights, thr=3, axis=1)
                    rf_pipelines.filter_stdv(intensity, weights, thr=3, axis=1)
                else:
                    pass
