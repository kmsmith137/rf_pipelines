# TODO import .c, print proper, pars as dict (if too many)
import numpy as np
from types import DictType
import rf_pipelines

class master_transform(rf_pipelines.py_wi_transform):
    """
    This transform utilizes HELPER FUNCTIONS to perform 
    an iterated chian of transforms, which include 
    detrending and clipping algorithms. The iteraion is 
    controlled (hence optimized) by user-defined parameters.

    Constructor syntax:

      t = master_transform(nt_chunk=1024, fdict=None, rms_cut=0., max_niter=1):
      
      'nt_chunk=1024' is the buffer size (in number of samples).

      'fdict=None' must be passed as a dictionary which contains the following lists of transforms:
       
       fdict = {'py':[], 'imitate_cpp':[], 'cpp':[]}
       
       where the key:value pairs are:
        
        'py': [python-based helper functions]
        'imitate_cpp': [python-based cpp-imitated helper functions]
        'cpp': [cpp-based helper functions]
       
       e.g., fdict = {'py' : [ "clip_fx(...)", "filter_stdv(...)" ], 
                      'imitate_cpp' : [ "filter_stdv(..., imitate_cpp=True)" ], 
                      'cpp' : [ "rf_pipelines_c.apply_intensity_clipper(...)" ]
                     }
      
      'rms_cut=0.' is the rms threshold for the entire chunk.
       If the chunk rms is above this threshold, then all weights 
       are set to zero and further iterations are terminated.
    
      'max_niter=1' is the maximum number of iterations for each chunk.
    """

    def __init__(self, nt_chunk=1024, fdict=None, rms_cut=0., max_niter=1):
        
        assert nt_chunk > 0
        assert (type(fdict) is DictType) and ({for i in fdict.keys()} == {'py', 'imitate_cpp', 'cpp'}),
            "master_transform: 'fdict' must be a dictionary with the following format:\n
            fdict = {'py':[], 'imitate_cpp':[], 'cpp':[]}"
        assert rms_cut >= 0., "master_transform: rms threshold must be >= 0."
        assert max_niter >= 1

        self.nt_chunk = nt_chunk
        self.fdict = fdict
        self.rms_cut = rms_cut
        self.max_niter = max_niter
        self.name = 'master_transform(nt_chunk=%d, rms_cut=%f, max_niter=%d)' % (nt_chunk, rms_cut, max_niter) # TODO compressed fdict (maybe)

    def set_stream(self, stream):
        
        self.nfreq = stream.nfreq
        
    def process_chunk(self, t0, t1, intensity, weights, pp_intensity, pp_weights):
        
        for ix in xrange(self.max_niter):
            
            (mean, rms) = rf_pipelines.weighted_mean_and_rms(intensity, weights, niter=6, sigma_clip=3) # TODO synch with updated cpp
            print "(mean, rms)=", (mean, rms)
            
            if rms > self.rms_cut:
                weights[:] = 0.
                print "master_clipper: weights are set to zero"
                break
            else:
                unmasked_percentage = np.count_nonzero(weights) / float(weights.size) * 100.
                print unmasked_percentage, "% not masked." # TODO break if constant
                
                # TODO fdict.values management

                # python-based helper functions
                if 'py' in fdict:
                    for pyf in fdict['py']:
                        exec(pyf)

                # python-based cpp-imitated helper functions
                if 'imitate_cpp' in fdict:
                    for ipyf in fdict['imitate_cpp']:
                        exec(ipyf)

                # cpp-based helper functions
                if 'cpp' in fdict:
                    for cppf in fdict['cpp']:
                        exec(cppf)
