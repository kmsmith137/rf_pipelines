import numpy as np
from types import DictType, ListType
import rf_pipelines
from rf_pipelines import rf_pipelines_c

class master_transform(rf_pipelines.py_wi_transform):
    """
    This transform utilizes HELPER FUNCTIONS to perform 
    an iterated chian of transforms, which include 
    detrending and clipping algorithms. The iteraion is 
    controlled (hence optimized) by user-defined parameters.

    Constructor syntax:

      t = master_transform(nt_chunk=1024, fdict=None, rms_cut=0., mask_cut=0.05, max_niter=1, test=False):
      
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
       are set to zero and further iterations are broken.

      'mask_cut=0.05' is the masking thershold for the entire chunk.
       Iterations proceed only if two consecutive iterations result
       in a fraction-of-unmasked difference greater than this value.
    
      'max_niter=1' is the maximum number of iterations for each chunk.

      'test=False' triggers the test flag.
    """

    def __init__(self, nt_chunk=1024, fdict=None, rms_cut=0., mask_cut=0.05, max_niter=1, test=False):
        
        assert nt_chunk > 0
        assert (type(fdict) is DictType) and ({key for key in fdict.keys()} == {'py', 'imitate_cpp', 'cpp'}),\
            "master_transform: 'fdict' must be a dictionary with the following format:\n fdict = {'py':[], 'imitate_cpp':[], 'cpp':[]}"
        for value in fdict.values():
            assert type(value) is ListType
        assert rms_cut >= 0., "master_transform: rms threshold must be >= 0."
        assert 0.0 < mask_cut < 0.1
        assert max_niter >= 1
        assert type(test) == bool

        self.nt_chunk = nt_chunk
        self.fdict = fdict
        self.rms_cut = rms_cut
        self.mask_cut = mask_cut
        self.max_niter = max_niter
        self.test = test
        self.name = 'master_transform(nt_chunk=%d, *(fdict)=%d, rms_cut=%f, mask_cut=%f, max_niter=%d)'\
            % (nt_chunk, sum(map(len, fdict.values())), rms_cut, mask_cut, max_niter)

    def set_stream(self, stream):
        
        self.nfreq = stream.nfreq
        
    def process_chunk(self, t0, t1, intensity, weights, pp_intensity, pp_weights):
        
        for ix in xrange(self.max_niter):

            (mean, rms) = rf_pipelines_c.weighted_mean_and_rms(intensity, weights)
            unmasked_before = np.count_nonzero(weights) / float(weights.size)

            if rms > self.rms_cut:
                weights[:] = 0.
                break
            
            elif (ix > 0) and (abs(unmasked_before - unmasked_after) < self.mask_cut):
                break

            else:
                for i in self.fdict.items():
                    key = i[0]
                    if self.test:
                        pass # FIXME compare keys.output
                    else:
                        if not key:
                            pass # FIXME empty value
                        else:
                            for fx in self.fdict[key]: # FIXME keys.order
                                exec(fx)

            unmasked_after = np.count_nonzero(weights) / float(weights.size)
