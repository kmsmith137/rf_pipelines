import numpy as np
from numpy import random
import rf_pipelines


class noise_filler(rf_pipelines.py_wi_transform):
    """
    
    """

    def __init__(self, var_files, n_varsamples, w_cutoff, nt_chunk):
        name = "mask_filler(w_cutoff=%d, nt_chunk=%d)" % (w_cutoff, nt_chunk)

        # Call base class constructor
        rf_pipelines.py_wi_transform.__init__(self, name)

        # Sort, load, stitch
        sorted_plots = sorted(var_files)
        arrays = map(np.load, var_files)
        concatenated = np.hstack((arrays))
        self.var = concatenated
        
        # Initialize some other values
        self.n_varsamples = n_varsamples   # number of samples per data point in variance array
        self.w_cutoff = w_cutoff
        self.nt_postpad = 0
        self.nt_prepad = 0

        # FOR NOW nt_chunk must be a multiple of n_varsamples
        assert nt_chunk % n_varsamples == 0, \
            'For now, nt_chunk(=%d) must be a multiple of n_varsamples(=%d). I might implement some buffering later to fix this!' % (nt_chunk, n_varsamples)
        self.nt_chunk = nt_chunk


    def set_stream(self, s):
         self.nfreq = s.nfreq


    def start_substream(self, isubstream, t0):
        # Can add some sort of buffer here later to remove need to nt_chunk % n_varsamples == 0 if desired
        pass


    def process_chunk(self, t0, t1, intensity, weights, pp_intensity, pp_weights):
        ivariance = 0
        imaxvar = self.var.shape[1]
        for frequency in xrange(self.nfreq):
            for i in xrange(intensity.shape[1]):
                if weights[frequency, i] > self.w_cutoff:
                    weights[frequency, i] = 2.0
                else:
                    sigma = (self.var[frequency, ivariance])**2
                    if sigma == 0:
                        weights[frequency, i] = 0
                    else:
                        intensity[frequency, i] = sigma * np.random.standard_normal()
                        weights[frequency, i] = 2.0
                if i % self.n_varsamples == 0 and i != 0:
                    ivariance += 1
            ivariance = 0

    def end_substream(self):
        # Write out plot
        pass

