import rf_pipelines
import numpy as np
from numpy import random


class mask_filler(rf_pipelines.py_wi_transform):
    """
    Modifies values in the intensity and weight arrays. If the weight is > w_cutoff, the weight is changed to 
    2.0 and the intensity is left unmodified. Otherwise the weight is changed to 2.0 AND the intensity 
    is replaced with gaussian random noise with standard deviation given by a previously calculated variance array.

    Note that if the entire frequency channel was masked in the variance array (variance = 0), it will remain 
    masked. 
    
    Constructor Arguments
    ----------------------
    var_file - the h5 variance file for the acquisition
    w_cutoff - weight cutoff above which the weight will not be replaced by random noise
    nt_chunk - the buffer size
    """

    def __init__(self, var_file, w_cutoff, nt_chunk=1024):
        name = "mask_filler(var_file=%s, w_cutoff=%d, nt_chunk=%d)" % (var_file, w_cutoff, nt_chunk)

        # Call base class constructor
        rf_pipelines.py_wi_transform.__init__(self, name)

        self.Variance = rf_pipelines.utils.Variance_Estimates(var_file)
        self.w_cutoff = w_cutoff
        self.nt_postpad = 0
        self.nt_prepad = 0
        self.nt_chunk = nt_chunk
        print 'WARNING nt_chunk should be less than or equal to v1_chunk * v2_chunk for the variance array.'

    def set_stream(self, s):
         self.nfreq = s.nfreq

    def process_chunk(self, t0, t1, intensity, weights, pp_intensity, pp_weights):
        
        # ---------------------------------------------------------
        # FIXME We need to normalize weights here so that
        # w_cutoff becomes more meaningful and less data-dependant.
        # e.g. weights[:] /= np.max(weights).
        # note: (1) Incoherent-beam acq --> w_max = 2.0
        #       (2) Any such changes may also need to be applied
        #           to other components (e.g. bonsai)
        # ---------------------------------------------------------
        var = self.Variance.eval((t0+t1)/2.)
        
        # 'intensity_valid' will be a 2D boolean-valued numpy array
        intensity_valid = (weights > self.w_cutoff)
        
        rand_intensity = np.random.standard_normal(size=intensity.shape)
        weights[:,:] = 0.0
        
        for (ifreq,v) in enumerate(var):
            if v > 0.0:
                rand_intensity[ifreq,:] *= v**0.5
                weights[ifreq,:] = 2.0

        intensity[:,:] = np.where(intensity_valid, intensity, rand_intensity)
