import numpy as np
from numpy import random

from rf_pipelines.rf_pipelines_c import pipeline_object, wi_transform
from rf_pipelines.utils import Variance_Estimates


class mask_filler(wi_transform):
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
        name = "mask_filler(var_file=%s, w_cutoff=%s, nt_chunk=%d)" % (var_file, w_cutoff, nt_chunk)

        # Call base class constructor
        wi_transform.__init__(self, name)

        self.Variance = Variance_Estimates(var_file)
        self.var_file = var_file
        self.w_cutoff = w_cutoff
        self.nt_chunk = nt_chunk


    def _bind_transform(self, json_attrs):
        if not json_attrs.has_key('t_initial'):
            raise RuntimeError("rf_pipelines.mask_filler: pipeline must contain a chime_file_stream (or another stream which defines the 't_initial' attribute)")
        
        self.t_initial = json_attrs['t_initial']
        self.dt_sample = json_attrs['dt_sample']


    def _process_chunk(self, intensity, weights, pos):
        # ---------------------------------------------------------
        # FIXME We need to normalize weights here so that
        # w_cutoff becomes more meaningful and less data-dependant.
        # e.g. weights[:] /= np.max(weights).
        # note: (1) Incoherent-beam acq --> w_max = 2.0
        #       (2) Any such changes may also need to be applied
        #           to other components (e.g. bonsai)
        # ---------------------------------------------------------
        
        t0 = self.t_initial + self.dt_sample * pos
        t1 = self.t_initial + self.dt_sample * (pos + nt_chunk)
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


    def jsonize(self):
        return { 'class_name': 'mask_filler',
                 'var_file': self.var_file,
                 'w_cutoff': self.w_cutoff,
                 'nt_chunk': self.nt_chunk }


    @staticmethod
    def from_json(j):
        return mask_filler(var_file = j['var_file'],
                           w_cutoff = j['w_cutoff'],
                           nt_chunk = j['nt_chunk'])


pipeline_object.register_json_constructor('mask_filler', mask_filler.from_json)
