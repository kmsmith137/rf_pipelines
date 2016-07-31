import numpy as np
import rf_pipelines

class clipper_transform(rf_pipelines.py_wi_transform):
    """
   Clips the intensity above a given stdv threshold. 
   Results are applied to the weights array for masking RFIs.

    Constructor syntax:

      t = clipper_transform(thr=3, axis=0, nt_chunk=1024)

      'thr=3.' is the multiplicative factor of maximum threshold,
        e.g., 3 * standard_deviation, meaning that anything above this
        limit is clipped.

      'axis=0' is the axis convention:
        0: along freq; constant time
        1: along time; constant freq

      'nt_chunk=1024' is the buffer size.
    """

    def __init__(self, thr=3., axis=0, nt_chunk=1024):
        
        self.thr = thr
        self.axis = axis
        self.nt_chunk = nt_chunk
        self.nt_prepad = 0
        self.nt_postpad = 0

        assert (self.axis == 0 or self.axis == 1), \
        'axis must be 0 (along freq; constant time) or 1 (along time; constant freq).'

    def set_stream(self, stream):

        self.nfreq = stream.nfreq

    def process_chunk(self, t0, t1, intensity, weights, pp_intensity, pp_weights):
        
        self.clip = np.std(intensity, axis=self.axis)
         
        if self.axis == 0:
            self.clip = np.tile(self.clip, (self.nfreq, 1))
        else:
            self.clip = np.transpose(np.tile(self.clip, (self.nt_chunk, 1)))
        
        assert np.shape(weights) == np.shape(intensity) == np.shape(self.clip)
        weights[ intensity > (self.thr*self.clip) ] = 0.
