import numpy as np
import rf_pipelines
from astropy.stats import median_absolute_deviation

class clipper_transform(rf_pipelines.py_wi_transform):
    """
   Clips the intensity above a given threshold. Results
   are applied to the weights array for masking RFIs.

    Constructor syntax:

      t = clipper_transform(thr=3, axis=0, func='stdv', nt_chunk=1024)

      'thr=3.' is the multiplicative factor of maximum threshold,
        e.g., 3 * standard_deviation, meaning that anything above this
        limit is clipped.

      'axis=0' is the axis convention:
        0: along freq; constant time
        1: along time; constant freq

      "func='stdv'" indicates the mathematical function that is used 
        for clipping:
          'stdv' : standard deviation (from the mean)
          'mad' : median absolute deviation

      'nt_chunk=1024' is the buffer size.
    """

    def __init__(self, thr=3., axis=0, func='stdv', nt_chunk=1024):
        
        self.thr = thr
        self.axis = axis
        self.func = func

        self.nt_chunk = nt_chunk
        self.nt_prepad = 0
        self.nt_postpad = 0

        assert (self.axis == 0 or self.axis == 1), \
        'Axis must be 0 (along freq; constant time) or 1 (along time; constant freq)'

    def set_stream(self, stream):

        self.nfreq = stream.nfreq

    def process_chunk(self, t0, t1, intensity, weights, pp_intensity, pp_weights):
        
        # FIXME: Here I was thinking of making a dictionary
        # of clipper functions! But is the inner func really 
        # necessary? The difficulty is that the intensity gets 
        # called in process_chunk().

        def func(intensity, axis, func):
            return {
                'stdv': np.std(intensity, axis=axis),
                'mad': median_absolute_deviation(intensity, axis=axis),
            } [func]
        
        self.clip = func(intensity, self.axis, self.func)
        
        if self.axis == 0:
            assert np.size(self.clip) == self.nt_chunk
            for time in xrange(self.nt_chunk):
                weights[ intensity[:,time] > (self.thr*self.clip[time]), time ] = 0.
        else:
            assert np.size(self.clip) == self.nfreq
            for freq in xrange(self.nfreq):
                weights[ freq, intensity[freq,:] > (self.thr*self.clip[freq]) ] = 0.
