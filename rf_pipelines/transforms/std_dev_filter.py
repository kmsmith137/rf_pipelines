import numpy as np
import rf_pipelines

class std_dev_filter(rf_pipelines.py_wi_transform):
    """
   FIXME: Masks weights array based on the weighted 
   (intensity) std.dev deviating by some sigma.   
   
    Constructor syntax:

      FIXME: t = std_dev_filter(thr=3., axis=None, nt_chunk=1024)

      FIXME: 'thr=3.' is the sigma value to clip (computed wrt the array of all channel
          standard deviations). Note: clipping based on absolute value of deviations

      'axis=0' is the axis convention:
        None: planar; freq and time. 
        0: along freq; constant time.
        1: along time; constant freq.

      'nt_chunk=1024' is the buffer size.
    """

    def __init__(self, thr=3., axis=None, nt_chunk=1024):
        
        assert thr >= 1., "threshold must be >= 1."
        assert (axis == None or axis == 0 or axis == 1),\
            "axis must be None (planar; freq and time), 0 (along freq; constant time), or 1 (along time; constant freq)."
        assert nt_chunk > 0

        self.thr = thr
        self.axis = axis
        self.nt_chunk  = nt_chunk
        self.nt_prepad = 0
        self.nt_postpad = 0

        name = 'std_dev_filter(thr=%f, axis=%s, nt_chunk=%d' % (thr, axis, nt_chunk)
        self.name = name

    def set_stream(self,stream):
        self.nfreq = stream.nfreq

    def process_chunk(self, t0, t1, intensity, weights, pp_intensity, pp_weights):
        
        num = np.asarray(np.sum(weights*(intensity)**2, axis=self.axis))
        den = np.asarray(np.sum(weights, axis=self.axis))

        np.putmask(den, den==0., 1.0)     # replace 0.0 by 1.0 to avoid divide-by-zero
        sd = np.sqrt(num/den)
        sd = rf_pipelines.tile_arr(sd, self.axis, self.nfreq, self.nt_chunk)
        
        if self.axis == 0:
            mask_axis = 1
        if self.axis == 1:
            mask_axis = 0
        
        mask = np.abs(sd-sd.mean(axis=mask_axis)) > (self.thr * rf_pipelines.tile_arr(sd.std(axis=mask_axis), mask_axis, self.nfreq, self.nt_chunk))
        assert mask.shape == weights.shape
        np.putmask(weights, mask, 0.)
