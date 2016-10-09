import numpy as np
import rf_pipelines


def expand_mask(weights, thr, axis):
    """Helper function for mask_expander.  Modifies 'weights' array in place."""

    (nfreq, nt_chunk) = weights.shape

    # Compute the mean of weights along the seleted axis.
    w_mean = np.mean(weights, axis=axis)

    if axis is None:
        if w_mean <= thr:
            weights[:] = 0.
    else:
        w_mean = rf_pipelines.tile_arr(w_mean, axis, nfreq, nt_chunk)
        np.putmask(weights, w_mean <= thr, 0.)


class mask_expander(rf_pipelines.py_wi_transform):
    """
   This transform expands the mask by filling the weights 
   array, provided that its mean is less than or equal to 
   a given threshold, with zeros along a selected axis.

   + Assumes that weights are between 0 and 1.

    Constructor syntax:

      t = mask_expander(thr=0.2, axis=None, nt_chunk=1024, test=False)
      
      'thr=0.2' should be between 0 and 1. Any selected
        weights sub-array with a mean value less than or 
        equal to this threshold is totally masked by the 
        transform.

      'axis=None' is the axis convention:
        None: planar; freq and time.
        0: along freq; constant time.
        1: along time; constant freq.
      
      'nt_chunk=1024' is the buffer size.
    """

    def __init__(self, thr=0.2, axis=None, nt_chunk=1024):
        
        assert (0 < thr < 1), "threshold must be between 0 and 1."
        assert (axis == None) or (axis == 0) or (axis == 1),\
                "axis must be None (planar; freq and time), 0 (along freq; constant time), or 1 (along time; constant freq)."

        self.thr = thr
        self.axis = axis
        self.nt_chunk = nt_chunk
        self.name = 'mask_expander(thr=%f, axis=%s, nt_chunk=%d)' % (thr, axis, nt_chunk)

    def set_stream(self, stream):        
        self.nfreq = stream.nfreq

    def process_chunk(self, t0, t1, intensity, weights, pp_intensity, pp_weights):
        expand_mask(weights, self.thr, self.axis)
