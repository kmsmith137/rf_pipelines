import numpy as np
import rf_pipelines

class std_dev_filter(rf_pipelines.py_wi_transform):
    """
   Masks weights array based on the weighted 
   (intensity) std.dev deviating by some sigma.   
   
    Constructor syntax:

      t = std_dev_filter(thr=6., axis=0, nt_chunk=1024)

      FIXME: 'thr=3.' is the sigma value to clip (computed wrt the array of all channel
          standard deviations). Note: clipping based on absolute value of deviations

      'axis=0' is the axis convention:
        None: planar; freq and time. 
        0: along freq; constant time.
        1: along time; constant freq.

      'nt_chunk=1024' is the buffer size.
    """

    def __init__(self, thr=6., axis=None, nt_chunk=1024):
        
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
        # TODO 1d time, 2d global; weighted stdv
        sd = np.ma.masked_where(weights==0.,intensity,copy=False).std(self.axis)
        weights[abs(sd-sd.mean()) > self.thr*sd.std()] = 0.
