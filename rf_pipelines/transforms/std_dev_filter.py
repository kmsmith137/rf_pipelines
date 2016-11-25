import numpy as np
import rf_pipelines

def filter_stdv(intensity, weights, thr, axis):
    """Helper function for std_dev_filter. Modifies 'weights' array in place."""

    (nfreq, nt_chunk) = intensity.shape
    
    # Compute the weighted standard deviation of the intensity
    # array along the selected axis.
    num = np.asarray(np.sum(weights*(intensity)**2, axis=axis))
    den = np.asarray(np.sum(weights, axis=axis))
    
    np.putmask(den, den==0., 1.0)
    sd = np.sqrt(num/den)

    # Tile 'sd' so that it matches with intensity.shape.
    sd = rf_pipelines.tile_arr(sd, axis, nfreq, nt_chunk)

    # This block of code creates a mask, and hence filters, 
    # based on the mean and stdv of 'sd' along the other axis.
    axis = np.abs(axis-1)
    sd_stdv = rf_pipelines.tile_arr(sd.std(axis=axis), axis, nfreq, nt_chunk)
    sd_mean = rf_pipelines.tile_arr(sd.mean(axis=axis), axis, nfreq, nt_chunk)

    np.putmask(weights, np.abs(sd-sd_mean) > (thr*sd_stdv), 0.)

class std_dev_filter(rf_pipelines.py_wi_transform):
    """
   Masks weights array based on the weighted (intensity) 
   standard deviation deviating by some sigma. 
   
    Constructor syntax:

      t = std_dev_filter(thr=3., axis=None, nt_chunk=1024)

      'thr=3.' is the sigma value to clip. 

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
        filter_stdv(intensity, weights, self.thr, self.axis)
