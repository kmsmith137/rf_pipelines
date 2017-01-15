import numpy as np
import rf_pipelines

def filter_stdv(intensity, weights, thr=3, axis=1, dsample_nfreq=None, dsample_nt=None, imitate_cpp=False):
    """Helper function for std_dev_clipper. Modifies 'weights' array in place."""
    
    (nfreq, nt_chunk) = intensity.shape
    
    # ------ Helper '__init__' calls ------
    assert axis in (0, 1), "axis must be 0 (along freq; constant time), or 1 (along time; constant freq)."
    assert thr >= 1., "threshold must be >= 1."
    assert nt_chunk > 0
    assert (dsample_nt is None or dsample_nt > 0), "Invalid downsampling number along the time axis!"
    assert (dsample_nfreq is None or dsample_nfreq > 0), "Invalid downsampling number along the freq axis!"

    # ------ Helper 'set_stream' calls ------
    coarse_grained = (dsample_nfreq < nfreq) or (dsample_nt < nt_chunk)

    if dsample_nfreq is None:
        dsample_nfreq = nfreq
    if dsample_nt is None:
        dsample_nt = nt_chunk

    if nfreq % dsample_nfreq != 0:
        raise RuntimeError("filter_stdv: current implementation requires 'dsample_nfreq' to be a divisor of stream nfreq.")
    if nt_chunk % dsample_nt != 0:
        raise RuntimeError("filter_stdv: current implementation requires 'dsample_nt' to be a divisor of 'nt_chunk'.")
    
    # ------ Helper 'process_chunk' calls ------
    # Let's make a ref to the original high-resolution weights.
    weights_hres = weights

    if coarse_grained:
        # Downsample the weights and intensity.
        (intensity, weights) = rf_pipelines.wi_downsample(intensity, weights, dsample_nfreq, dsample_nt)

    # Pass 1: Compute the weighted standard deviation of the intensity array along the 
    # selected axis.  We also build up a weights array 'sd_weights' which is 0 or 1.

    if imitate_cpp:
        (mean, rms) = rf_pipelines.weighted_mean_and_rms(intensity, weights, axis=axis)
        sd = rms**2
        sd_weights = (rms > 0)

    else:
        num = np.asarray(np.sum(weights*(intensity)**2, axis=axis))
        den = np.asarray(np.sum(weights, axis=axis))
        np.putmask(den, den==0., 1.0)
        sd = np.sqrt(num/den)
        sd_weights = np.ones_like(sd)

    # Pass 2: Compute the mean and rms of the 'sd' array, and clip based on it.
    # The result is a 2D boolean array 'mask' (clipped values are represented by True).

    if imitate_cpp:

        # Code cleanup: compute the mean/rms by calling weighted_mean_and_rms().
        #
        # Note: it would be trivial to implement iterated clipping, since
        # weighted_mean_and_rms() already has a 'niter' argument.  If this
        # turns out to be useful, then it would have tiny computational cost,
        # since the iteration is done when the array is 1D!
        
        (sd_mean, sd_rms) = rf_pipelines.weighted_mean_and_rms(sd, sd_weights, sigma_clip=thr)

        # 1D boolean mask (clipped values are represented by True)
        mask = np.logical_or((sd_weights == 0.0), (np.abs(sd-sd_mean) >= thr*sd_rms))

        # 2D boolean mask (still needs upsampling)
        mask = rf_pipelines.tile_arr(mask, axis, dsample_nfreq, dsample_nt)

    else:
        # Tile 'sd' so that it matches with the shape of intensity.
        sd = rf_pipelines.tile_arr(sd, axis, dsample_nfreq, dsample_nt)

        # This block of code creates a mask, and hence filters, 
        # based on the mean and stdv of 'sd' along THE OTHER axis.
        axis = np.abs(axis-1)
        sd_stdv = rf_pipelines.tile_arr(sd.std(axis=axis), axis, dsample_nfreq, dsample_nt)
        sd_mean = rf_pipelines.tile_arr(sd.mean(axis=axis), axis, dsample_nfreq, dsample_nt)

        assert sd.shape == sd_mean.shape == sd_stdv.shape

        # Boolean array which is True for masked values
        mask = np.abs(sd-sd_mean) > (thr*sd_stdv)

    # Upsample to original resolution
    if coarse_grained:
        mask = rf_pipelines.upsample(mask, nfreq, nt_chunk)

    # Apply mask to original hi-res weights array
    np.putmask(weights_hres, mask, 0.)

class std_dev_clipper(rf_pipelines.py_wi_transform):
    """
    Masks weights array based on the weighted (intensity) 
    standard deviation deviating by some sigma. 
   
    Constructor syntax:

      t = std_dev_clipper(thr=3., axis=1, nt_chunk=1024)

      'thr=3.' is the sigma value to clip. 

      'axis=1' is the axis convention:
        0: along freq; constant time.
        1: along time; constant freq.

      'nt_chunk=1024' is the buffer size.
      
      'dsample_nfreq' and 'dsample_nt' are the downsampled
       number of pixles along the freq and time axes, respectively.
    """
    
    def __init__(self, thr=3., axis=1, nt_chunk=1024, dsample_nfreq=None, dsample_nt=None):
        
        self.thr = thr
        self.axis = axis
        self.nt_chunk  = nt_chunk
        self.nt_prepad = 0
        self.nt_postpad = 0
        self.dsample_nfreq = dsample_nfreq
        self.dsample_nt = dsample_nt

        name = 'std_dev_clipper(thr=%f, axis=%s, nt_chunk=%d' % (thr, axis, nt_chunk)
        if dsample_nfreq is not None:
            name += ', dsample_nfreq=%d' % dsample_nfreq
        if dsample_nt is not None:
            name += ', dsample_nt=%d' % dsample_nt
        name += ')'
        self.name = name

    def set_stream(self, stream):
        self.nfreq = stream.nfreq

    def process_chunk(self, t0, t1, intensity, weights, pp_intensity, pp_weights):
        filter_stdv(intensity, weights, self.thr, self.axis, self.dsample_nfreq, self.dsample_nt)
