import sys
import numpy as np

from rf_pipelines.rf_pipelines_c import wi_transform
from rf_pipelines.utils import tile_arr, upsample, weighted_mean_and_rms, wi_downsample


def filter_stdv(intensity, weights, thr=3, axis=1, dsample_nfreq=None, dsample_nt=None):
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
        (intensity, weights) = wi_downsample(intensity, weights, dsample_nfreq, dsample_nt)

    # Pass 1: Compute the weighted standard deviation of the intensity array along the 
    # selected axis.  We also build up a weights array 'sd_weights' which is 0 or 1.

    (mean, rms) = weighted_mean_and_rms(intensity, weights, axis=axis)
    sd = rms**2 # Compute the variance for saving computational cost in C++
    sd_weights = (rms > 0)

    # Pass 2: Compute the mean and rms of the 'sd' array, and clip based on it.
    # The result is a 2D boolean array 'mask' (clipped values are represented by True).
    #
    # Note: it would be trivial to implement iterated clipping, since
    # weighted_mean_and_rms() already has a 'niter' argument.  If this
    # turns out to be useful, then it would have tiny computational cost,
    # since the iteration is done when the array is 1D!
    
    (sd_mean, sd_rms) = weighted_mean_and_rms(sd, sd_weights, sigma_clip=thr)

    # 1D boolean mask (clipped values are represented by True)
    mask = np.logical_or((sd_weights == 0.0), (np.abs(sd-sd_mean) >= thr*sd_rms))

    # 2D boolean mask (still needs upsampling)
    mask = tile_arr(mask, axis, dsample_nfreq, dsample_nt)

    # Upsample to original resolution
    if coarse_grained:
        mask = upsample(mask, nfreq, nt_chunk)

    # Apply mask to original hi-res weights array
    np.putmask(weights_hres, mask, 0.)


class std_dev_clipper(wi_transform):
    """
    Masks weights array based on the weighted (intensity) 
    standard deviation deviating by some sigma. 
   
    Constructor syntax:

      t = std_dev_clipper(nt_chunk=1024, sigma=3, axis=1, Df=1, Dt=1, two_pass=False)

      'nt_chunk=1024' is the buffer size.
      
      'sigma=3.' is the sigma value to clip. 

      'axis=1' is the axis convention:
        0: along freq; constant time.
        1: along time; constant freq.

      (Df,Dt)=(1,1) are the downsampling factors in frequency, time.

      The 'two_pass' argument is ignored by the python implementation, but is included
      as a dummy constructor argument, so that the C++ and python clippers will have
      the same syntax.
    """


class std_dev_clipper(wi_transform):
    def __init__(self, nt_chunk=1024, sigma=3, axis=1, Df=1, Dt=1, two_pass=False):
        name = 'std_dev_clipper_python(thr=%f, axis=%s, nt_chunk=%d, Df=%d, Dt=%d)' % (thr, axis, nt_chunk, Df, Dt)
        wi_transform.__init__(self, name)
        
        self.thr = sigma
        self.axis = axis
        self.nt_chunk  = nt_chunk
        self.Df = Df
        self.Dt = Dt

        # self.dsample_nt can be initialized here
        # self.dsample_nfreq will be initialized in _bind_transform()
        self.dsample_nt = nt_chunk // Dt


    def _bind_transform(self, json_attrs):
        assert self.nfreq % self.Df == 0
        self.dsample_nfreq = self.nfreq // self.Df


    def _process_chunk(self, intensity, weights, pos):
        filter_stdv(intensity, weights, self.thr, self.axis, self.dsample_nfreq, self.dsample_nt)
