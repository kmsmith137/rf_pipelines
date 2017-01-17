import sys
import numpy as np
import rf_pipelines

def clip_fx(intensity, weights, thr=3, n_internal=1, axis=None, dsample_nfreq=None, dsample_nt=None, imitate_cpp=False):
    """
    Helper function for intensity_clipper. Modifies 'weights' array in place.
    """
    
    (nfreq, nt_chunk) = intensity.shape

    # ------ Helper '__init__' calls ------
    assert axis in (None, 0, 1), "axis must be None (planar; freq and time), 0 (along freq; constant time), or 1 (along time; constant freq)."
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
        raise RuntimeError("clip_fx: current implementation requires 'dsample_nfreq' to be a divisor of stream nfreq.")
    if nt_chunk % dsample_nt != 0:
        raise RuntimeError("clip_fx: current implementation requires 'dsample_nt' to be a divisor of 'nt_chunk'.")
    
    # ------ Helper 'process_chunk' calls ------
    # Let's make a ref to the original high-resolution weights.
    weights_hres = weights
    
    if coarse_grained:
        # Downsample the weights and intensity.
        (intensity, weights) = rf_pipelines.wi_downsample(intensity, weights, dsample_nfreq, dsample_nt)

    if axis == None: 
        # In the 2D case, the behavior with and without imitate_cpp is the same.
        # Compute (mean,rms) values after 'n_internal' iterations while clipping at the level of 'thr'
        (mean, rms) = rf_pipelines.weighted_mean_and_rms(intensity, weights, n_internal, thr)
        clip = rf_pipelines.tile_arr(np.asarray(rms), axis, dsample_nfreq, dsample_nt)
        
        # Boolean array which is True for masked values
        mask = np.abs(intensity-mean) >= (thr * clip)

    elif imitate_cpp:
        # 1D case, with imitate_cpp=True
        (mean, rms) = rf_pipelines.weighted_mean_and_rms(intensity, weights, n_internal, thr, axis)

        mean = rf_pipelines.tile_arr(mean, axis, dsample_nfreq, dsample_nt)
        clip = rf_pipelines.tile_arr(rms, axis, dsample_nfreq, dsample_nt)

        # Boolean array which is True for masked values
        mask = np.abs(intensity-mean) >= (thr * clip)

    else:
        # 1D case, with imitate_cpp=False
        # Compute (sum_i W_i I_i^2) and (sum_i W_i)
        num = np.asarray(np.sum(weights*(intensity)**2, axis=axis))
        den = np.asarray(np.sum(weights, axis=axis))

        np.putmask(den, den==0., 1.0)     # replace 0.0 by 1.0 to avoid divide-by-zero

        clip = np.sqrt(num/den)
        clip = rf_pipelines.tile_arr(clip, axis, dsample_nfreq, dsample_nt)

        mask = np.abs(intensity) > (thr * clip)

    if coarse_grained:
        mask = rf_pipelines.upsample(mask, nfreq, nt_chunk)

    # Assign zero weights to those elements that have an
    # intensity value beyond the threshold limit.
    np.putmask(weights_hres, mask, 0.)

class intensity_clipper(rf_pipelines.py_wi_transform):
    """
    This transform clips the intensity along a selected 
    axis -- also works in planar (2D) mode -- and above 
    a given threshold. Results are applied to the weights 
    array (i.e., weights[clipped] = 0.) for masking 
    extreme values.

    Limitations:
        - 'n_internal' is only supported in 2D (axis=None) if 'imitate_cpp=False'
        - In 2D (axis=None), the same threshold value 'thr=3' is used for the internal 
          and external clipping loops. See 'intensity_clipper_cpp.py' for an advanced 
          implementation (in 2D and 1D) where the two 'thr' values need not be the same. 
        - FIXME Assumes zero mean (i.e., the intensity has already been detrended along 
          the selected axis).

    Constructor syntax:

      t = intensity_clipper(thr=3, n_internal=1, axis=None, nt_chunk=1024, dsample_nfreq=None, dsample_nt=None, test=False, imitate_cpp=False)

      'thr=3.' is the multiplicative factor of maximum threshold,
       e.g., 3 * standard_deviation, meaning that (the absolute
       value of) any sample above this limit is clipped.

      'n_internal=1' is the number of iterations used when calculating the weighted mean/rms before the final clip.

      'axis=None' is the axis convention:
        None: planar; freq and time. 
        0: along freq; constant time.
        1: along time; constant freq.

      'nt_chunk=1024' is the buffer size.

      'dsample_nfreq' and 'dsample_nt' are the downsampled 
       number of pixles along the freq and time axes, respectively.

      'test=False' enables a test mode.

      'imitate_cpp=False' enables an imitated (Python) version of C++ algorithms.
    """
    
    def __init__(self, thr=3., n_internal=1, axis=None, nt_chunk=1024, dsample_nfreq=None, dsample_nt=None, test=False, imitate_cpp=False):

        name = 'intensity_clipper(thr=%f, n_internal=%d, axis=%s, nt_chunk=%d' % (thr, n_internal, axis, nt_chunk)
        if dsample_nfreq is not None:
            name += ', dsample_nfreq=%d' % dsample_nfreq
        if dsample_nt is not None:
            name += ', dsample_nt=%d' % dsample_nt
        name += ')'

        self.thr = thr
        self.n_internal = n_internal
        self.axis = axis
        self.name = name
        self.nt_chunk = nt_chunk
        self.nt_prepad = 0
        self.nt_postpad = 0
        self.dsample_nfreq = dsample_nfreq
        self.dsample_nt = dsample_nt
        self.test = test
        self.imitate_cpp = imitate_cpp

    def set_stream(self, stream):
        self.nfreq = stream.nfreq

    def process_chunk(self, t0, t1, intensity, weights, pp_intensity, pp_weights):
        # If 'test' is specified, this will be a pseudo-transform which doesn't modify data
        # in the pipeline, but simulates some fake Gaussian data and reports the clipped fraction
        if self.test:
            intensity = np.random.normal(0, 1, size=intensity.shape)
            weights = np.ones(weights.shape)

        # Using clip_fx() mask the 'weights' in place.
        clip_fx(intensity, weights, self.thr, self.n_internal, self.axis, self.dsample_nfreq, self.dsample_nt, self.imitate_cpp)

        if self.test: 
            unmasked_percentage = np.count_nonzero(weights_hres) / float(weights_hres.size) * 100.
            print unmasked_percentage, "% not masked."
