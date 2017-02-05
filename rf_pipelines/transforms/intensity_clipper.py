import sys
import numpy as np

import rf_pipelines
from rf_pipelines import rf_pipelines_c


def intensity_clipper(nt_chunk=1024, sigma=3., axis=None, niter=1, iter_sigma=0., Df=1, Dt=1, two_pass=False, cpp=True, test=False):
    """
    This transform clips the intensity along a selected 
    axis -- also works in planar (2D) mode -- and above 
    a given threshold. Results are applied to the weights 
    array (i.e., weights[clipped] = 0.) for masking 
    extreme values.

    Constructor syntax:

      t = intensity_clipper(nt_chunk=1024, sigma=3., axis=None, niter=1, iter_sigma=0., Df=1, Dt=1, two_pass=False, cpp=True, test=False)

      'nt_chunk=1024' is the buffer size.

      'sigma=3.' is the multiplicative factor of maximum threshold,
       e.g., 3 * standard_deviation, meaning that (the absolute
       value of) any sample above this limit is clipped.

      'axis=None' is the axis convention:
        None: planar; freq and time. 
        0: along freq; constant time.
        1: along time; constant freq.

      'niter=1' is the number of iterations used when calculating the weighted mean/rms before the final clip.

      'iter_sigma=0' is the internal threshold used when iterating to the final weight mean/rms.
         If iter_sigma=0 is specified, then it defaults to 'sigma'
    
      (Df,Dt)=(1,1) are the downsampling factors in frequency, time.

      'cpp=True' will use fast C++ transforms
      'cpp=False' will use reference python transforms

      If 'two_pass=True' then a more numerically stable but slightly slower clipping algorithm
      will be used (only meaningful if cpp=True).

      'test=False' enables a test mode (only meaningful if cpp=False)
    """

    if cpp:
        return rf_pipelines_c.make_intensity_clipper(nt_chunk, axis, sigma, niter, iter_sigma, Df, Dt, two_pass)

    if (iter_sigma != 0) and (iter_sigma != sigma):
        print >>sys.stderr, 'rf_pipelines intensity_clipper(): warning: iter_sigma argument is currently ignored by python transform'

    return intensity_clipper_python(sigma, niter, axis, nt_chunk, Df, Dt, test)


def clip_fx(intensity, weights, thr=3, n_internal=1, axis=None, dsample_nfreq=None, dsample_nt=None):
    """Helper function for intensity_clipper_python. Modifies 'weights' array in place."""
    
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
        # 2D case: Compute (mean,rms) values after 'n_internal' iterations while clipping at the level of 'thr'
        (mean, rms) = rf_pipelines.weighted_mean_and_rms(intensity, weights, n_internal, thr)
        clip = rf_pipelines.tile_arr(np.asarray(rms), axis, dsample_nfreq, dsample_nt)
        
        # Boolean array which is True for masked values
        mask = np.abs(intensity-mean) >= (thr * clip)

    else:
        # 1D case
        (mean, rms) = rf_pipelines.weighted_mean_and_rms(intensity, weights, n_internal, thr, axis)

        mean = rf_pipelines.tile_arr(mean, axis, dsample_nfreq, dsample_nt)
        clip = rf_pipelines.tile_arr(rms, axis, dsample_nfreq, dsample_nt)

        # Boolean array which is True for masked values
        mask = np.abs(intensity-mean) >= (thr * clip)

    if coarse_grained:
        mask = rf_pipelines.upsample(mask, nfreq, nt_chunk)

    # Assign zero weights to those elements that have an
    # intensity value beyond the threshold limit.
    np.putmask(weights_hres, mask, 0.)


class intensity_clipper_python(rf_pipelines.py_wi_transform):
    """
    Limitation: In 2D (axis=None), the same threshold value 'thr=3' is used for the internal 
                and external clipping loops. See 'intensity_clippers.cpp' for an advanced 
                implementation (in 2D and 1D) where the two 'thr' values need not be the same. 
    """
    
    def __init__(self, thr=3., n_internal=1, axis=None, nt_chunk=1024, Df=1, Dt=1, test=False):
        name = 'intensity_clipper_python(thr=%f, n_internal=%d, axis=%s, nt_chunk=%d, Df=%d, Dt=%d)' % (thr, n_internal, axis, nt_chunk, Df, Dt)
        rf_pipelines.py_wi_transform.__init__(self, name)

        self.thr = thr
        self.n_internal = n_internal
        self.axis = axis
        self.name = name
        self.nt_chunk = nt_chunk
        self.nt_prepad = 0
        self.nt_postpad = 0
        self.Df = Df
        self.Dt = Dt
        self.test = test
        
        assert Df > 0
        assert Dt > 0
        assert nt_chunk > 0
        assert nt_chunk % Dt == 0

        # self.dsample_nt can be initialized here
        # self.dsample_nfreq will be initialized in set_stream()
        self.dsample_nt = nt_chunk // Dt


    def set_stream(self, stream):
        assert stream.nfreq % self.Df == 0

        self.nfreq = stream.nfreq
        self.dsample_nfreq = stream.nfreq // self.Df


    def process_chunk(self, t0, t1, intensity, weights, pp_intensity, pp_weights):
        # If 'test' is specified, this will be a pseudo-transform which doesn't modify data
        # in the pipeline, but simulates some fake Gaussian data and reports the clipped fraction
        if self.test:
            intensity = np.random.normal(0, 1, size=intensity.shape)
            weights = np.ones(weights.shape)

        # Using clip_fx() mask the 'weights' in place.
        clip_fx(intensity, weights, self.thr, self.n_internal, self.axis, self.dsample_nfreq, self.dsample_nt)

        if self.test: 
            unmasked_percentage = np.count_nonzero(weights_hres) / float(weights_hres.size) * 100.
            print unmasked_percentage, "% not masked."
