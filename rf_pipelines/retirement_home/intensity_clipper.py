import sys
import numpy as np

from rf_pipelines.rf_pipelines_c import wi_transform
from rf_pipelines.utils import tile_arr, upsample, weighted_mean_and_rms, wi_downsample


def clip_fx(intensity, weights, thr=3, n_internal=1, axis=None, dsample_nfreq=None, dsample_nt=None):
    """Helper function for python intensity_clipper. Modifies 'weights' array in place."""
    
    (nfreq, nt_chunk) = intensity.shape

    assert axis in (None, 0, 1), "axis must be None (planar; freq and time), 0 (along freq; constant time), or 1 (along time; constant freq)."
    assert thr >= 1., "threshold must be >= 1."
    assert nt_chunk > 0
    assert (dsample_nt is None or dsample_nt > 0), "Invalid downsampling number along the time axis!"
    assert (dsample_nfreq is None or dsample_nfreq > 0), "Invalid downsampling number along the freq axis!"


    coarse_grained = (dsample_nfreq < nfreq) or (dsample_nt < nt_chunk)

    if dsample_nfreq is None:
        dsample_nfreq = nfreq
    if dsample_nt is None:
        dsample_nt = nt_chunk

    if nfreq % dsample_nfreq != 0:
        raise RuntimeError("clip_fx: current implementation requires 'dsample_nfreq' to be a divisor of stream nfreq.")
    if nt_chunk % dsample_nt != 0:
        raise RuntimeError("clip_fx: current implementation requires 'dsample_nt' to be a divisor of 'nt_chunk'.")


    # Let's make a ref to the original high-resolution weights.
    weights_hres = weights
    
    if coarse_grained:
        # Downsample the weights and intensity.
        (intensity, weights) = wi_downsample(intensity, weights, dsample_nfreq, dsample_nt)

    # Compute (mean, rms) values after 'n_internal' iterations while clipping at the level of 'thr'
    (mean, rms) = weighted_mean_and_rms(intensity, weights, n_internal, thr, axis)
    
    mean = tile_arr(np.asarray(mean), axis, dsample_nfreq, dsample_nt)
    clip = tile_arr(np.asarray(rms), axis, dsample_nfreq, dsample_nt)

    # Boolean array which is True for masked values
    mask = np.abs(intensity-mean) >= (thr * clip)

    if coarse_grained:
        mask = upsample(mask, nfreq, nt_chunk)

    # Assign zero weights to those elements that have an
    # intensity value beyond the threshold limit.
    np.putmask(weights_hres, mask, 0.)


class intensity_clipper(wi_transform):
    """
    This transform clips the intensity along a selected 
    axis -- also works in planar (2D) mode -- and above 
    a given threshold. Results are applied to the weights 
    array (i.e., weights[clipped] = 0.) for masking 
    extreme values.

    Constructor syntax:

      t = intensity_clipper(nt_chunk=1024, sigma=3., axis=None, niter=1, iter_sigma=0., Df=1, Dt=1, two_pass=False)

      'nt_chunk=1024' is the buffer size.

      'sigma=3.' is the multiplicative factor of maximum threshold,
       e.g., 3 * standard_deviation, meaning that (the absolute
       value of) any sample above this limit is clipped.

      'axis=None' is the axis convention:
        None: planar; freq and time. 
        0: along freq; constant time.
        1: along time; constant freq.

      'niter=1' is the number of iterations used when calculating the weighted mean/rms before the final clip.
    
      (Df,Dt)=(1,1) are the downsampling factors in frequency, time.

      The 'iter_sigma' and 'two_pass' arguments are ignored by the python reference implementation,
      but are included as dummy arguments, so that the C++ and python intensity_clippers will
      have the same syntax.

    Limitation: In 2D (axis=None), the same threshold value 'thr=3' is used for the internal 
                and external clipping loops. See 'intensity_clippers.cpp' for an advanced 
                implementation (in 2D and 1D) where the two 'thr' values need not be the same. 
    """
    
    def __init__(nt_chunk=1024, sigma=3., axis=None, niter=1, iter_sigma=0., Df=1, Dt=1, two_pass=False):
        if (iter_sigma != 0) and (iter_sigma != sigma):
            print >>sys.stderr, 'rf_pipelines intensity_clipper(): warning: iter_sigma argument is currently ignored by python transform'
            
        name = 'intensity_clipper_python(thr=%f, n_internal=%d, axis=%s, nt_chunk=%d, Df=%d, Dt=%d)' % (thr, n_internal, axis, nt_chunk, Df, Dt)
        wi_transform.__init__(self, name)

        self.thr = thr
        self.n_internal = niter
        self.axis = axis
        self.name = name
        self.nt_chunk = nt_chunk
        self.Df = Df
        self.Dt = Dt
        
        assert Df > 0
        assert Dt > 0
        assert nt_chunk > 0
        assert nt_chunk % Dt == 0

        # self.dsample_nt can be initialized here
        # self.dsample_nfreq will be initialized in _bind_transform()
        self.dsample_nt = nt_chunk // Dt

        
    def _bind_transform(self, json_attrs):
        assert self.nfreq % self.Df == 0
        self.dsample_nfreq = self.nfreq // self.Df


    def _process_chunk(self, intensity, weights, pos):
        # Using clip_fx() mask the 'weights' in place.
        clip_fx(intensity, weights, self.thr, self.n_internal, self.axis, self.dsample_nfreq, self.dsample_nt)
