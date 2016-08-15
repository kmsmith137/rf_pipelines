import numpy as np
import rf_pipelines

class clipper_transform(rf_pipelines.py_wi_transform):
    """
   This transform clips the intensity along a selected 
   axis -- also works in planar (2d) mode -- and above 
   a given threshold. Results are applied to the weights 
   array (i.e., weights[clipped] = 0.) for masking 
   extreme values.

   + Assumes zero mean (i.e., the intensity has already 
   been detrended along the selected axis).
   + Currently based on the weighted standard deviation 
   as explained in "chime_zerodm_notes".
   + Available in a coarse-grained mode by using 
   'dsample_nfreq', and 'dsample_nt'.
    
    Constructor syntax:

      t = clipper_transform(thr=3, axis=2, nt_chunk=1024,\
          dsample_nfreq=None, dsample_nt=None, test=False)

      'thr=3.' is the multiplicative factor of maximum threshold,
        e.g., 3 * standard_deviation, meaning that (the absolute
        value of) any sample above this limit is clipped.

      'axis=0' is the axis convention:
        0: along freq; constant time.
        1: along time; constant freq.
        2: planar; freq and time.

      'nt_chunk=1024' is the buffer size.

      'dsample_nfreq' and 'dsample_nt' are the downsampled 
       number of pixles along the freq and time axes, respectively.

      'test=False' enables a test mode.
    """

    def __init__(self, thr=3., axis=2, nt_chunk=1024, dsample_nfreq=None, dsample_nt=None, test=False):

        assert (axis == 0 or axis == 1 or axis == 2),\
            "axis must be 0 (along freq; constant time), 1 (along time; constant freq), or 2 (planar; freq and time)."
        assert thr >= 1., "threshold must be >= 1."
        assert nt_chunk > 0

        assert (dsample_nt is None or dsample_nt > 0), "Invalid downsampling number along the time axis!"
        assert (dsample_nfreq is None or dsample_nfreq > 0), "Invalid downsampling number along the freq axis!"

        if test:
            assert (dsample_nfreq is None and dsample_nt is None), "Set dsample_nfreq and dsample_nt to None before running the test!"

        self.thr = thr
        self.axis = axis
        self.nt_chunk = nt_chunk
        self.nt_prepad = 0
        self.nt_postpad = 0
        self.dsample_nfreq = dsample_nfreq
        self.dsample_nt = dsample_nt
        self.test = test

    def set_stream(self, stream):
        if self.dsample_nfreq is None:
            self.dsample_nfreq = stream.nfreq
        if self.dsample_nt is None:
            self.dsample_nt = self.nt_chunk

        if stream.nfreq % self.dsample_nfreq != 0:
            raise RuntimeError("plotter_transform: current implementation requires 'dsample_nfreq' to be a divisor of stream nfreq.")
        if self.nt_chunk % self.dsample_nt != 0:
            raise RuntimeError("clipper_transform: current implementation requires 'dsample_nt' to be a divisor of 'nt_chunk'.")

        self.coarse_grained = (self.dsample_nfreq < stream.nfreq) or (self.dsample_nt < self.nt_chunk)
        self.nfreq = stream.nfreq

    def process_chunk(self, t0, t1, intensity, weights, pp_intensity, pp_weights):
        
        if self.test:
            weights, intensity = self._clipper_transform__test(weights, intensity)

        # Let's make a ref to the original high-resolution weights.
        weights_hres = weights
        
        if self.coarse_grained:
            # Downsample the weights and intensity.
            (intensity, weights) = rf_pipelines.wi_downsample(intensity, weights,\
                    self.dsample_nfreq, self.dsample_nt)

        # The purpose of this next block of code is to compute 'clip', a 2D array of
        # shape (dsample_nfreq, dsample_nt) which contains the estimated rms.

        if self.axis != 2: # 1d mode
            # Compute (sum_i W_i I_i^2) and (sum_i W_i)
            num = np.sum(weights*(intensity)**2, axis=self.axis)
            den = np.sum(weights, axis=self.axis)
            np.putmask(den, den==0., 1.0)     # replace 0.0 by 1.0 to avoid divide-by-zero

            clip = np.sqrt(num/den)      # 1D array
            clip = rf_pipelines.tile_arr(clip, self.axis, self.dsample_nfreq, self.dsample_nt)   # 2D array

        else:
            # 2d mode
            num = np.sum(weights * intensity**2)
            den = np.sum(weights)

            if den == 0.:
                return

            clip = np.zeros((self.dsample_nfreq, self.dsample_nt))
            clip[:] = np.sqrt(num/den)

        assert weights.shape == intensity.shape == clip.shape

        # Boolean array which is True for masked values
        mask = np.abs(intensity) > (self.thr * clip)
        if self.coarse_grained:
            mask = rf_pipelines.upsample(mask, self.nfreq, self.nt_chunk)

        # Assign zero weights to those elements that have an
        # intensity value beyond the threshold limit.
        np.putmask(weights_hres, mask, 0.)

        if self.test: 
            print np.count_nonzero(weights) /\
                float(self.nfreq*self.nt_chunk) * 100, "% not masked."

    def __test(self, weights, intensity):
        # Let's replace the intensity array with gaussian noise
        # centered at 0 with std=1. Also set all weights to 1.
        # Masking elements beyond thr=1 should result in keeping 
        # ~68% of all elements in weights (i.e., non-zero and 
        # within 1std).
        intensity[:] = np.random.normal(0, 1, intensity.size).reshape(intensity.shape)
        weights[:] = 1.
        self.thr = 1.
        return weights, intensity
