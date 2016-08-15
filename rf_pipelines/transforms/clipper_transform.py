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
   'course_grained', 'dsample_nfreq', and 'dsample_nt'.
    
    Constructor syntax:

      t = clipper_transform(thr=3, axis=0, nt_chunk=1024,\
          dsample_nfreq=512, dsample_nt=512,\
          course_grained=False, test=False)

      'thr=3.' is the multiplicative factor of maximum threshold,
        e.g., 3 * standard_deviation, meaning that (the absolute
        value of) any sample above this limit is clipped.

      'axis=0' is the axis convention:
        0: along freq; constant time.
        1: along time; constant freq.
        2: planar; freq and time.

      'nt_chunk=1024' is the buffer size.

      'dsample_nfreq=512' and 'dsample_nt=512' are the downsampled 
       number of pixles along the freq and time axes, respectively.

      'coarse_grained=False' enables the coarse-grained clipper by
      downsampling the arrays. The weights array is upsampled after
      the clipping process.

      'test=False' enables a test mode.
    """

    def __init__(self, thr=3., axis=0, nt_chunk=1024, dsample_nfreq=None, dsample_nt=None, test=False):

        assert (axis == 0 or axis == 1 or axis == 2),\
            "axis must be 0 (along freq; constant time), 1 (along time; constant freq), or 2 (planar; freq and time)."
        assert thr >= 1., "threshold must be >= 1."
        assert nt_chunk > 0

        assert (dsample_nt is None or dsample_nt > 0), "Invalid downsampling number along the time axis!"
        assert (dsample_nfreq is None or dsample_nfreq > 0), "Invalid downsampling number along the freq axis!"

        if test:
            assert coarse_grained != True, "Set coarse_grained to False before running the test!"

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

        # This 2d array will be used as a boolean mask
        # for selecting the intensity elements beyond
        # the threshold.
        self.clip = np.zeros([self.dsample_nfreq, self.dsample_nt])

    def process_chunk(self, t0, t1, intensity, weights, pp_intensity, pp_weights):
        
        if self.test:
            weights, intensity = self._clipper_transform__test(weights, intensity)

        # Let's make a ref to the original high-resolution weights.
        weights_hres = weights
        
        if self.coarse_grained:
            # Downsample the weights and intensity.
            (intensity, weights) = rf_pipelines.wi_downsample(intensity, weights,\
                    self.dsample_nfreq, self.dsample_nt)

        if self.axis != 2: # 1d mode
            # This 1d array holds the sum of weights along
            # the selected axis. Its size is equal to the
            # unselected axis.  We then call rf_pipelines.tile_arr()
            # so that its dimensions match the original.
            den = np.sum(weights, axis=self.axis)
            den = rf_pipelines.tile_arr(den, self.axis, self.dsample_nfreq, self.dsample_nt)

            # Here is an element-by-element operation (2d).
            # Note that np.sum(2d) results in a 1d array. Therefore,
            # we have to use rf_pipelines.tile_arr() to make the 2d elements 
            # one-to-one.
            indy, indx = np.where(den > 0.)

            # Compute sum_i W_i I_i^2, and use tile_arr() to get an array of the original dimensions
            num = np.sum(weights*(intensity)**2, axis=self.axis)
            num = rf_pipelines.tile_arr(num, self.axis, self.dsample_nfreq, self.dsample_nt)

            self.clip[indy,indx] = np.sqrt(num[indy,indx]/den[indy,indx])

        else:
            # 2d mode
            num = np.sum(weights * intensity**2)
            den = np.sum(weights)

            if den == 0.:
                return

            self.clip[:] = np.sqrt(num/den)

        assert weights.shape == intensity.shape == self.clip.shape

        # Boolean array which is True for masked values
        mask = np.abs(intensity) > (self.thr*self.clip)
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
