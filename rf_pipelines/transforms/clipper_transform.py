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
   'upsample_nfreq' and 'upsample_nt'
    
    Constructor syntax:

      t = clipper_transform(thr=3, axis=0, nt_chunk=1024,\ 
          upsample_nfreq=1, upsample_nt=1, test=False)

      'thr=3.' is the multiplicative factor of maximum threshold,
        e.g., 3 * standard_deviation, meaning that (the absolute
        value of) any sample above this limit is clipped.

      'axis=0' is the axis convention:
        0: along freq; constant time
        1: along time; constant freq
        2: planar; freq and time

      'nt_chunk=1024' is the buffer size.

      'upsample_nfreq' and 'upsample_nt=1' are the upsampling 
      factors along the freq and time axes, respectively.

      'test=False' enables a test mode.
    """

    def __init__(self, thr=3., axis=0, nt_chunk=1024, upsample_nfreq=1, upsample_nt=1, test=False):

        assert (axis == 0 or axis == 1 or axis == 2),\
            "axis must be 0 (along freq; constant time), 1 (along time; constant freq), or 2 (planar; freq and time)"
        assert thr >= 1., "threshold must be >= 1."
        assert upsample_nt > 0, "invalid upsampling factor"
        assert upsample_nfreq > 0, "invalid upsampling factor"

        if upsample_nt % nt_chunk != 0:
            raise RuntimeError("clipper_transform: current implementation requires 'upsample_nt' to be a multiple of 'nt_chunk'.")

        self.thr = thr
        self.axis = axis
        self.nt_chunk = nt_chunk
        self.nt_prepad = 0
        self.nt_postpad = 0
        self.upsample_nfreq = upsample_nfreq
        self.upsample_nt = upsample_nt
        self.test = test

    def set_stream(self, stream):
 
        if self.upsample_nfreq % stream.nfreq != 0:
                raise RuntimeError("plotter_transform: current implementation requires 'upsample_nfreq' to be a multiple of stream nfreq")

        self.nfreq = stream.nfreq

        # This 2d array will be used as a boolean mask
        # for selecting the intensity elements beyond
        # the threshold.
        self.clip = np.zeros([self.nfreq,self.nt_chunk])

    def process_chunk(self, t0, t1, intensity, weights, pp_intensity, pp_weights):
        
        if self.test:
            weights, intensity = self._clipper_transform__test(weights, intensity)
        
        if self.axis != 2: # 1d mode
            # This 1d array holds the sum of weights along
            # the selected axis. Its size is equal to the
            # unselected axis.
            self.sum_weights = np.sum(weights, axis=self.axis)

            # Let's tile our 1d array by using 
            # rf_pipelines.tile_arr() so that its dimensions
            # (now in 2d; [nfreq, nt_chunk]) match 
            # with intensity, weights, and self.clip.
            self.sum_weights = rf_pipelines.tile_arr(self.sum_weights, self.axis,\
                    self.nfreq, self.nt_chunk)
        
            # Here is an element-by-element operation (2d).
            # Note that np.sum(2d) results in a 1d array. Therefore,
            # we have to use rf_pipelines.tile_arr() to make the 2d elements 
            # one-to-one.
            indx, indy = np.where(self.sum_weights > 0.)
            self.clip[indx,indy] = np.sqrt(rf_pipelines.tile_arr(\
                np.sum(weights*(intensity)**2,\
                axis=self.axis), self.axis, self.nfreq,\
                self.nt_chunk)[indx,indy]/\
                self.sum_weights[indx,indy])
        else:
            # 2d mode
            self.sum_weights = np.sum(weights)
            if self.sum_weights != 0.:
                self.clip[:] = np.sqrt(np.sum(weights*(intensity)**2)/self.sum_weights)

        # Assign zero weights to those elements that have an
        # intensity value beyond the threshold limit.
        assert weights.shape == intensity.shape == self.clip.shape
        np.putmask(weights, np.abs(intensity) > (self.thr*self.clip), 0.)
        
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
