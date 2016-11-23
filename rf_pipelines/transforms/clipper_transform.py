import numpy as np
import rf_pipelines

class clipper_transform(rf_pipelines.py_wi_transform):
    """
   This transform clips the intensity along a selected 
   axis -- also works in planar (2D) mode -- and above 
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

      t = clipper_transform(thr=3, axis=None, nt_chunk=1024,\
          dsample_nfreq=None, dsample_nt=None, test=False)

      'thr=3.' is the multiplicative factor of maximum threshold,
        e.g., 3 * standard_deviation, meaning that (the absolute
        value of) any sample above this limit is clipped.

      'axis=None' is the axis convention:
        None: planar; freq and time. 
        0: along freq; constant time.
        1: along time; constant freq.

      'nt_chunk=1024' is the buffer size.

      'dsample_nfreq' and 'dsample_nt' are the downsampled 
       number of pixles along the freq and time axes, respectively.

      'test=False' enables a test mode.
    """

    def __init__(self, thr=3., axis=None, nt_chunk=1024, dsample_nfreq=None, dsample_nt=None, test=False):

        assert (axis == None or axis == 0 or axis == 1),\
            "axis must be None (planar; freq and time), 0 (along freq; constant time), or 1 (along time; constant freq)."
        assert thr >= 1., "threshold must be >= 1."
        assert nt_chunk > 0

        assert (dsample_nt is None or dsample_nt > 0), "Invalid downsampling number along the time axis!"
        assert (dsample_nfreq is None or dsample_nfreq > 0), "Invalid downsampling number along the freq axis!"

        name = 'clipper_transform(thr=%f, axis=%s, nt_chunk=%d' % (thr, axis, nt_chunk)
        if dsample_nfreq is not None:
            name += ', dsample_nfreq=%d' % dsample_nfreq
        if dsample_nt is not None:
            name += ', dsample_nt=%d' % dsample_nt
        name += ')'

        self.thr = thr
        self.axis = axis
        self.name = name
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
        # If 'test' is specified, this will be a pseudo-transform which doesn't modify data
        # in the pipeline, but simulates some fake Gaussian data and reports the clipped fraction
        if self.test:
            intensity = np.random.normal(0, 1, size=intensity.shape)
            weights = np.ones(weights.shape)

        # Let's make a ref to the original high-resolution weights.
        weights_hres = weights
        
        if self.coarse_grained:
            # Downsample the weights and intensity.
            (intensity, weights) = rf_pipelines.wi_downsample(intensity, weights,\
                    self.dsample_nfreq, self.dsample_nt)

        # The purpose of this next block of code is to compute 'clip', a 2D array of
        # shape (dsample_nfreq, dsample_nt) which contains the estimated rms.

        # Compute (sum_i W_i I_i^2) and (sum_i W_i)
        num = np.asarray(np.sum(weights*(intensity)**2, axis=self.axis))
        den = np.asarray(np.sum(weights, axis=self.axis))
        
        np.putmask(den, den==0., 1.0)     # replace 0.0 by 1.0 to avoid divide-by-zero

        clip = np.sqrt(num/den)
        clip = rf_pipelines.tile_arr(clip, self.axis, self.dsample_nfreq, self.dsample_nt)

        assert weights.shape == intensity.shape == clip.shape
        
        if self.axis == None:
            (mean, rms) = rf_pipelines.weighted_mean_and_rms(intensity, weights, 6, 3)
            clip[:] = rms

            # Boolean array which is True for masked values
            mask = np.abs(intensity-mean) > (self.thr * clip)
        
        if self.axis != None:
            mask = np.abs(intensity) > (self.thr * clip)

        if self.coarse_grained:
            mask = rf_pipelines.upsample(mask, self.nfreq, self.nt_chunk)

        # Assign zero weights to those elements that have an
        # intensity value beyond the threshold limit.
        np.putmask(weights_hres, mask, 0.)

        if self.test: 
            unmasked_percentage = np.count_nonzero(weights_hres) / float(weights_hres.size) * 100.
            print unmasked_percentage, "% not masked."
