import numpy as np
import rf_pipelines

class clipper_transform(rf_pipelines.py_wi_transform):
    """
   Clips the intensity along a selected axis and 
   above a given threshold. Results are applied to 
   the weights array (i.e., weights[clipped] = 0.) 
   for masking extreme values.

   Currently based on the standard deviation from
   the mean (np.std)

    Constructor syntax:

      t = clipper_transform(thr=3, axis=0, nt_chunk=1024)

      'thr=3.' is the multiplicative factor of maximum threshold,
        e.g., 3 * standard_deviation, meaning that any sample
        above/below this limit is clipped.

      'axis=0' is the axis convention:
        0: along freq; constant time
        1: along time; constant freq

      'nt_chunk=1024' is the buffer size.
    """

    def __init__(self, thr=3., axis=0, nt_chunk=1024):
        
        self.thr = thr
        self.axis = axis
        self.nt_chunk = nt_chunk
        self.nt_prepad = 0
        self.nt_postpad = 0

        assert (self.axis == 0 or self.axis == 1), \
        'axis must be 0 (along freq; constant time) or 1 (along time; constant freq).'
        assert self.thr >= 1., 'threshold must be >= 1.'

    def set_stream(self, stream):

        self.nfreq = stream.nfreq
        
    def process_chunk(self, t0, t1, intensity, weights, pp_intensity, pp_weights):
        
        # Compute the standard deviation of the intensity 
        # array (2d) along the selected axis. The result is 
        # a 1d array with a size equivalent to the number of 
        # elements along the other axis.
        self.clip = np.std(intensity, axis=self.axis)

        # The idea here is to dimensionally match the self.clip 
        # array (1d) with the intensity (2d) and weights (2d)
        # so that their elements pick up a one-to-one logical 
        # correspondance. To this end, we tile (in-place) self.clip 
        # by copying the 1d array along the selected axis.
        if self.axis == 0:
            self.clip = np.tile(self.clip, (self.nfreq, 1))
        else:
            self.clip = np.transpose(np.tile(self.clip, (self.nt_chunk, 1)))
        
        # Verify that all the arrays have the same dimensions.
        assert np.shape(weights) == np.shape(intensity) == np.shape(self.clip)

        # Assign zero weights to those elements that have an
        # intensity value above/below the threshold limit.
        # The statement in [ ] creates a boolean array which
        # picks up elements in the weights array.
        weights[ np.abs(intensity) > (self.thr*self.clip) ] = 0.
