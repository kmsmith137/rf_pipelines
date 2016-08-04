import numpy as np
import rf_pipelines

class clipper_transform(rf_pipelines.py_wi_transform):
    """
   This transform clips the intensity along a selected 
   axis and above a given threshold. Results are applied 
   to the weights array (i.e., weights[clipped] = 0.) 
   for masking extreme values.

   Currently based on the standard deviation from
   the mean (np.std)

    Constructor syntax:

      t = clipper_transform(thr=3, axis=0, nt_chunk=1024)

      'thr=3.' is the multiplicative factor of maximum threshold,
        e.g., 3 * standard_deviation, meaning that (the absolute
        value of) any sample above this limit is clipped.

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

        assert (self.axis == 0 or self.axis == 1),\
            'axis must be 0 (along freq; constant time) or 1 (along time; constant freq).'
        assert self.thr >= 1., 'threshold must be >= 1.'

    def set_stream(self, stream):

        self.nfreq = stream.nfreq
        
        # This 2d array will be used as a boolean mask
        # for selecting the intensity elements beyond
        # the threshold.
        self.clip = np.zeros([self.nfreq,self.nt_chunk])

    def process_chunk(self, t0, t1, intensity, weights, pp_intensity, pp_weights):

        # This 1d array holds the sum of weights along
        # the selected axis. Its size is equal to the
        # unselected axis.
        self.sum_weights = np.sum(weights, axis=self.axis)

        # This 1d array has to match (in dimensions)
        # with self.sum_weights during the first element-
        # by-element operation (see the next line).
        self.weighted_mean = np.zeros([np.array(self.sum_weights.shape)])

        # The first element-by-element operation (1d).
        # Here we find the weighted mean values
        # (for those that have non-zero weights)
        # along the selected axis. The results has 
        # the size of the unselected axis.
        self.weighted_mean[self.sum_weights > 0.] = \
            np.sum(weights*intensity, axis=self.axis)[self.sum_weights > 0.]/\
            self.sum_weights[self.sum_weights > 0.]

        # Let's tile our 1d arrays by using 
        # self.tile_arr() so that their dimensions
        # (now in 2d; [nfreq, nt_chunk]) match 
        # with intensity, weights, and self.clip.
        self.sum_weights = self.tile_arr(self.sum_weights)
        self.weighted_mean = self.tile_arr(self.weighted_mean) 

        # Here is the second element-by-element operation (2d).
        # Note that np.sum(2d) results in a 1d array. Therefore,
        # we have to use self.tile_arr() to make the 2d elements 
        # one-to-one.
        self.clip[self.sum_weights > 0.] = np.sqrt(self.tile_arr(\
            np.sum(weights*(intensity-self.weighted_mean)**2,\
            axis=self.axis))[self.sum_weights > 0.]/\
            self.sum_weights[self.sum_weights > 0.])

        assert weights.shape == intensity.shape == self.clip.shape

        # Assign zero weights to those elements that have an
        # intensity value above/below the threshold limit.
        # The statement in [ ] creates a boolean array which
        # picks up elements in the weights array.
        weights[ np.abs(intensity) > (self.thr*self.clip) ] = 0.

    def tile_arr(self, arr):
        if self.axis == 0:
            return np.tile(arr, (self.nfreq, 1))
        else:
            return np.transpose(np.tile(arr, (self.nt_chunk, 1)))
