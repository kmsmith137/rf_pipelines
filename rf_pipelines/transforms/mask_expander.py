import numpy as np
import rf_pipelines

class mask_expander(rf_pipelines.py_wi_transform):
    """
   This transform expands the mask by filling the weights 
   array, provided that its mean is less than or equal to 
   a given threshold, with zeros along a selected axis.

   + Assumes that weights are between 0 and 1.

    Constructor syntax:

      t = mask_expander(thr=0.2, axis=None, nt_chunk=1024, test=False)
      
      'thr=0.2' should be between 0 and 1. Any selected
        weights sub-array with a mean value less than or 
        equal to the threshold is totally masked by this 
        transform.

      'axis=None' is the axis convention:
        None: planar; freq and time.
        0: along freq; constant time.
        1: along time; constant freq.
      
      'nt_chunk=1024' is the buffer size.

      'test=False' enables a test mode.
    """

    def __init__(self, thr=0.2, axis=None, nt_chunk=1024, test=False):
        
        assert (0 < thr < 1), "threshold must be between 0 and 1."
        assert (axis == None) or (axis == 0) or (axis == 1),\
                "axis must be None (planar; freq and time), 0 (along freq; constant time), or 1 (along time; constant freq)."

        self.thr = thr
        self.axis = axis
        self.nt_chunk = nt_chunk
        self.test = test

    def set_stream(self, stream):
        
        self.nfreq = stream.nfreq

        if self.test:
            # 1D mode: 'self.nmax' specifies the max-1 number 
            # of elements along the unselected axis.
            if self.axis == 0:
                self.nmax = self.nt_chunk
            elif self.axis == 1:
                self.nmax = self.nfreq
            # 2D mode: 
            # 'self.stats' += [# of processed chunks, fraction of unmasked weights]
            else:
                self.stats = [1, 0]

    def process_chunk(self, t0, t1, intensity, weights, pp_intensity, pp_weights):
        # If 'test' is specified, this will be a pseudo-transform which doesn't modify weights
        # in the pipeline, but simulates some fake weights and reports the masked fraction.
        if self.test:
            # 1D mode: replace the weights with equally-spaced
            # (and linearly increasing) numbers between 0 and 1.
            if self.axis != None:
                weights = np.linspace(0., 1., self.nmax)
                weights = rf_pipelines.tile_arr(weights, self.axis, self.nfreq, self.nt_chunk)
            # 2D mode: replace the weights with some Gaussian random
            # numbers with (mean, stdv) = (self.thr, 0.01).
            else:
                weights = np.random.normal(self.thr, 0.01, size=weights.shape)

        # Let's make a ref to the weights.
        weights_origin = weights

        # Compute the mean of weights along the seleted axis.
        w_mean = np.mean(weights, axis=self.axis)
        if self.axis == None:
            w_mean = np.array([w_mean])

        # Based on the given threshold, replace all elements 
        # with 0 and 1.
        np.putmask(w_mean, w_mean <= self.thr, 0.) # Target elements (to be tiled in 1D mode)
        np.putmask(w_mean, w_mean > self.thr, 1.)
        
        # 1D mode: Tile the array so that it matches (in dimension) with the original weights.
        if self.axis != None:
            w_mean = rf_pipelines.tile_arr(w_mean, self.axis, self.nfreq, self.nt_chunk)
        
        # Expand the mask by zeroing out the target elements in the weights.
        weights_origin[:] = weights_origin[:] * w_mean[:]

        if self.test:
            unmasked_fraction = np.count_nonzero(weights_origin) / float(weights_origin.size)
            if self.axis != None:
                print unmasked_fraction * 100., "% of the current weights are not masked."
            else:
                self.stats[0] += 1 
                self.stats[1] += unmasked_fraction
                print "===================================="
                print "[# of processed chunks, fraction of unmasked weights] =", self.stats
                print "This must approach 50%; try smaller nt_chunk for better statistics:"
                print self.stats[1]/self.stats[0] * 100., "% of the processed chunks have been unmasked"
                print "===================================="
