import numpy as np
import rf_pipelines

class badchannel_mask(rf_pipelines.py_wi_transform):
    """
   This transform sets bad freq channels of a weights array to 0.

    Constructor syntax:

      t = badchannel_mask(maskpath, freq_lo_MHz=400.0, freq_hi_MHz=800.0, nt_chunk=1024)

      'maskpath' is the full path to a mask file that contains affected freq 
      intervals, written in rows with the following format: e.g., 420.02,423.03

      'freq_lo_MHz=400.0' and 'freq_hi_MHz=800.0' are the lowest and the highest 
      possible feq for the mask, respectively.

      'nt_chunk=1024' is the buffer size, which is allowed to have a different 
      value in a chain of transforms.
    """

    def __init__(self, maskpath, freq_lo_MHz=400.0, freq_hi_MHz=800.0, nt_chunk=1024):
        
        self.maskpath = maskpath
        self.freq_lo_MHz = freq_lo_MHz
        self.freq_hi_MHz = freq_hi_MHz
        
        self.nt_chunk = nt_chunk
        self.nt_prepad = 0
        self.nt_postpad = 0
        
        # Read the frequency mask into a numpy array.
        self.freq_mask = np.genfromtxt(self.maskpath, delimiter=',')
        
        # Now let's make sure all input arguments are valid.
        assert self.freq_hi_MHz > self.freq_lo_MHz, "freq_hi_MHz must be greater than freq_lo_MHz."
        i = 0
        for k in self.freq_mask:
            assert k[0] <= k[1], "The frequency interval %r in '%s' is not a valid input" % (k, self.maskpath)
            # Remove useless intervals.
            if ((k[0] < self.freq_lo_MHz) and (k[1] < self.freq_lo_MHz)) \
            or ((k[0] > self.freq_hi_MHz) and (k[1] > self.freq_hi_MHz)):
                self.freq_mask = np.delete(self.freq_mask, i, 0)
            i += 1

    def set_stream(self, s):

        # self.nfreq corresponds to the number of freq channels in the chunk.
        self.nfreq = s.nfreq

        # First we need to scale our frequency mask so that it matches with
        # the number of channels in the chunk. Then we make a copy of
        # the mask that will be processed and used as an array of indexes.        
        scale = self.nfreq / (self.freq_hi_MHz - self.freq_lo_MHz)
        self.index_mask = self.freq_mask * scale

        # Subtracting the max value from the mask (which runs from low
        # to high values) leaves us with an array that runs in reverse.
        self.index_mask = (self.freq_hi_MHz * scale) - self.index_mask

        # This part is a bit tricky! First we make sure that we include
        # any non-zero overlap with the mask. To this end, we take the
        # ceiling of the first element -- remember this is in a reversed
        # order, which means, e.g., a[0,1] has become a[1,0] -- and the
        # floor of the second element. Next, we convert them to integers
        # so they can be feed as indexes into the weights array.
        self.index_mask[:,0] = (np.ceil(self.index_mask[:,0])).astype(int)
        self.index_mask[:,1] = (np.floor(self.index_mask[:,1])).astype(int)
        
        # This line is to make sure that we don't use negative
        # indexes in the lower bound. Numpy arrays however accept
        # numbers beyond the maximum index -- so no need to constrain
        # the upper bound.
        self.index_mask[self.index_mask < 0.] = int(0)

    def process_chunk(self, t0, t1, intensity, weights, pp_intensity, pp_weights):
        
        # Here we loop over bad frequency intervals. Note that index 
        # values have to be used in the reversed order since we already 
        # subtracted the max value from the mask.
        for k in self.index_mask:
            weights[k[1]:k[0],:] = 0.
