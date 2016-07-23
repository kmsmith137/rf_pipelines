import numpy as np
import rf_pipelines

class badchannel_mask(rf_pipelines.py_wi_transform):
    """
   This transform sets bad freq channels of a weights array to 0.

    Constructor syntax:

      t = badchannel_mask(maskpath, nt_chunk=1024)

      'maskpath' is the full path to a mask file that contains affected freq 
      intervals, written in rows with the following format: e.g., 420.02,423.03

      'nt_chunk=1024' is the buffer size, which is allowed to have a different 
      value in a chain of transforms.
    """

    def __init__(self, maskpath, nt_chunk=1024):
        
        self.maskpath = maskpath 
        self.nt_chunk = nt_chunk
        self.nt_prepad = 0
        self.nt_postpad = 0
        
        # Read the frequency mask into a numpy array.
        self.freq_mask = np.genfromtxt(self.maskpath, delimiter=',')

    def set_stream(self, stream):

        # self.nfreq corresponds to the number of freq channels in the chunk.
        self.nfreq = stream.nfreq
        
        # freq_lo_MHz and freq_hi_MHz are the lowest and the highest 
        # possible feq for the mask, respectively.
        self.freq_lo_MHz = stream.freq_lo_MHz
        self.freq_hi_MHz = stream.freq_hi_MHz
        
        # Check whether all freq intervals are valid and whether they overlap 
        # with the freq range of the stream. Keep only the ones that do overlap.
        trimmed_mask = []
        for (freq0,freq1) in self.freq_mask:
            assert freq0 <= freq1, "The frequency interval %r in '%s' is not a valid input" \
                % ([freq0,freq1], self.maskpath)
            if (self.freq_lo_MHz <= freq0 <= self.freq_hi_MHz) \
                or (self.freq_lo_MHz <= freq1 <= self.freq_hi_MHz):
                    trimmed_mask.append([freq0,freq1])
        self.freq_mask = np.array(trimmed_mask)

        # First we need to scale our frequency mask so that it matches with
        # the number of channels in the chunk. Then we make a copy of
        # the mask that will be processed and used as an array of indexes.        
        self.scale = self.nfreq / (self.freq_hi_MHz - self.freq_lo_MHz)
        self.index_mask = self.freq_mask * self.scale

        # Subtracting the max value from the mask (which runs from low
        # to high values) leaves us with an array that runs in reverse.
        self.index_mask = (self.freq_hi_MHz * self.scale) - self.index_mask

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

        # The following is a private method for testing the class.
        self._badchannel_mask__test()

    def __test(self):
        pass

    def process_chunk(self, t0, t1, intensity, weights, pp_intensity, pp_weights):
        
        # Here we loop over bad frequency intervals. Note that index 
        # values have to be used in the reversed order since we already 
        # subtracted the max value from the mask.
        for (freq1,freq0) in self.index_mask:
            weights[freq0:freq1,:] = 0.
