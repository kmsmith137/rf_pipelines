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
        
        # Reading the mask file
        self.mask = np.genfromtxt(self.maskpath, delimiter=',')

    def set_stream(self, s):
        
        self.nfreq = s.nfreq
        
        # Scaling bad frequencies according to the stream size
        scale = self.nfreq / (self.freq_hi_MHz - self.freq_lo_MHz)
        self.mask = self.mask * scale

        # Computing extreme index values for the scaled frequencies
        fmax = self.freq_hi_MHz * scale
        fmin = self.freq_lo_MHz * scale

        # Constraining the mask within the lowest and the highest index
        self.mask[self.mask < fmin] = fmin
        self.mask[self.mask > fmax] = fmax - 1. # Subtracting 1 for indexing (see below)

        # The following 3 lines allow us to use the mask values as 
        # array indexes that run from the highest to the lowest freq
        self.mask = fmax - 1. - self.mask
        self.mask[:,0] = (np.ceil(self.mask[:,0])).astype(int)
        self.mask[:,1] = (np.floor(self.mask[:,1])).astype(int)

    def process_chunk(self, t0, t1, intensity, weights, pp_intensity, pp_weights):

        # Looping over bad freq intervals
        for k in self.mask:
            (ifreq, jfreq) = (k[1], k[0])
            if ifreq != jfreq:
                weights[ifreq:jfreq,:] = 0.
            else:
                weights[ifreq,:] = 0. 
