import sys
import numpy as np
import rf_pipelines

class rfi_transform(rf_pipelines.py_wi_transform):
    """
    a simple transform for removing rfi

    Constructor syntax:

      t = rfi_transform(filepath, freq_lo_MHz=400.0, freq_hi_MHz=800.0, ntchunk=0)

      'filepath' is the full path to an rfi mask file that contains affected freq 
      intervals, written in rows with the following format: e.g., 420.02,423.03 
    """

    def __init__(self, filepath, freq_lo_MHz=400.0, freq_hi_MHz=800.0, nt_chunk=1024):
        self.filepath = filepath
        self.freq_lo_MHz = freq_lo_MHz
        self.freq_hi_MHz = freq_hi_MHz
        self.nt_chunk = nt_chunk
        self.nt_prepad = 0
        self.nt_postpad = 0
        
    def set_stream(self, s):
        self.nfreq = s.nfreq

    def process_chunk(self, t0, t1, intensity, weights, pp_intensity, pp_weights):

        rfi = np.genfromtxt(self.filepath, delimiter=',')
        scale = self.nfreq / (self.freq_hi_MHz - self.freq_lo_MHz)
        
        rfi[rfi < self.freq_lo_MHz] = self.freq_lo_MHz
        rfi[rfi > self.freq_hi_MHz] = self.freq_hi_MHz
        rfi = (self.freq_hi_MHz - rfi) * scale

        rfi[:,0] = (np.ceil(rfi[:,0])).astype(int)
        rfi[:,1] = (np.floor(rfi[:,1])).astype(int)

        for k in rfi:
            (ifreq, jfreq) = (k[1], k[0])
            if ifreq != jfreq:
                weights[ifreq:jfreq,:] = 0.
            else:
                weights[ifreq,:] = 0. 
