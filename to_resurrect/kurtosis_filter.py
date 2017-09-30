import numpy as np
from scipy.stats import kurtosis
import rf_pipelines

class kurtosis_filter(rf_pipelines.py_wi_transform):
    """
    This transform masks channels on the basis of excess kurtosis.

    For Gaussian data, the excess kurtosis will be zero.
    For data from a chi-squared distribution, the ex.kurtosis is 12/df.
    For single valued data, the ex.kurtosis is returned as -3.

    Negative ex.kurtosis -> a broader than Gaussian distribution (leptokurtic)
    Positive ex.kurtosis -> a more sharply peaked distribution   (platykurtic)

    Constructor syntax:
        
        t = kurtosis_filter(thr=(-1,1),nt_chunk=1024)

        'thr=(-1,1)' gives the acceptable range for data to be unmasked

        'nt_chunk=512' is the buffer size.
    """ 

    def __init__(self,thr=(-1,1),nt_chunk=1024):
        rf_pipelines.py_wi_transform.__init__(self)

        assert (type(thr) is tuple) and (thr[0] < thr[1]), "Bad threshold choice! See docstring."
        self.lo_cut, self.hi_cut = thr 
        self.nt_chunk = nt_chunk

    def set_stream(self,stream):
        self.nfreq = stream.nfreq

    def process_chunk(self, t0, t1, intensity, weights, pp_intensity, pp_weights):
        k = kurtosis(np.ma.masked_where(weights == 0,intensity,copy=False),1)
        k[~np.isfinite(k)] = -4
        weights[(k < self.lo_cut) | (k > self.hi_cut)] = 0.
