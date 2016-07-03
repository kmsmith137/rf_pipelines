import sys
import numpy as np
import rf_pipelines

try:
    import PIL.Image
except ImportError:
    print >>sys.stderr, "rf_pipelines: import PIL.Image failed; many things will work but plotting will fail"


def write_png(filename, arr, weights=None, transpose=False, ytop_to_bottom=False):
    arr = np.array(arr, dtype=np.float)
    assert arr.ndim == 2

    if weights is not None:
        weights = np.array(weights, dtype=np.float)
        assert weights.shape == arr.shape
    
    if not transpose:
        arr = np.transpose(arr)
        weights = np.transpose(weights) if (weights is not None) else None
    
    if not ytop_to_bottom:
        arr = arr[::-1]
        weights = weights[::-1] if (weights is not None) else None

    (mean, rms) = (np.mean(arr), np.var(arr)**0.5)

    # color in range [0,1].
    color = 0.5 + 0.16*(arr-mean)/rms    # factor 0.16 preserves convention from some old code
    color = np.maximum(color, 0.0)
    color = np.minimum(color, 0.999999)  # 0.99999 instead of 1.0, to make roundoff-robust
    
    # rgb in range [0,1]
    red = 256. * color
    blue = 256. * (1-color)

    if weights is not None:
        (wmin, wmax) = (np.min(weights), np.max(weights))
        if wmin < 0:
            raise RuntimeError('write_png: negative weights are currently treated as an error')

        # avoid NaN
        if wmax == 0.0:
            wmax = 1.0

        weights = weights/wmax
        red *= weights
        blue *= weights

    rgb = np.zeros((arr.shape[0],arr.shape[1],3), dtype=np.uint8)
    rgb[:,:,0] = red
    rgb[:,:,2] = blue

    img = PIL.Image.fromarray(rgb)
    img.save(filename)
    print >>sys.stderr, 'wrote %s' % filename


def _downsample_2d(arr, new_nfreq, new_ntime):
    """Helper for wi_downsample."""

    assert arr.ndim == 2
    assert new_nfreq > 0
    assert new_ntime > 0

    (nfreq, ntime) = arr.shape
    assert nfreq % new_nfreq == 0
    assert ntime % new_ntime == 0

    arr = np.reshape(arr, (new_nfreq, nfreq//new_nfreq, new_ntime, ntime//new_ntime))
    arr = np.sum(arr, axis=3)
    arr = np.sum(arr, axis=1)
    return arr
    

def wi_downsample(intensity, weights, new_nfreq, new_ntime):
    """Downsamples a pair of 2D arrays (intensity, weights), returning a new pair (ds_intensity, ds_weights)."""

    wi = _downsample_2d(weights * intensity, new_nfreq, new_ntime)
    w = _downsample_2d(weights, new_nfreq, new_ntime)

    mask = (w > 0.0)
    wi = wi / np.where(mask, w, 1.0)
    wi = np.where(mask, wi, 0.0)

    return (wi, w)


class plotter_transform(rf_pipelines.wi_transform):
    def __init__(self, img_prefix, img_nfreq, img_nt, downsample_nt=1, nt_chunk=0):
        assert img_nt > 0
        assert img_nfreq > 0
        assert downsample_nt > 0
        
        if nt_chunk == 0:
            nt_chunk = min(img_nt, 1024//downsample_nt+1) * downsample_nt    # default value

        assert nt_chunk > 0

        if nt_chunk % downsample_nt != 0:
            raise RuntimeError("plotter_transform: current implementation requires 'nt_chunk' to be a multiple of 'downsample_nt'")

        self.img_prefix = img_prefix
        self.img_nfreq = img_nfreq
        self.img_nt = img_nt
        self.nt_chunk_ds = nt_chunk // downsample_nt

        self.nt_chunk = nt_chunk
        self.nt_prepad = 0
        self.nt_postpad = 0


    def set_stream(self, s):
        self.nfreq = s.nfreq

        if s.nfreq % self.img_nfreq != 0:
            raise RuntimeError("plotter_transform: current implementation requires 'img_nfreq' to be a divisor of stream nfreq")


    def start_substream(self, isubstream, t0):
        self.intensity_buf = np.zeros((self.img_nfreq,self.img_nt), dtype=np.float32)
        self.weight_buf = np.zeros((self.img_nfreq,self.img_nt), dtype=np.float32)
        self.isubstream = isubstream
        self.ifile = 0
        self.ipos = 0


    def process_chunk(self, t0, t1, intensity, weights, pp_intensity, pp_weights):
        (intensity, weights) = wi_downsample(intensity, weights, self.img_nfreq, self.nt_chunk_ds)

        ichunk = 0
        while ichunk < self.nt_chunk_ds:
            n = min(self.nt_chunk_ds - ichunk, self.img_nt - self.ipos)
            
            self.intensity_buf[:,self.ipos:(self.ipos+n)] = intensity[:,ichunk:(ichunk+n)]
            self.weight_buf[:,self.ipos:(self.ipos+n)] = weights[:,ichunk:(ichunk+n)]
            self.ipos += n
            ichunk += n

            if self.ipos == self.img_nt:
                self._write_file()


    def end_substream(self):
        if self.ipos > 0:
            self._write_file()


    def _write_file(self):
        filename = self.img_prefix
        if self.isubstream > 0:
            filename += str(isubstream+1)
        filename += ('_%s.png' % self.ifile)

        write_png(filename, self.intensity_buf[:,:self.ipos], weights=self.weight_buf[:,:self.ipos], transpose=True, ytop_to_bottom=True)

        self.ifile += 1
        self.ipos = 0