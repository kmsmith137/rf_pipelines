import sys
import numpy as np

try:
    import PIL.Image
except:
    print >>sys.stderr, "rf_pipelines: import PIL.Image failed; many things will work but plotting will fail"

from rf_pipelines_c import wi_stream, wi_transform, make_psrfits_stream, make_simple_detrender


####################################################################################################


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
        (wmin, wmax) = (np.minimum(weights), np.maximum(weights))
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
