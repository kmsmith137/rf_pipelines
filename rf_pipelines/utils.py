import sys
import numpy as np

try:
    import PIL.Image
except:
    pass  # warning message has already been printed in rf_pipelines/__init__.py


def weighted_mean_and_rms(arr, weights, niter=1, sigma_clip=3.0):
    """
    If niter > 1, then the calculation will be iterated, "clipping" outlier samples which
    deviate from the mean by more than 3 sigma (or a different threshold, if the sigma_clip
    parameter is specified).
    """

    assert weights is not None
    assert arr.shape == weights.shape
    assert niter >= 1
    assert sigma_clip >= 2.0   # lower than this really wouldn't make sense

    assert np.all(weights >= 0.0)
    weights = np.copy(weights)

    (mean, rms) = (0.0, 0.0)

    for iter in xrange(niter):
        wsum = np.sum(weights)
        if wsum <= 0.0:
            return (mean, 0.0)

        mean = np.sum(weights*arr) / np.sum(weights)
        var = np.sum((weights*(arr-mean))**2) / np.sum(weights**2)

        if var <= 0.0:
            return (mean, rms)

        rms = np.sqrt(var)
        mask = np.abs(arr-mean) <= sigma_clip*rms
        weights = np.where(mask, weights, 0.0)

    return (mean, rms)


def write_png(filename, arr, weights=None, transpose=False, ytop_to_bottom=False, clip_niter=3, sigma_clip=3.0):
    """This is a quick-and-dirty plotting routine that I cut-and-paste everywhere."""

    arr = np.array(arr, dtype=np.float)
    assert arr.ndim == 2

    if weights is None:
        weights = np.ones(arr.shape, dtype=np.float)
    else:
        weights = np.array(weights, dtype=np.float)
        assert weights.shape == arr.shape
    
    if not transpose:
        arr = np.transpose(arr)
        weights = np.transpose(weights) if (weights is not None) else None
    
    if not ytop_to_bottom:
        arr = arr[::-1]
        weights = weights[::-1] if (weights is not None) else None

    (wmin, wmax) = (np.min(weights), np.max(weights))
    if wmin < 0:
        raise RuntimeError('write_png: negative weights are currently treated as an error')

    # A corner case..
    if wmax == 0.0:
        print >>sys.stderr, '%s: array was completely masked, writing all-black image' % filename
        rgb = np.zeros((arr.shape[0], arr.shape[1], 3), dtype=np.uint8)
        img = PIL.Image.fromarray(rgb)
        img.save(filename)
        return

    (mean, rms) = weighted_mean_and_rms(arr, weights, clip_niter, sigma_clip)

    # Anorther corner case: if rms is zero then use 1.0 and fall through.  
    # This will plot an image with constant color values.
    if rms <= 0.0:
        rms = 1.0

    # normalize weights to [0,1]
    weights = weights/wmax

    # color in range [0,1].
    color = 0.5 + 0.16*(arr-mean)/rms    # factor 0.16 preserves convention from some old code
    color = np.maximum(color, 0.0001)    # 0.0001 instead of 0.0, to make roundoff-robust
    color = np.minimum(color, 0.9999)    # 0.9999 instead of 1.0, to make roundoff-robust
    
    # rgb in range [0,1]
    red = 256. * color * weights
    blue = 256. * (1-color) * weights

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
