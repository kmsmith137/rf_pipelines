import sys
import numpy as np

try:
    import PIL.Image
except:
    pass  # warning message has already been printed in rf_pipelines/__init__.py


def expand_array(arr, new_shape, axis):
    arr = np.array(arr)
    ret = np.empty(new_shape, dtype=arr.dtype)

    if axis is None:
        assert arr.ndim == 0
        ret.fill(arr)
        return ret

    assert 0 <= axis < len(new_shape)

    expected_shape = new_shape[:axis] + new_shape[(axis+1):]
    assert arr.shape == expected_shape
    
    # adds new length-1 axis
    arr = np.expand_dims(arr, axis)

    # numpy array assignment will correctly "promote" the length-1 axis to length new_shape[axis].
    ret[:] = arr[:]
    return ret


def weighted_mean_and_rms(arr, weights, niter=1, sigma_clip=3.0, axis=None):
    """
    If niter > 1, then the calculation will be iterated, "clipping" outlier samples which
    deviate from the mean by more than 3 sigma (or a different threshold, if the sigma_clip
    parameter is specified).

    If axis=None, returns a pair of scalars (mean, rms).

    If axis is not None, returns a pair of arrays (mean, rms) whose dimension is
    one less than the input arrays.
    """

    assert weights is not None
    assert arr.shape == weights.shape
    assert niter >= 1
    assert sigma_clip >= 1.0   # lower than this really wouldn't make sense
    assert np.all(weights >= 0.0)

    input_shape = arr.shape
    weights = np.copy(weights)   # make copy, since we'll modify the weights as we iterate

    for iter in xrange(niter):
        wsum = np.sum(weights, axis=axis)

        if iter == 0:
            # Some initializations in the first iteration of the loop
            output_shape = wsum.shape
            mean = np.zeros(output_shape)
            rms = np.zeros(output_shape)

        # The mask has the dimensions of the output array.
        # Outputs where the corresponding weight sum is zero (or negative) are masked.
        # Masked outputs will be represented by rms=0 and arbitrary mean.
        mask = (wsum > 0.0)

        # "Regulate" wsum: replaced masked entries by 1.0 to prevent divide-by zero.
        # These entries will be masked in the end anyway, so the value we use is arbitrary; we just
        # want python to process these array elements without throwing a divide-by-zero exception.
        wsum = np.where(mask, wsum, 1.0)

        # Compute the mean, and expand it to have the dimensions of the input array.
        mean = np.sum(weights*arr, axis=axis) / wsum
        mean_e = expand_array(mean, input_shape, axis)

        var = np.sum(weights*(arr-mean_e)**2, axis=axis) / wsum

        # If the variance is very small compared to the mean, then the result isn't 
        # numerically stable, so just mask it by setting the variance to zero.
        mask = np.logical_and(mask, var > 1.0e-10 * mean**2)
        var = np.where(mask, var, 0.0)

        # Compute rms and expand it to have dimensions of the input array
        rms = np.sqrt(var)
        rms_e = expand_array(rms, input_shape, axis)

        # After this line, 'mask' now has dimensions of the input array.
        # Note that the previous mask is implicitly included, since the rms will be zero
        # for masked entries, and the comparison will evaluate to False.
        mask = np.abs(arr-mean_e) < sigma_clip * rms_e

        # Expand the mask and proceed to the next iteration
        weights = np.where(mask, weights, 0.0)

    # Sometimes useful for debugging
    # print >>sys.stderr, 'weighted_mean_and_rms:', (mean,rms)

    return (mean, rms)


def write_png(filename, arr, weights=None, transpose=False, ytop_to_bottom=False, clip_niter=3, sigma_clip=3.0):
    """
    Writes a 2D floating-point array as a png image.  Currently we use a simple blue-purple-red colormap.

       'arr': A 2D array to be plotted

       'weights': If specified, elements with zero/low weight will be black/greyed out.
 
       'transpose': If set, array axis ordering will be (y,x) rather than the default (x,y).

       'ytop_to_bottom': If set, the array y-axis will run from top->bottom in the image, rather than the default bottom->top.

       'clip_niter', 'sigma_clip': By default, colors are assigned by computing the mean and rms after clipping 3-sigma 
           outliers using three masking iterations.  These arguments override the defaults.
    """

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

    # Another corner case: if rms is zero then use 1.0 and fall through.  
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

    (nfreq, ntime) = intensity.shape
    w = w / (nfreq//new_nfreq * ntime//new_ntime)
    return (wi, w)

def upsample(arr, new_nfreq, new_nt):
    """Upsamples a 2d array"""

    (old_nfreq, old_nt) = arr.shape
    assert new_nfreq % old_nfreq == 0
    assert new_nt % old_nt == 0

    (r_nfreq, r_nt) = (new_nfreq // old_nfreq, new_nt // old_nt)
    ret = np.zeros((old_nfreq, r_nfreq, old_nt, r_nt), dtype=arr.dtype)

    for i in xrange(r_nfreq):
        for j in xrange(r_nt):
            ret[:,i,:,j] = arr[:,:]
    
    return np.reshape(ret, (new_nfreq, new_nt))

def tile_arr(arr, axis, nfreq, nt_chunk):
    """tiles (i.e., copies) a scalar or a 1d array to a 2d array.
    It's used for matching 1d and 2d arrays in element-by-element 
    operations. It can also be useful in creating 2d simulations.
    
    Axis convention:
    None: planar; freq and time.
    0: tile along freq; constant time
    1: tile along time; constant freq
    """
    
    assert (arr.ndim == 0 or arr.ndim ==1)
    assert (axis == None) or (axis == 0) or (axis == 1),\
            "axis must be None (planar; freq and time), 0 (along freq; constant time), or 1 (along time; constant freq)."
    if axis == 0:
        return np.tile(arr, (nfreq,1))
    elif axis == 1:
        return np.transpose(np.tile(arr, (nt_chunk,1)))
    else:
        return np.tile(arr, (nfreq,nt_chunk))

####################################################################################################


def json_show(obj, depth=1):
    """
    Prints a partially-expanded summary of json object 'obj' to stdout.
    The 'depth' parameter controls the amount of expansion.
    """

    print json_str(obj, depth, indent='')

    
def json_str(obj, depth=1, indent=''):
    """
    Returns a partially-expanded summary of json object 'obj' as a string.

    The 'depth' parameter has the following meaning:
       depth < 0:   one-word summary (e.g. if obj is a list then 'list' will be returned)
       depth = 0:   one-line summary, long lists/dicts will be abbreviated
       depth = 1:   multi-line summary, all entries of lists/dicts will be shown
       depth > 1:   multi-line summary, sublists/subdicts will be partially expanded to (depth-1).
    """

    if isinstance(obj, basestring):
        return '"%s"' % obj

    if isinstance(obj,int) or isinstance(obj,float) or isinstance(obj,bool):
        return repr(obj)

    if isinstance(obj, list):
        if depth < 0:
            return 'list'

        if depth > 0:
            x = [ '%s    %s\n' % (indent, json_str(t,depth-1,indent+'    ')) for t in obj ]
            return '[\n%s%s]' % (''.join(x), indent)

        if len(obj) > 4:
            return '[ %s, %s, ..., %s ]' % (json_str(obj[0],-1), json_str(obj[1],-1), json_str(obj[-1],-1))

        x = [ json_str(t,-1) for t in obj ]
        return '[ %s ]' % (', '.join(x))

    if isinstance(obj, dict):
        if depth < 0:
            return 'dict'

        if depth > 0:
            x = [ ((isinstance(v,dict) or isinstance(v,list)), k) for (k,v) in obj.iteritems() ]
            x = [ (k,obj[k]) for (t,k) in sorted(x) ]
            x = [ '%s    "%s": %s\n' % (indent, k, json_str(v,depth-1,indent+'    ')) for (k,v) in x ]
            return '{\n%s%s}' % (''.join(x), indent)

        x = sorted(obj.keys())
        
        if len(obj) > 4:
            return '{ "%s":%s, "%s":%s, ..., "%s":%s }' % (x[0], json_str(obj[x[0]],-1), x[1], json_str(obj[x[1]],-1), x[-1], json_str(obj[x[-1]],-1))

        x = [ '"%s":%s' % (k,json_str(obj[k],-1)) for k in x ]
        return '{ %s }' % (', '.join(x))
    
    raise RuntimeError('rf_pipelines.json_str(): unrecognized object')

