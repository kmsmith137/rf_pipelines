import os
import sys
import time
import json
import glob
import h5py
import traceback
import subprocess
import numpy as np
import rf_pipelines_c

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
    If axis=None, returns a pair of scalars (mean, rms).

    If axis is not None, returns a pair of arrays (mean, rms) whose dimension is
    one less than the input arrays.

    If niter > 1, then the calculation will be iterated, "clipping" outlier samples which
    deviate from the mean by more than 3 sigma (or a different threshold, if the sigma_clip
    parameter is specified).

    Sometimes, too much data has been masked to compute the weighted mean and rms.
    In this case, the output arrays will contain elements with rms=0 and arbitrary mean.
    """

    assert weights is not None
    assert arr.shape == weights.shape
    assert niter >= 1
    assert sigma_clip >= 1.0   # lower than this really wouldn't make sense
    assert np.all(weights >= 0.0)

    for iter in xrange(niter):
        if iter > 0:
            # The 'mean' and 'rms' arrays have been computed in a previous iteration of the loop.
            mean_e = expand_array(mean, arr.shape, axis)
            rms_e = expand_array(rms, arr.shape, axis)
            mask = (np.abs(arr-mean_e) < sigma_clip * rms_e)
            weights = weights * mask     # make copy (using the *= operator here would be a bug)

        wsum = np.sum(weights, axis=axis)

        # Start building up a mask, which keeps track of locations where there is not enough
        # data to compute a weighted mean/rms.  The mask has shape 'output_shape'.
        mask = (wsum > 0.0)
        wsum = np.where(mask, wsum, 1.0)    # avoids divide-by-zero

        # Mean (array of shape 'output_shape')
        mean = np.sum(weights*arr, axis=axis) / wsum

        # Variance (array of shape 'output_shape')
        mean_e = expand_array(mean, arr.shape, axis)
        var = np.sum(weights*(arr-mean_e)**2, axis=axis) / wsum

        # Expand the mask to include entries where the variance is very small
        # compared to the mean (in these locations, the variance calculation is
        # numerically unstable).
        mask = np.logical_and(mask, var > 1.0e-10 * mean**2)

        # The multiplication ensures that rms=0 for masked entries
        rms = np.sqrt(mask * var)

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

    # Note: PIL's conventions are reversed relative to ours in both cases:
    #   - PIL axis ordering is (y,x)
    #   - PIL y-axis direction is top-to-bottom
    #
    # so both instances of "not" below are intentional!
    
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
        print '%s: array was completely masked, writing all-black image' % filename
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
    print 'wrote %s' % filename


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
    """
    Downsamples a pair of 2D arrays (intensity, weights), returning a new pair (ds_intensity, ds_weights).

    This python function is morally equivalent to the C++ function rf_pipelines.wi_downsample(),
    but note that the python version takes arguments (new_nfreq, new_nt), whereas the C++ version
    takes (Df, Dt), and the normalization of the weights also differs by a factor (Df*Dt).
    """

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



def _json_compare(j1, j2, name1=None, name2=None):
    """
    Helper function for json_assert_equal().
    Checks recursively whether json objects 'j1' and 'j2' are equal.

    The return value is a quadruple (b, s, j1, j2, name1, name2), where:
        b = True if equal, False if unequal (boolean)
        s = One-line string diagnosing reason for inequality (empty string if b==False)
        (j1,name1) = Mismatched argument 1 (None if b==False)
        (j2,name2) = Mismatched argument 2 (None if b==False)
    """

    if name1 is None:
        name1 = 'argument1'
    if name2 is None:
        name2 = 'argument2'

    if j1.__class__ != j2.__class__:
        s = 'type(%s)=%s, type(%s)=%s' % (name1, j1.__class__, name2, j2.__class__)
        
    if isinstance(j1, dict):
        k1 = set(j1.keys())
        k2 = set(j2.keys())
        k12 = k1.difference(k2)
        k21 = k2.difference(k1)

        if (len(k12) > 0):
            s = ', '.join(["'%s'" % x for x in k12])
            s = '%s contains key(s) %s which are absent in %s' % (name1, s, name2)
            return (False, s, j1, j2, name1, name2)

        if (len(k21) > 0):
            s = ', '.join(["'%s'" % x for x in k21])
            s = '%s contains key(s) %s which are absent in %s' % (name2, s, name1)
            return (False, s, j1, j2, name1, name2)

        for k in k1:
            n1 = '%s[%s]' % (name1, k)
            n2 = '%s[%s]' % (name2, k)
            t = _json_compare(j1[k], j2[k], n1, n2)
            if not t[0]:
                return t

    elif isinstance(j1, list):
        for i in xrange(len(j1)):
            n1 = '%s[%d]'% (name1, i)
            n2 = '%s[%d]'% (name2, i)
            t = _json_compare(j1[i], j2[i], n1, n2)
            if not t[0]:
                return t

    elif (j1 != j2):
        s = '%s and %s are unequal' % (name1, name2)
        return (False, s, j1, j2, name1, name2)

    return (True, '', None, None, None, None)


def json_assert_equal(j1, j2, name1=None, name2=None, verbose=True):
    """Checks that json values j1,j2 are equal, in a verbose way."""

    (b, s, j1, j2, name1, name2) = _json_compare(j1, j2, name1, name2)
    
    if b:
        return

    if verbose:
        print '%s = %s' % (name1, json_str(j1))
        print '%s = %s' % (name2, json_str(j2))
        s += ', see diagnostic print-statements above'

    raise RuntimeError(s)


def json_read(filename):
    """
    Returns a json value.

    I use this wrapper because the built-in python json.load() doesn't give a
    useful error message if it fails.  Using the wrapper, at least we get a
    filename!  (FIXME: I think switching to the 'simplejson' library will help.)
    """

    f = open(filename, 'r')

    try:
        return json.load(f)
    except:
        raise RuntimeError("%s: couldn't parse json file" % filename)


def json_write(filename, p, clobber=False, verbose=True):
    """
    This helper function is sometimes used to write json files.
    The argument 'p' must be an object of class rf_pipelines.pipeline_object.
    
    If 'filename' already exists, and clobber=True, it will be overwritten.

    If 'filename' already exists, and clobber=False, then we check that
    the json content of the file agrees with 'p'.
    """

    assert isinstance(filename, basestring)
    assert isinstance(p, rf_pipelines_c.pipeline_object)
    assert filename.endswith('.json')
    
    d = os.path.dirname(filename)
    if (len(d) > 0) and (not os.path.exists(d)):
        raise RuntimeError("Directory '%s' does not exist" % d)

    j = p.jsonize()

    if clobber or not os.path.exists(filename):
        f = open(filename, 'w')
        json.dump(j, f, indent=4)
        print >>f, ''  # extra newline
        del f          # close file

        if verbose:
            print 'wrote %s' % filename

        return

    # If we get here, then the file exists, and clobber=False.
    # In this case, we check that the json content of the file agrees with 'p'.

    jfile = json.load(open(filename))
    (is_equal, s, j1, j2, name1, name2) = _json_compare(j, jfile, 'new_json', 'old_json')

    if is_equal:
        if verbose:
            print '%s: not updated' % filename
        return

    if verbose:
        print '%s = %s' % (name1, json_str(j1))
        print '%s = %s' % (name2, json_str(j2))
        s += ', see diagnostic print-statements above'

    raise RuntimeError("'%s' already exists, and %s" % (filename, s))



####################################################################################################


def var_comparison_png(name, arr, min=0.5, max=2):
    """
    Plots the ratio of variance files. Variances between min and max are on a greenscale. Lower 
    than min is saturated blue and higher than max is saturated red. Zero variance is white. 
    """
    low = np.where(var < min)
    high = np.where(var > max)
    zero = np.where(var == 0.)
    factor = 255. / high
    
    rgb = np.zeros((var.shape[0], var.shape[1], 3), dtype=np.uint8)
    rgb[:,:,1] = var * green_factor     # Set green values
    rgb[low[0], low[1], 2] = 255.       # Blue
    rgb[low[0], low[1], 1] = 0.
    rgb[high[0], high[1], 0] = 255.     # Red (high values)
    rgb[high[0], high[1], 1:3] = 0.
    rgb[zero[0], zero[1], :] = 255.     # Finally, if anything is 0, make white

    img = PIL.Image.fromarray(rgb)
    img.save(name)
    print 'wrote %s' % name    


def var_to_png(name, var, max=0.01):
    """
    Plots variance h5 files. Variances below max are on a greenscale and above are saturated red. 
    Zero variance is white. 
    """
    high = np.where(var > max)
    zero = np.where(var == 0.)
    factor = 255. / max

    rgb = np.zeros((var.shape[0], var.shape[1], 3), dtype=np.uint8)
    rgb[:,:,1] = var * factor          # Set green values
    rgb[high[0], high[1], 0] = 200.    # Then red (high values)
    rgb[high[0], high[1], 1] = 0.
    rgb[zero[0], zero[1], :] = 255     # Finally, if anything is 0, make white

    img = PIL.Image.fromarray(rgb)
    img.save(name)
    print 'wrote %s' % name    


def triggers_png(name, arr, threshold1=6, threshold2=10, transpose=False, ytop_to_bottom=False):
    assert 0 < threshold1 < threshold2

    if not transpose:
        arr = np.transpose(arr)
        
    if not ytop_to_bottom:
        arr = arr[::-1]

    # 2D boolean arrays
    below_threshold1 = (arr < threshold1)
    below_threshold2 = (arr < threshold2)

    # below threshold1: scale range [0,threshold1] -> [0,1]
    t0 = arr / threshold1
    t0 = np.maximum(t0, 0.0001)    # 0.0001 instead of 0.0, to make roundoff-robust
    t0 = np.minimum(t0, 0.9999)    # 0.9999 instead of 1.0, to make roundoff-robust

    # below threshold1: use (dark blue) -> (dark red) scale, same as write_png()
    # between threshold1 and threshold2: yellow (r=g=255, b=51)
    # above threshold2: green (r=b=100, g=255)

    rgb = np.zeros((arr.shape[0], arr.shape[1], 3), dtype=np.uint8)
    rgb[:,:,0] = np.where(below_threshold1, 256*t0,     np.where(below_threshold2, 255, 100))
    rgb[:,:,1] = np.where(below_threshold1, 0,          np.where(below_threshold2, 255, 255))
    rgb[:,:,2] = np.where(below_threshold1, 256*(1-t0), np.where(below_threshold2, 51, 100))

    img = PIL.Image.fromarray(rgb)
    img.save(name)
    print 'wrote %s' % name    


class Variance_Estimates():
    def __init__(self, h5):
        self.var = self._read_h5(h5, 'variance')
        self.t = self._read_h5(h5, 'time')[0]  # [0] required because t was stored as 2-D array
        assert self.var.shape[1] == self.t.shape[0]
        size = (self.t[1] - self.t[0]) / 2.

        # Interpolate zeros
        x = len(self.var[0])
        indices = np.arange(x)
        for frequency in xrange(len(self.var)):
            nonzero = np.nonzero(self.var[frequency])[0]
            if len(nonzero) < 0.25 * x:
                self.var[frequency] = np.zeros((x))
            else:
                self.var[frequency] = np.interp(indices, nonzero, self.var[frequency, nonzero])

        # We print a one-time warning if the variance is requested outside the range (warn_tmin, warn_tmax).
        # (See eval() below.)
        self.warn_tmin = self.t[0] - size
        self.warn_tmax = self.t[-1] + size
        self.warn_flag = False

        
    def eval(self, t):
        if (not self.warn_flag) and ((t < self.warn_tmin) or (t > self.warn_tmax)):
            self.warn_flag = True  # only warn once
            print "Variance_Estimate: This variance file ranges approximately from times", self.warn_tmin, 'to', self.warn_tmax
            print ("Variance_Estimate: Requesting variances outside of this time range will result in eval() returning"
                   + " the endpoints of the variance array.")

        ret = []
        for f in range(self.var.shape[0]):
            ret += [ np.interp(t, self.t, self.var[f]) ] 
        return ret


    def _read_h5(self, fname, dset):
        with h5py.File(fname, 'r') as hf:
            return hf[dset][:]


####################################################################################################


class web_viewer_context_manager:
    """Helper class for run_for_web_viewer(): manages I/O redirection to log files, and renaming tmp_dir -> final_dir."""

    def __init__(self, tmp_dir, final_dir, redirect_to_log_files):
        assert sys.stdout.fileno() == 1
        assert sys.stderr.fileno() == 2

        self.tmp_dir = tmp_dir
        self.final_dir = final_dir
        self.redirect_to_log_files = redirect_to_log_files

        if os.path.exists(self.tmp_dir):
            raise RuntimeError("web viewer temporary directory '%s' already exists" % self.tmp_dir)
        if os.path.exists(self.final_dir):
            raise RuntimeError("web viewer output directory '%s' already exists" % self.tmp_dir)

        os.makedirs(self.tmp_dir)

        if self.redirect_to_log_files:
            self.stdout_log_filename = os.path.join(self.tmp_dir, 'rf_pipeline_stdout.txt')
            self.stderr_log_filename = os.path.join(self.tmp_dir, 'rf_pipeline_stderr.txt')

            self.tee_subprocess = subprocess.Popen(['tee', self.stderr_log_filename], stdin=subprocess.PIPE, stdout=sys.stderr.fileno())
            self.stdout_log_fileobj = open(self.stdout_log_filename, 'w')   # write to stdout logfile directly
            self.stderr_log_fileobj = self.tee_subprocess.stdin             # write to stderr logfile through tee


    def __enter__(self):
        if self.redirect_to_log_files:
            sys.stdout.flush()
            sys.stderr.flush()

            self.stdout_save = os.dup(sys.stdout.fileno())
            self.stderr_save = os.dup(sys.stderr.fileno())
            
            os.dup2(self.stdout_log_fileobj.fileno(), sys.stdout.fileno())
            os.dup2(self.stderr_log_fileobj.fileno(), sys.stderr.fileno())


    def __exit__(self, etype, value, tb):
        if self.redirect_to_log_files:
            sys.stdout.flush()
            sys.stderr.flush()
            os.dup2(self.stdout_save, sys.stdout.fileno())
            os.dup2(self.stderr_save, sys.stderr.fileno())
            os.close(self.stdout_save)
            os.close(self.stderr_save)

            self.stdout_log_fileobj.close()
            self.stderr_log_fileobj.close()
            self.tee_subprocess.wait()

            if etype is not None:
                self._try_to_append_traceback(self.stderr_log_filename, etype, value, tb)

        if self._maybe_rename():
            what = 'json file and log files' if self.redirect_to_log_files else 'json file'
            print 'rf_pipelines: pipeline %s successfully written to %s' % (what, self.final_dir)
        elif self.redirect_to_log_files:
            print 'rf_pipelines: pipeline json file could not be written, suggest inspecting log files rf_pipeline_std*.txt in the temporary directory', self.tmp_dir
        else:
            print 'rf_pipelines: pipeline json file could not be written, there may be some useful info in the temporary directory', self.tmp_dir


    def _try_to_append_traceback(self, filename, etype, value, tb):
        try:
            f = open(filename, 'a')
            traceback.print_exception(etype, value, tb, file=f)
        except:
            pass

            
    def _maybe_rename(self):
        """Returns true if pipeline output json file exists, is parseable, and directory rename succeeds."""

        try:
            filename = os.path.join(self.tmp_dir, 'rf_pipeline_0.json')
            json.load(open(filename, 'r'))
            os.rename(self.tmp_dir, self.final_dir)
        except:
            return False

        return True
    

def run_for_web_viewer(run_name, p, verbosity=2, redirect_to_log_files=False, extra_attrs=None):
    """
    Runs a pipeline, with output directory chosen appropriately for the web viewer
    at frb1.physics.mcgill.ca.

    The 'run_name' argument should be a short descriptive string.  The pipeline rundir 
    will look schematically like "(username)/(run_name)_(time)".

    The 'p' argument should be an rf_pipelines.pipeline_object.

    For the meanings of the 'verbosity' and 'extra_attrs' arguments, see the
    rf_pipelines.pipeline.run() docstring.

    If 'redirect_to_log_files' is True, then stdout/stderr of the running pipeline
    will be redirected to log files 'rf_pipeline_stdout.txt' and 'rf_pipeline_stderr.txt'
    in the web viewer directory.

    FIXME: what's the best way to generalize this function so that it makes sense
    on machines other than frb1.physics.mcgill.ca?
    """
    
    assert isinstance(run_name, basestring)
    assert isinstance(p, rf_pipelines_c.pipeline_object)

    # Directory names beginning with underscore are pipeline runs in progress.
    basename = '%s-%s' % (run_name, time.strftime('%y-%m-%d-%X'))
    dirname = os.path.join('/data2/web_viewer', os.environ['USER'])
    temp_dir = os.path.join(dirname, '_' + basename)
    final_dir = os.path.join(dirname, basename)

    print 'rf_pipelines: starting web_viewer run, temp_dir=%s, final_dir=%s' % (temp_dir, final_dir)

    if redirect_to_log_files:
        print 'rf_pipelines: redirecting output to log files in temp_dir, there will be little or no output while pipeline is running!'

    with web_viewer_context_manager(temp_dir, final_dir, redirect_to_log_files):
        p.run(outdir=temp_dir, clobber=False, verbosity=verbosity, extra_attrs=extra_attrs)
