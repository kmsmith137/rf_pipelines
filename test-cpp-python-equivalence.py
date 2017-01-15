#!/usr/bin/env python
#
# This tests equivalence of the C++ and python implementations of
# polynomial_detrender, intensity_clipper, and std_dev_clipper.

import sys
import numpy as np
import numpy.random as rand

import rf_pipelines
from rf_pipelines import rf_pipelines_c

# It's useful for debugging to have the same random data realizations every time.
rand.seed(1)


####################################################################################################
#
# General utils


def copy_array(arr, tame=False):
    """
    Make a float32 copy of an array, in a way which is artificially
    designed to make trouble for the python-to-C++ layer.
    """

    assert arr.ndim == 2

    pad = 0
    step = 1
    transpose = False

    if (not tame) and (rand.uniform() < 0.5):
        pad = rand.randint(0, 100)

    if (not tame) and (rand.uniform() < 0.5):
        step = rand.randint(2, 5)
    
    if (not tame) and (rand.uniform() < 0.5):
        arr = np.transpose(arr)
        transpose = True

    (nx, ny) = arr.shape

    ret = np.zeros((nx,ny*step + pad), dtype=np.float32)
    ret[:] = rand.uniform(-1.0e10, 1.0e10, size=ret.shape)

    ret = ret[:,:ny*step:step]
    ret[:,:] = arr[:,:]
    
    if transpose:
        ret = np.transpose(ret)

    return ret


def random_divisor(n):
    m = 2
    ret = 1

    while m**2 <= n:
        p = 0
        while (n % m) == 0:
            n = n / m
            p += 1

        ret = ret * m**rand.randint(0,p+1)
        m += 1

    if rand.uniform() < 0.5:
        ret = ret * n
        
    return ret


####################################################################################################


def test_utils():
    for iter in xrange(1000):
        Df = 2**rand.randint(0,6)
        Dt = 2**rand.randint(0,6)
        nfreq = Df * rand.randint(8,16)
        nt = Dt * 8 * rand.randint(1,8)

        sys.stderr.write('.')
        # print >>sys.stderr, '(Df,Dt,nfreq,nt)=(%d,%d,%d,%d)' % (Df,Dt,nfreq,nt)

        intensity = rand.uniform(size=(nfreq,nt))
        weights = rand.uniform(size=(nfreq,nt))

        #
        # Test 1: compare rf_pipelines.wi_downsample() and rf_pipelines_c.wi_downsample()
        #

        (ds_int, ds_wt) = rf_pipelines.wi_downsample(intensity, weights, nfreq//Df, nt//Dt)
        (ds_int2, ds_wt2) = rf_pipelines_c.wi_downsample(copy_array(intensity), copy_array(weights), Df, Dt)

        # Different weights convention used in python/C++ wi_downsample().
        ds_wt *= (Df*Dt)

        epsilon_w = np.max(np.abs(ds_wt - ds_wt2))
        epsilon_i = np.max(np.abs(ds_int - ds_int2))

        assert epsilon_w < 1.0e-3
        assert epsilon_i < 1.0e-3

        #
        # Test 2: compare rf_pipelines.weighted_mean_and_rms() and rf_pipelines_c.weighted_mean_and_rms(),
        # with (niter, Df, Dt, axis) = (1, 1, 1, None).
        #

        (mean1, rms1) = rf_pipelines.weighted_mean_and_rms(intensity, weights)
        (mean2, rms2) = rf_pipelines_c.weighted_mean_and_rms(copy_array(intensity), copy_array(weights), 3.0)

        epsilon_m = np.abs(mean1-mean2)
        epsilon_r = np.abs(rms1-rms2)

        assert epsilon_m < 1.0e-4
        assert epsilon_r < 1.0e-4

    print >>sys.stderr, 'test_utils: pass'


####################################################################################################
#
# Test polynomial detrender


def apply_reference_detrender(intensity, weights, axis, polydeg):
    """
    The python implementation of the polynomial_detrender is available as a transform, but
    there's no interface via a standalone function.  This little hack wraps a standalone-function
    interface around the transform.
    """

    (nfreq, nt) = intensity.shape
    assert weights.shape == (nfreq, nt)
    
    class fake_stream:
        def __init__(self, nfreq):
            self.nfreq = nfreq

    t = rf_pipelines.polynomial_detrender(deg=polydeg, axis=axis, nt_chunk=nt)
    t.set_stream(fake_stream(nfreq))
    t.process_chunk(0, 0, intensity, weights, None, None)


def random_sparse_vector(num_elts, num_nonzero):
    assert num_elts > num_nonzero

    ret = np.zeros(num_elts)

    for i in xrange(num_nonzero):
        ret[rand.randint(0,num_elts)] = rand.uniform()

    return ret


def make_detrender_test_data(nfreq, nt, axis, polydeg):
    """Helper function for detrender unit tests.  Returns 4-tuple (intensity, undoctored_weights, doctored_weights, zeroed_weights)."""

    # The python reference detrender includes hardcoded special behavior
    # if the sum of the weights is < 20, so we choose a large scale for
    # the weights, to avoid triggering this.

    intensity0 = rand.standard_normal(size=(nfreq,nt))
    weights0 = rand.uniform(100.0, 200.0, size=(nfreq,nt))

    doctored_weights = copy_array(weights0)
    zeroed_weights = copy_array(weights0)

    # In a few locations, we "doctor" the weights to make the polynomial
    # fit poorly behaved, and we "zero" the weights as well.

    if axis == 0:
        for it in xrange(nt):
            if rand.uniform() < 0.5:
                doctored_weights[:,it] = random_sparse_vector(nfreq, polydeg)
                zeroed_weights[:,it] = 0.

    elif axis == 1:
        for ifreq in xrange(nfreq):
            if rand.uniform() < 0.5:
                doctored_weights[ifreq,:] = random_sparse_vector(nt, polydeg)
                zeroed_weights[ifreq,:] = 0.

    return (intensity0, weights0, doctored_weights, zeroed_weights)


def test_polynomial_detrenders():
    for iter in xrange(1000):
        axis = rand.randint(0,2)
        polydeg = rand.randint(0,10)

        if axis == 0:
            nfreq = rand.randint(8*polydeg+16, 16*polydeg+32)
            nt = 8 * rand.randint(1, 10)
        elif axis == 1:
            nfreq = rand.randint(2, 20)
            nt = 8 * rand.randint(polydeg+2, 2*polydeg+4)

        sys.stderr.write('.')
        # print >>sys.stderr, 'axis=', axis, 'polydeg=', polydeg, 'nfreq=', nfreq, 'nt=', nt

        (intensity0, weights0, doctored_weights, zeroed_weights) = make_detrender_test_data(nfreq, nt, axis, polydeg)
        
        intensity1 = copy_array(intensity0)
        rf_pipelines_c.apply_polynomial_detrender(intensity1, doctored_weights, axis, polydeg)

        intensity2 = np.copy(intensity0)
        apply_reference_detrender(intensity2, np.copy(weights0), axis, polydeg)

        # We require that the arrays be equal where unmasked
        epsilon = np.max(np.abs(zeroed_weights*intensity1 - zeroed_weights*intensity2))

        # Debug
        # print 'weights0 =', weights0
        # print 'zeroed_weights =', zeroed_weights
        # print 'intensity0 =', intensity0
        # print 'intensity1 =', intensity1
        # print 'intensity2 =', intensity2
        # print '    ', epsilon

        assert np.array_equal(doctored_weights, zeroed_weights)
        assert epsilon < 1.0e-3

    print >>sys.stderr, 'test_polynomial_detrenders: pass'


####################################################################################################
#
# Test intensity_clipper and std_dev_clipper
#
# Note: test_clippers() only tests the intensity clipper in the case niter = 1.
# See later in the file for a unit test which addresses the niter > 1 case.


def make_clipper_test_data(nfreq, nt, axis, Df, Dt):
    """Helper function for clipper unit tests.  Returns (intensity, weights) pair."""

    intensity = rand.standard_normal(size=(nfreq,nt))
    weights = rand.uniform(size=(nfreq,nt))

    if rand.uniform() < 0.02:
        # Test a corner case by masking all elements of the array.
        weights[:,:] = 0.0
        
    elif rand.uniform() < 0.02:
        # Test a corner case by masking all elements of the array except one.
        ifreq = rand.randint(0, nfreq)
        it = rand.randint(0, nt)
        weights[:,:] = 0.0
        weights[ifreq,it] = rand.uniform()
        
    elif (axis == 0) and (rand.uniform() < 0.02):
        # Test a corner case by masking all columns of the array except one
        it = rand.randint(0, nt)
        weights[:,it] = rand.uniform(size=nfreq)
        
    elif (axis == 1) and (rand.uniform() < 0.02):
        # Test a corner case by masking all rows of the array except one
        ifreq = rand.randint(0, nfreq)
        weights[ifreq,:] = rand.uniform(size=nt)
        
    elif axis == 0:
        # Test a corner case by making a few columns "sparse"
        for j in xrange(nt//Dt):
            if rand.uniform() < 0.8:
                continue
            
            weights[:,(j*Dt):((j+1)*Dt)] = 0.
            if rand.uniform() < 0.5:
                continue
            
            i = rand.randint(0,nfreq//Df)
            weights[(i*Df):((i+1)*Df),(j*Dt):((j+1)*Dt)] = rand.uniform(size=(Df,Dt))

    elif axis == 1:
        # Test a corner case by making a few rows "sparse"
        for i in xrange(nfreq//Df):
            if rand.uniform() < 0.8:
                continue
            
            weights[(i*Df):((i+1)*Df),:] = 0.
            if rand.uniform() < 0.5:
                continue
            
            j = rand.randint(0,nt//Dt)
            weights[(i*Df):((i+1)*Df),(j*Dt):((j+1)*Dt)] = rand.uniform(size=(Df,Dt))

    return (intensity, weights)


def test_clippers():
    for iter in xrange(1000):
        sys.stderr.write('.')

        axis = rand.randint(0,2) if (rand.uniform() < 0.66) else None
        Df = 2**rand.randint(0,6)
        Dt = 2**rand.randint(0,6)
        nfreq = Df * rand.randint(8,16)
        nt = Dt * 8 * rand.randint(1,8)
        thresh = rand.uniform(1.1, 1.3)

        # Debug
        # print >>sys.stderr, '(Df,Dt,axis,nfreq,nt,thresh)=(%d,%d,%s,%d,%d,%s)' % (Df,Dt,axis,nfreq,nt,thresh)

        (intensity, weights0) = make_clipper_test_data(nfreq, nt, axis, Df, Dt)

        weights1 = copy_array(weights0, tame=True)
        rf_pipelines.clip_fx(intensity, weights1, thr = 0.999 * thresh, n_internal=1, axis=axis, dsample_nfreq=nfreq//Df, dsample_nt=nt//Dt, imitate_cpp=True)

        weights2 = copy_array(weights0, tame=True)
        rf_pipelines.clip_fx(intensity, weights2, thr = 1.001 * thresh, n_internal=1, axis=axis, dsample_nfreq=nfreq//Df, dsample_nt=nt//Dt, imitate_cpp=True)
            
        weights3 = copy_array(weights0)
        rf_pipelines_c.apply_intensity_clipper(copy_array(intensity), weights3, axis, thresh, Df=Df, Dt=Dt)

        ok = np.logical_and(weights1 <= weights3, weights3 <= weights2)

        if not np.all(ok):
            print >>sys.stderr, 'intensity_clipper failed for (Df,Dt,axis,nfreq,nt,thresh)=(%d,%d,%s,%d,%d,%s)' % (Df,Dt,axis,nfreq,nt,thresh)

            t = np.argmax(np.logical_not(ok))
            (ifreq, it) = np.unravel_index(t, ok.shape)
            print >>sys.stderr, 'failure at (ifreq,it)=', (ifreq,it)
            print >>sys.stderr, 'intensity:', intensity[ifreq,it]
            print >>sys.stderr, 'weights:', weights1[ifreq,it], weights3[ifreq,it], weights2[ifreq,it]
            sys.exit(1)
            
        if axis is None:
            continue   # std_dev clipper is not defined for axis=None

        weights1 = np.array(weights0, dtype=np.float32)
        rf_pipelines.filter_stdv(intensity, weights1, thr = 0.999 * thresh, axis = axis, dsample_nfreq = nfreq//Df, dsample_nt = nt//Dt, imitate_cpp = True)

        weights2 = np.array(weights0, dtype=np.float32)
        rf_pipelines.filter_stdv(intensity, weights2, thr = 1.001 * thresh, axis = axis, dsample_nfreq = nfreq//Df, dsample_nt = nt//Dt, imitate_cpp = True)
        
        weights3 = np.array(weights0, dtype=np.float32)
        rf_pipelines_c.apply_std_dev_clipper(intensity, weights3, axis, thresh, Df, Dt)

        ok = np.logical_and(weights1 <= weights3, weights3 <= weights2)
        
        if not np.all(ok):
            print >>sys.stderr, 'std_dev_clipper failed for (Df,Dt,axis,nfreq,nt,thresh)=(%d,%d,%s,%d,%d,%s)' % (Df,Dt,axis,nfreq,nt,thresh)

            t = np.argmax(np.logical_not(ok))
            (ifreq, it) = np.unravel_index(t, ok.shape)
            print >>sys.stderr, 'failure at (ifreq,it)=', (ifreq,it)
            print >>sys.stderr, 'intensity:', intensity[ifreq,it]
            print >>sys.stderr, 'weights:', weights1[ifreq,it], weights3[ifreq,it], weights2[ifreq,it]
            sys.exit(1)
                    
    print 'test_clippers: pass'


####################################################################################################
#
# The previously-defined test_clippers() tests the intensity_clipper for niter == 1.
#
# This test covers the niter > 1 case.  The logic here is a little tricky and the unit
# test is organized as a "correctness proof".



def test_iterated_intensity_clippers():
    for iter in xrange(1000):
        Df = 2**rand.randint(0,6)
        Dt = 2**rand.randint(0,6)
        nfreq = Df * rand.randint(8,16)
        nt = Dt * 8 * rand.randint(1,8)
        sigma = rand.uniform(1.1, 1.3)
        niter = rand.randint(1, 6)
        iter_sigma = rand.uniform(1.1, 1.3)
        
        sys.stderr.write('.')
        # print 'iteration %d: (Df,Dt,nfreq,nt,sigma,niter,iter_sigma) = %s' % (iter, (Df,Dt,nfreq,nt,sigma,niter,iter_sigma))

        (intensity0, weights0) = make_clipper_test_data(nfreq, nt, axis=None, Df=Df, Dt=Dt)

        intensity = copy_array(intensity0)

        # Test 1: AXIS_TIME iterated intensity_clipper is equivalent to 
        # looping over row blocks and running the AXIS_NONE clipper.

        weights1 = copy_array(weights0)
        weights2 = copy_array(weights0)

        # AXIS_TIME
        rf_pipelines_c.apply_intensity_clipper(intensity, weights1, 1, sigma, niter=niter, iter_sigma=iter_sigma, Df=Df, Dt=Dt)

        # AXIS_NONE
        for ifreq in xrange(nfreq//Df):
            iblock = intensity[(ifreq*Df):((ifreq+1)*Df),:]
            wblock = weights2[(ifreq*Df):((ifreq+1)*Df),:]
            rf_pipelines_c.apply_intensity_clipper(iblock, wblock, None, sigma, niter=niter, iter_sigma=iter_sigma, Df=Df, Dt=Dt)

        assert np.array_equal(weights1, weights2)

        # Test 2: AXIS_FREQ iterated intensity_clipper is equivalent to 
        # looping over column blocks and running the AXIS_NONE clipper.

        weights1 = copy_array(weights0)
        weights2 = copy_array(weights0)

        # AXIS_FREQ
        rf_pipelines_c.apply_intensity_clipper(intensity, weights1, 0, sigma, niter=niter, iter_sigma=iter_sigma, Df=Df, Dt=Dt)

        for it in xrange(nt//Dt):
            # need a little hacking to satisfy simd-derived rf_pipelines_c divisibility requirements...
            iblock = np.zeros((nfreq, 8*Dt))
            wblock = np.zeros((nfreq, 8*Dt))

            iblock[:,:Dt] = intensity[:,(it*Dt):((it+1)*Dt)]
            wblock[:,:Dt] = weights2[:,(it*Dt):((it+1)*Dt)]
        
            rf_pipelines_c.apply_intensity_clipper(iblock, wblock, None, sigma, niter=niter, iter_sigma=iter_sigma, Df=Df, Dt=Dt)

            weights2[:,(it*Dt):((it+1)*Dt)] = wblock[:,:Dt]

        assert np.array_equal(weights1, weights2)

        # At this point in the code, the task of proving correctness of the iterated
        # intensity_clipper has been reduced to the case axis=AXIS_NONE
        #
        # Test 3: intensity_clipper with downsampling factors (Df,Dt) is equivalent
        # to downsampling the array, and runnning intensity_clipper with (Dt,Dt)=(1,1).
        #
        # Note: this test depends on correctness of rf_pipelines_c.wi_downsample(),
        # which is independently unit-tested in test_utils() above.

        weights1 = copy_array(weights0)
        weights2 = copy_array(weights0)

        # AXIS_NONE
        rf_pipelines_c.apply_intensity_clipper(intensity, weights1, None, sigma, niter=niter, iter_sigma=iter_sigma, Df=Df, Dt=Dt)

        (ds_int, ds_wt) = rf_pipelines_c.wi_downsample(intensity, weights0, Df, Dt)
        rf_pipelines_c.apply_intensity_clipper(ds_int, ds_wt, None, sigma, niter=niter, iter_sigma=iter_sigma, Df=1, Dt=1)

        # Apply upsampled mask to weights2
        t = rf_pipelines.upsample(ds_wt, nfreq, nt)
        weights2 = np.where(t > 0, weights2, 0)

        assert np.array_equal(weights1, weights2)

        # At this point in the code, the task of proving correctness of the iterated
        # intensity_clipper has been reduced to the case (axis,Df,Dt) = (AXIS_NONE,1,1).
        #
        # Test 4: this test shows that correctness of intensity_clipper(niter) implies
        # correctness of weighted_mean_rms(niter+1).

        weights1 = copy_array(weights0)

        # (axis, Df, Dt, iter_sigma) = (None, 1, 1, sigma)
        rf_pipelines_c.apply_intensity_clipper(intensity, weights1, None, sigma, niter=niter, iter_sigma=sigma)
        
        (mean1, rms1) = rf_pipelines_c.weighted_mean_and_rms(intensity, weights1, sigma, 1)
        (mean2, rms2) = rf_pipelines_c.weighted_mean_and_rms(intensity, weights0, sigma, niter+1)
        (epsilon_m, epsilon_r) = (np.abs(mean1-mean2), np.abs(rms1-rms2))

        assert epsilon_m < 1.0e-6
        assert epsilon_r < 1.0e-6

        # Test 5: this test shows that correctness of weighted_mean_rms(niter) implies
        # correctness of intensity_clipper(niter).
        #
        # Taken together with test 4, this gives an inductive proof of correctness for
        # all niter > 1, which completes the test!

        weights1 = copy_array(weights0)
        weights2 = copy_array(weights0)

        rf_pipelines_c.apply_intensity_clipper(intensity, weights1, None, sigma, niter=niter, iter_sigma=iter_sigma)

        (mean, rms) = rf_pipelines_c.weighted_mean_and_rms(intensity, weights2, iter_sigma, niter)
        
        w32 = copy_array(weights0, tame=True)
        z32 = np.array(0.0, dtype=np.float32)

        weights_lo = np.where(np.abs(intensity-mean) < 0.999 * sigma * rms, w32, z32)
        weights_hi = np.where(np.abs(intensity-mean) < 1.001 * sigma * rms, w32, z32)
        
        assert np.all(weights_lo <= weights1)
        assert np.all(weights1 <= weights_hi)


    print >>sys.stderr, 'test_iterated_intensity_clippers: pass'


####################################################################################################
#
# Test transforms: tests equivalence between apply_* functions and transform objects.


def make_weird_data(nfreq, nt):
    """Returns (intensity, weights) pair, designed to be as weird as possible."""

    if rand.uniform() < 0.5:
        axis = rand.randint(0,2)
        polydeg = rand.randint(0, min(nt,10))
        (intensity, weights, w1, w2) = make_detrender_test_data(nfreq, nt, axis, polydeg)
        return (intensity, weights)

    else:
        axis = rand.randint(0,2) if (rand.uniform() < 0.66) else None
        Df = random_divisor(nfreq)
        Dt = random_divisor(nt)
        return make_clipper_test_data(nfreq, nt, axis, Df, Dt)
    

def make_random_transform():
    """Returns (transform, f_apply) pair."""

    transform_type = rand.randint(0,3)

    if transform_type == 0:
        # poly detrender
        axis = rand.randint(0,2)
        polydeg = rand.randint(0,10)
        epsilon = rand.uniform(0.01, 0.1)
        nt_chunk = 8 * rand.randint(polydeg+2, 2*polydeg+4)

        t = rf_pipelines_c.make_polynomial_detrender(nt_chunk, axis, polydeg, epsilon)
        f = lambda intensity, weights: rf_pipelines_c.apply_polynomial_detrender(intensity, weights, axis, polydeg, epsilon)

    elif transform_type == 1:
        # intensity_clipper
        axis = rand.randint(0,2) if (rand.uniform() < 0.66) else None
        Df = 2**rand.randint(0,6)
        Dt = 2**rand.randint(0,6)
        sigma = rand.uniform(1.1, 1.3)
        niter = rand.randint(1,5)
        iter_sigma = rand.uniform(1.8, 2.0)
        nt_chunk = Dt * 8 * rand.randint(1,8)

        t = rf_pipelines_c.make_intensity_clipper(nt_chunk, axis, sigma, niter, iter_sigma, Df, Dt)
        f = lambda intensity, weights: rf_pipelines_c.apply_intensity_clipper(intensity, weights, axis, sigma, niter=niter, iter_sigma=iter_sigma, Df=Df, Dt=Dt)

    else:
        # std_dev_clipper
        axis = rand.randint(0,2)
        Df = 2**rand.randint(0,6)
        Dt = 2**rand.randint(0,6)
        sigma = rand.uniform(1.1, 1.3)
        nt_chunk = Dt * 8 * rand.randint(1,8)

        t = rf_pipelines_c.make_std_dev_clipper(nt_chunk, axis, sigma, Df, Dt)
        f = lambda intensity, weights: rf_pipelines_c.apply_std_dev_clipper(intensity, weights, axis, sigma, Df=Df, Dt=Dt)

    assert t.nt_chunk == nt_chunk
    assert t.nt_prepad == 0
    assert t.nt_postpad == 0

    return (t, f)


class test_initializer(rf_pipelines.py_wi_transform):
    def __init__(self, f_apply, nt_chunk):
        rf_pipelines.py_wi_transform.__init__(self)

        self.f_apply = f_apply
        self.nt_chunk = nt_chunk
        self.nt_prepad = 0
        self.nt_postpad = 0

        self.expected_intensity = [ ]
        self.expected_weights = [ ]


    def set_stream(self, s):
        self.nfreq = s.nfreq


    def process_chunk(self, t0, t1, intensity, weights, pp_intensity, pp_weights):
        (intensity0, weights0) = make_weird_data(self.nfreq, self.nt_chunk)

        intensity[:,:] = intensity0[:,:]
        weights[:,:] = weights0[:,:]

        self.expected_intensity.append(copy_array(intensity0))
        self.expected_weights.append(copy_array(weights0))

        self.f_apply(self.expected_intensity[-1], self.expected_weights[-1])


class test_finalizer(rf_pipelines.py_wi_transform):
    def __init__(self, initializer):
        rf_pipelines.py_wi_transform.__init__(self)

        self.initializer = initializer
        self.nt_chunk = initializer.nt_chunk
        self.nt_prepad = 0
        self.nt_postpad = 0
        self.ichunk = 0

    def set_stream(self, s):
        self.nfreq = s.nfreq

    def process_chunk(self, t0, t1, intensity, weights, pp_intensity, pp_weights):
        assert np.array_equal(intensity, self.initializer.expected_intensity[self.ichunk])
        assert np.array_equal(weights, self.initializer.expected_weights[self.ichunk])
        self.ichunk += 1


def test_transforms():
    for iouter in xrange(20):
        sys.stderr.write('.')

        transform_chain = [ ]
        for itransform in xrange(10):
            (t, f) = make_random_transform()
            ti = test_initializer(f, t.nt_chunk)
            tf = test_finalizer(ti)
            transform_chain += [ ti, t, tf ]

        nfreq = 32 * rand.randint(1, 8)
        nt_tot = rand.randint(1000, 5000)
        freq_lo_MHz = 400.   # arbitrary
        freq_hi_MHz = 800.
        dt_sample = 1.0e-3
        stream = rf_pipelines.gaussian_noise_stream(nfreq, nt_tot, freq_lo_MHz, freq_hi_MHz, dt_sample)

        stream.run(transform_chain, outdir=None, noisy=False)

    print >>sys.stderr, 'test_iterated_intensity_clippers: pass'


####################################################################################################


test_utils()
test_polynomial_detrenders()
test_clippers()
test_iterated_intensity_clippers()
test_transforms()
