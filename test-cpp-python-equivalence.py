#!/usr/bin/env python

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

    t = rf_pipelines.legendre_detrender(polydeg, axis, nt)
    t.set_stream(fake_stream(nfreq))
    t.process_chunk(0, 0, intensity, weights, None, None)


def random_sparse_vector(num_elts, num_nonzero):
    assert num_elts > num_nonzero

    ret = np.zeros(num_elts)

    for i in xrange(num_nonzero):
        ret[rand.randint(0,num_elts)] = rand.uniform()

    return ret


def test_detrender():
    for iter in xrange(1000):
        axis = rand.randint(0,2)
        polydeg = rand.randint(0,10)

        if axis == 0:
            nfreq = rand.randint(8*polydeg+16, 16*polydeg+32)
            nt = 8 * rand.randint(1, 10)
        elif axis == 1:
            nfreq = rand.randint(2, 20)
            nt = 8 * rand.randint(polydeg+2, 2*polydeg+4)

        # Debug
        # (axis, polydeg, nfreq, nt) = (xx, xx, xx, xx)
        # print 'axis=', axis, 'polydeg=', polydeg, 'nfreq=', nfreq, 'nt=', nt

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

    print 'test_detrender: pass'

test_detrender()
sys.exit(0)


####################################################################################################


def remove_axis(arr, ax):
    if ax == 0:
        return arr[0,:]
    elif ax == 1:
        return arr[:,0]
    else:
        raise RuntimeError('bad call to remove_axis')


nfreq = 256
nt = 512
thresh = 1.2
axis = 0



for iter in xrange(10):
    print >>sys.stderr, 'starting clipper testing round', (iter+1)

    for Df in [ 2**i for i in xrange(6) ]:
        for Dt in [ 2**i for i in xrange(6) ]:
            for axis in [ None, 0, 1 ]:
                intensity = rand.standard_normal(size=(nfreq,nt))
                weights0 = rand.uniform(size=(nfreq,nt))

                weights1 = np.array(weights0, dtype=np.float32)
                rf_pipelines.clip_fx(intensity, weights1, thr = 0.999 * thresh, n_internal=1, axis=axis, dsample_nfreq=nfreq//Df, dsample_nt=nt//Dt, imitate_cpp=True)

                weights2 = np.array(weights0, dtype=np.float32)
                rf_pipelines.clip_fx(intensity, weights2, thr = 1.001 * thresh, n_internal=1, axis=axis, dsample_nfreq=nfreq//Df, dsample_nt=nt//Dt, imitate_cpp=True)
            
                weights3 = np.array(weights0, dtype=np.float32)
                rf_pipelines_c.apply_intensity_clipper(intensity, weights3, axis, thresh, Df=Df, Dt=Dt)

                ok = np.logical_and(weights1 <= weights3, weights3 <= weights2)

                if not np.all(ok):
                    print >>sys.stderr, 'intensity_clipper failed for (Df,Dt,axis)=(%d,%d,%s)' % (Df,Dt,axis)

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

                if (not np.all(weights1 <= weights3)) or (not np.all(weights3 <= weights2)):
                    weights1 = remove_axis(weights1, axis)
                    weights2 = remove_axis(weights2, axis)
                    weights3 = remove_axis(weights3, axis)
                    
                    for (w1, w3, w2) in zip(weights1, weights3, weights2):
                        print w1, w3, w2
                    
                    raise RuntimeError('sd_clipper failed for (Df,Dt,axis)=(%d,%d,%s)' % (Df,Dt,axis))
                

####################################################################################################



