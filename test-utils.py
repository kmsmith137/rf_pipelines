#!/usr/bin/env python

import sys
import numpy as np
import numpy.random as rand

import rf_pipelines


def take(arr, i, axis):
    """Replacement for numpy.take(), which is broken in some versions of numpy."""

    assert 0 <= axis < len(arr.shape)

    a = np.reshape(arr, (np.prod(arr.shape[:axis], dtype=np.int),
                         arr.shape[axis], 
                         np.prod(arr.shape[(axis+1):], dtype=np.int)))

    a = a[:,i,:]
    
    return np.reshape(a, arr.shape[:axis] + arr.shape[(axis+1):])
    


def test_expand_array():
    for iter in xrange(1000):
        out_ndim = rand.randint(1, 5)
        out_shape = tuple(rand.randint(1, 8, size=out_ndim))

        if rand.uniform() < 0.3:
            # axis=None
            t = rand.uniform()
            a = rf_pipelines.utils.expand_array(t, out_shape, axis=None)
            assert np.all(a == t)

        else:
            axis = rand.randint(0, out_ndim)
            in_shape = tuple(out_shape[:axis]) + tuple(out_shape[(axis+1):])

            t = rand.uniform(size=in_shape)
            u = rf_pipelines.utils.expand_array(t, out_shape, axis=axis)
            
            for i in xrange(out_shape[axis]):
                assert np.array_equal(t, take(u,i,axis=axis))

    print >>sys.stderr, 'test_expand_array: pass'


def test_weighted_mean_rms():
    """
    Here we just test consistency between the (axis == None) and (axis != None) cases.

    There is also a unit test (in test-cpp-python-equivalence.cpp) which tests the
    (axis == None) and (niter == 1) case, by comparing to the C++ equivalent.

    FIXME strictly speaking, in order to delare weighted_mean_and_rms() completely
    unit-tested, we would need a unit test for the case (axis == None) and (niter > 1).
    This is omitted for now, since it's nontrivial, and a quick-and-dirty empirical
    test looked fine.
    """
    
    for iter in xrange(1000):
        sys.stderr.write('.')

        nfreq = rand.randint(100, 200)
        nt = rand.randint(100, 200)
        niter = rand.randint(1, 5)
        axis = rand.randint(0, 2)
        sigma = rand.uniform(1.6, 1.8)

        intensity = rand.standard_normal(size=(nfreq,nt))
        weights = rand.uniform(size=(nfreq,nt))
        
        (mean, rms) = rf_pipelines.weighted_mean_and_rms(intensity, weights, niter, sigma, axis)

        mean2 = np.zeros_like(mean)
        rms2 = np.zeros_like(rms)

        if axis == 0:
            assert mean.shape == rms.shape == (nt,)
            for it in xrange(nt):
                (mean2[it], rms2[it]) = rf_pipelines.weighted_mean_and_rms(intensity[:,it], weights[:,it], niter, sigma, axis=None)

        else:
            assert mean.shape == rms.shape == (nfreq,)
            for ifreq in xrange(nfreq):
                (mean2[ifreq], rms2[ifreq]) = rf_pipelines.weighted_mean_and_rms(intensity[ifreq,:], weights[ifreq,:], niter, sigma, axis=None)
        
        assert np.array_equal(mean, mean2)
        assert np.array_equal(rms, rms2)

    print >>sys.stderr, 'test_weighted_mean_rms: pass'


test_expand_array()
test_weighted_mean_rms()
