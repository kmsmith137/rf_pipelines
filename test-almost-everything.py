#!/usr/bin/env python

import sys
import numpy as np
import numpy.random as rand

import rf_pipelines


###########################################  helpers  ##############################################


def maxdiff(a1, a2):
    assert a1.shape == a2.shape
    return np.max(np.abs(a1-a2))


def fact2(n, kmin=1):
    """Factorizes n as (2^m k), where k is odd, and returns m."""

    assert n >= kmin

    m = 0
    while ((n % 2) == 0) and (n >= 2*kmin):
        (m, n) = (m+1, n//2)

    return m


def wi_copy(intensity, weights, nt_chunk=None):
    assert intensity.shape == weights.shape
    assert intensity.ndim == 2

    (nfreq, nt) = intensity.shape
    assert nfreq > 0
    assert nt > 0

    if nt_chunk is None:
        return (np.copy(intensity), np.copy(weights))

    assert nt_chunk > 0
    nt_new = ((nt + nt_chunk - 1) // nt_chunk) * nt_chunk

    i_copy = np.empty((nfreq, nt_new), dtype=np.float32)
    i_copy[:,:nt] = intensity[:,:]
    i_copy[:,nt:] = rand.uniform(size=(nfreq,nt_new-nt))

    w_copy = np.empty((nfreq, nt_new), dtype=np.float32)
    w_copy[:,:nt] = weights[:,:]
    w_copy[:,nt:] = 0.

    return (i_copy, w_copy)

    

####################################################################################################
#
# make_random_*
#
# These functions return json-seralized random pipeline_objects (not the pipeline_objects themselves)


def make_random_polynomial_detrender(nfreq_ds, nds):
    assert nfreq_ds > 0
    assert nds > 0

    if (rand.randint(0,2)) and (nfreq_ds >= 16):
        axis = 'AXIS_FREQ'
        maxdeg = (nfreq_ds // 8) - 2
        maxdeg = min(maxdeg, 4)
        polydeg = rand.randint(0, maxdeg+1)
        nt_chunk = 8 * nds * max(rand.randint(-5,5), 0)
    else:
        axis = 'AXIS_TIME'
        polydeg = rand.randint(0,5)
        nt_chunk = 8 * nds * rand.randint(polydeg+2, 2*polydeg+4)

    return {
        'class_name': 'polynomial_detrender',
        'nt_chunk': nt_chunk,
        'axis': axis,
        'polydeg': polydeg,
        'epsilon': rand.uniform(1.0e-2, 2.0e-2)
    }


def make_random_spline_detrender(nfreq_ds, nds):
    assert nfreq_ds >= 32
    max_nbins = nfreq_ds // 32

    return {
        'class_name': 'spline_detrender',
        'nt_chunk': 8 * nds * max(rand.randint(-5,10),0),
        'axis': 'AXIS_FREQ',
        'nbins': rand.randint(1, min(max_nbins+1,5)),
        'epsilon': rand.uniform(3.0e-4, 6.0e-4)
    }


def make_random_intensity_clipper(nfreq_ds, nds):
    assert nfreq_ds >= 32

    a = rand.randint(0,3)
    Dt = 2**max(rand.randint(-3,5),0)

    if (a == 0) and (nfreq_ds >= 32):
        axis = 'AXIS_FREQ'
        m = fact2(nfreq_ds, kmin=32)
        Df = 2**max(rand.randint(-3,m+1),0)
        nt_chunk = 8 * nds * Dt * rand.randint(1, 8)
    elif a == 1:
        axis = 'AXIS_TIME'
        m = fact2(nfreq_ds)
        Df = 2**max(rand.randint(-3,m+1),0)
        nt_chunk = 8 * nds * Dt * rand.randint(8, 16)
    else:
        axis = 'AXIS_NONE'
        m = fact2(nfreq_ds)
        Df = 2**max(rand.randint(-3,m+1),0)
        p = (8 * Df) // nfreq_ds
        nt_chunk = 8 * nds * Dt * rand.randint(p+1, 10)

    return {
        'class_name': 'intensity_clipper',
        'axis': axis,
        'Df': Df,
        'Dt': Dt,
        'nt_chunk': nt_chunk,
        'sigma': rand.uniform(1.7, 2.0),
        'iter_sigma': rand.uniform(1.7, 2.0),
        'niter': max(rand.randint(-3,4), 1),
        'two_pass': True if rand.randint(0,2) else False
    }


def make_random_std_dev_clipper(nfreq_ds, nds):
    assert nfreq_ds >= 32
    assert nfreq_ds % 8 == 0

    axis = 'AXIS_FREQ' if rand.randint(0,2) else 'AXIS_TIME'
    m = fact2(nfreq_ds, kmin=4)
    Df = 2**max(rand.randint(-3,m-2),0)
    Dt = 2**max(rand.randint(-3,5),0)
    nt_chunk = 8 * nds * Dt * rand.randint(8,16)

    return {
        'class_name': 'std_dev_clipper',
        'axis': axis,
        'Df': Df,
        'Dt': Dt,
        'nt_chunk': nt_chunk,
        'sigma': rand.uniform(1.7, 2.0),
        'two_pass': True if rand.randint(0,2) else False
    }



def make_random_transform_list(nfreq_ds, nds, nelements):
    assert nelements > 0

    ret = [ ]

    while len(ret) < nelements:
        r = rand.randint(0,3)

        if (r == 0):
            ret.append(make_random_polynomial_detrender(nfreq_ds, nds))
        elif (r == 1) and (nfreq_ds >= 32):
            ret.append(make_random_spline_detrender(nfreq_ds, nds))
        elif (r == 2):
            ret.append(make_random_intensity_clipper(nfreq_ds, nds))
        elif (r == 3) and (nfreq_ds % 8 == 0) and (nfreq_ds >= 32):
            ret.append(make_random_std_dev_clipper(nfreq_ds, nds))

    return ret


def make_random_pipeline(nfreq_ds, nds, nelements):
    t = make_random_transform_list(nfreq_ds, nds, nelements)
    return rf_pipelines.pipeline(t) if (len(t) > 1) else t[0]


####################################################################################################
#
# emulate_pipeline


def emulate_pipeline(pipeline_json, intensity, weights):
    """Returns new (intensity, weights) pair."""

    assert intensity.ndim == 2
    assert intensity.shape == weights.shape
    (nfreq, nt) = intensity.shape

    if pipeline_json['class_name'] == 'pipeline':
        for q in pipeline_json['elements']:
            (intensity, weights) = emulate_pipeline(q, intensity, weights)
        return (intensity, weights)

    if pipeline_json['class_name'] == 'polynomial_detrender':
        axis = pipeline_json['axis']
        polydeg = pipeline_json['polydeg']
        nt_chunk = pipeline_json['nt_chunk']
        epsilon = pipeline_json['epsilon']

        if axis == 'AXIS_FREQ':
            (i_copy, w_copy) = wi_copy(intensity, weights, 8)
            rf_pipelines.apply_polynomial_detrender(i_copy, w_copy, axis, polydeg, epsilon)
            return (i_copy, w_copy)

        if axis == 'AXIS_TIME':
            (i_copy, w_copy) = wi_copy(intensity, weights, nt_chunk)
            for it in xrange(0, nt, nt_chunk):
                rf_pipelines.apply_polynomial_detrender(i_copy[:,(it):(it+nt_chunk)], w_copy[:,(it):(it+nt_chunk)], axis, polydeg, epsilon)
            return (i_copy, w_copy)

        raise RuntimeError('emulate_pipeline: unsupported polynomial_detrender axis "%s"' % axis)


    if pipeline_json['class_name'] == 'spline_detrender':
        axis = pipeline_json['axis']
        nbins = pipeline_json['nbins']
        epsilon = pipeline_json['epsilon']

        if axis == 'AXIS_FREQ':
            (i_copy, w_copy) = wi_copy(intensity, weights, 8)
            rf_pipelines.apply_spline_detrender(i_copy, w_copy, axis, nbins, epsilon)
            return (i_copy, w_copy)

        raise RuntimeError('emulate_pipeline: unsupported polynomial_detrender axis "%s"' % axis)


    if pipeline_json['class_name'] == 'intensity_clipper':
        axis = pipeline_json['axis']
        Df = pipeline_json['Df']
        Dt = pipeline_json['Dt']
        niter = pipeline_json['niter']
        sigma = pipeline_json['sigma']
        iter_sigma = pipeline_json['iter_sigma']
        nt_chunk = pipeline_json['nt_chunk']
        two_pass = pipeline_json['two_pass']

        if axis == 'AXIS_FREQ':
            (i_copy, w_copy) = wi_copy(intensity, weights, 8*Dt)
            rf_pipelines.apply_intensity_clipper(i_copy, w_copy, axis, sigma, niter, iter_sigma, Df, Dt, two_pass)
            return (i_copy, w_copy)

        (i_copy, w_copy) = wi_copy(intensity, weights, nt_chunk)

        for it in xrange(0, nt, nt_chunk):
            rf_pipelines.apply_intensity_clipper(i_copy[:,(it):(it+nt_chunk)], w_copy[:,(it):(it+nt_chunk)], axis, sigma, niter, iter_sigma, Df, Dt, two_pass)

        return (i_copy, w_copy)

    raise RuntimeError('emulate_pipeline: unsupported class_name "%s"' % pipeline_json['class_name'])


####################################################################################################


class initial_stream(rf_pipelines.wi_stream):
    def __init__(self):
        rf_pipelines.wi_stream.__init__(self, 'initial_stream')
        self.nfreq = 2**rand.randint(0, 10)
        self.nfreq *= rand.randint(128//self.nfreq + 1, 2048//self.nfreq + 1)
        self.nt_tot = rand.randint(1000, 2000)
        self.intensity = rand.standard_normal(size=(self.nfreq,self.nt_tot))
        self.weights = rand.uniform(0.5, 1.0, size=(self.nfreq,self.nt_tot))
        self.nt_chunk = rand.randint(20, 50)


    def _fill_chunk(self, intensity, weights, pos):
        intensity[:,:] = 0.
        weights[:,:] = 0.

        if pos >= self.nt_tot:
            return False

        n = min(self.nt_tot - pos, self.nt_chunk)
        intensity[:,:n] = self.intensity[:,pos:(pos+n)]
        weights[:,:n] = self.weights[:,pos:(pos+n)]
        return True            


class final_transform(rf_pipelines.wi_transform):
    def __init__(self, nt_chunk=None):
        if nt_chunk is None:
            nt_chunk = rand.randint(10,20)

        rf_pipelines.wi_transform.__init__(self, "final_transform")
        self.nt_chunk = nt_chunk
        self.intensity_chunks = [ ]
        self.weight_chunks = [ ]

    def _process_chunk(self, intensity, weights, pos):
        self.intensity_chunks.append(np.copy(intensity))
        self.weight_chunks.append(np.copy(weights))

    def get_results(self):
        intensity = np.concatenate(self.intensity_chunks, axis=1)
        weights = np.concatenate(self.weight_chunks, axis=1)
        return (intensity, weights)


def run_test():
    s = initial_stream()
    u = final_transform()
    tj = make_random_transform_list(s.nfreq, 1, nelements=rand.randint(10,20))
    t = [ rf_pipelines.pipeline_object.from_json(j) for j in tj ]

    p = rf_pipelines.pipeline([s] + t + [u])
    p.bind()

    # Check jsonization (test is slightly stronger if bind() comes first)
    tj2 = [ x.jsonize() for x in t ]
    assert tj == tj2

    # First run
    p.run(outdir=None, verbosity=0)
    (i0,w0) = u.get_results()

    nt0 = i0.shape[1]
    assert nt0 >= s.nt_tot
    assert i0.shape == w0.shape == (s.nfreq, nt0)
    assert np.all(w0[:,s.nt_tot:] == 0.0)
    
    # Second run
    pj = { 'class_name': 'pipeline', 'elements': tj }
    (i1,w1) = emulate_pipeline(pj, s.intensity, s.weights)
    
    nt1 = i1.shape[1]
    assert nt1 >= s.nt_tot
    assert i1.shape == w1.shape == (s.nfreq, nt1)
    assert np.all(w1[:,s.nt_tot:] == 0.0)

    # Compare
    eps_i = maxdiff((i0*w0)[:,:s.nt_tot], (i1*w1)[:,:s.nt_tot])
    eps_w = maxdiff(w0[:,:s.nt_tot], w1[:,:s.nt_tot])

    assert eps_i < 1.0e-5
    assert eps_w < 1.0e-5


####################################################################################################


niter = 100

for iter in xrange(niter):
    print 'iteration %d/%d' % (iter, niter)
    run_test()
