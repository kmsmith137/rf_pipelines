#!/usr/bin/env python
#
# Tests wi_sub_pipeline, in special case Dt=1 for now.
# Also indirectly tests jsonize/from_json() for a few transforms.
#
# FIXME cleanup: combine with test-cpp-python-equivalence.py

import numpy as np
import numpy.random as rand
import rf_pipelines


def make_random_transform():
    transform_type = rand.randint(0,3)

    if transform_type == 0:
        axis = 'freq'  # FIXME generalize later
        nbins = rand.randint(1, 5)
        nt_chunk = 8 * rand.randint(5, 11)
        epsilon = rand.uniform(3.0e-4, 1.0e-3)

        return rf_pipelines.spline_detrender(nt_chunk, axis, nbins, epsilon)

    elif transform_type == 1:
        # intensity_clipper
        axis = rand.randint(0,2) if (rand.uniform() < 0.66) else None
        Df = 2**rand.randint(0,4)
        Dt = 2**rand.randint(0,4)
        sigma = rand.uniform(1.3, 1.7)
        niter = rand.randint(1,5)
        iter_sigma = rand.uniform(1.8, 2.0)
        nt_chunk = Dt * 8 * rand.randint(1,8)
        two_pass = True if rand.randint(0,2) else False

        return rf_pipelines.intensity_clipper(nt_chunk, axis, sigma, niter, iter_sigma, Df, Dt, two_pass)

    else:
        # std_dev_clipper
        axis = rand.randint(0,2)
        Df = 2**rand.randint(0,4)
        Dt = 2**rand.randint(0,4)
        sigma = rand.uniform(1.3, 1.7)
        nt_chunk = Dt * 8 * rand.randint(1,8)
        two_pass = True if rand.randint(0,2) else False

        return rf_pipelines.std_dev_clipper(nt_chunk, axis, sigma, Df, Dt, two_pass)


def make_random_pipeline():
    n = rand.randint(1, 5)
    return rf_pipelines.pipeline([ make_random_transform() for i in xrange(n) ])


def make_random_pipeline_json():
    p = make_random_pipeline()
    j = p.jsonize()

    # throw in this test of jsonize()/from_json()
    jj = rf_pipelines.pipeline_object.from_json(j).jsonize()
    assert j == jj

    return j


####################################################################################################


class initial_stream(rf_pipelines.wi_stream):
    def __init__(self, intensity_arr, weights_arr, nt_chunk=None):
        assert intensity_arr.ndim == 2
        assert intensity_arr.shape == weights_arr.shape

        if nt_chunk is None:
            nt_chunk = rand.randint(10,20)

        rf_pipelines.wi_stream.__init__(self, 'initial_stream')

        self.nfreq = intensity_arr.shape[0]
        self.nt_chunk = nt_chunk
        self.nt_tot = intensity_arr.shape[1]
        self.intensity_arr = intensity_arr
        self.weights_arr = weights_arr


    def _fill_chunk(self, intensity, weights, pos):
        intensity[:,:] = 0.
        weights[:,:] = 0.

        if pos >= self.nt_tot:
            return False

        n = min(self.nt_tot - pos, self.nt_chunk)
        intensity[:,:n] = self.intensity_arr[:,pos:(pos+n)]
        weights[:,:n] = self.weights_arr[:,pos:(pos+n)]
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


def run_pipeline(pipeline_json, intensity_arr, weights_arr):
    # Just for fun, randomize 'nt_chunk'.
    p0 = initial_stream(intensity_arr, weights_arr)
    p1 = rf_pipelines.pipeline_object.from_json(pipeline_json)
    p2 = final_transform()
    
    p = rf_pipelines.pipeline([p0,p1,p2])
    p.run(outdir=None, verbosity=0, debug=True)

    (intensity, weights) = p2.get_results()
    return (intensity, weights)


####################################################################################################


def maxdiff(a1, a2):
    assert a1.shape == a2.shape
    return np.max(np.abs(a1-a2))


def run_test():
    Df = 2**rand.randint(0,5)
    nfreq = Df * 8 * rand.randint(10, 20)
    nt_tot = 8 * rand.randint(150, 500)
    input_intensity = rand.standard_normal(size=(nfreq,nt_tot))
    input_weights = rand.uniform(0.5, 1.0, size=(nfreq,nt_tot))
    p0_json = make_random_pipeline_json()
    p1_json = make_random_pipeline_json()
    p2_json = make_random_pipeline_json()
    
    # First run
    (i0,w0) = run_pipeline(p0_json, input_intensity, input_weights)
    (i0,w0) = (i0[:,:nt_tot], w0[:,:nt_tot])
    (i1,w1) = rf_pipelines.wi_downsample(i0, w0, Df, 1)
    (i2,w2) = run_pipeline(p1_json, i1, w1)
    (i2,w2) = (i2[:,:nt_tot], w2[:,:nt_tot])
    rf_pipelines.weight_upsample(w0, w2)
    (i3,w3) = run_pipeline(p2_json, i0, w0)
    
    # Second run

    si = initial_stream(input_intensity, input_weights)
    p0 = rf_pipelines.pipeline_object.from_json(p0_json)
    p1 = rf_pipelines.pipeline_object.from_json(p1_json)
    ps = rf_pipelines.wi_sub_pipeline(p1, Df=Df, Dt=1)
    p2 = rf_pipelines.pipeline_object.from_json(p2_json)
    tf = final_transform()

    p = rf_pipelines.pipeline([ si, p0, ps, p2, tf ])
    p.run(outdir=None, verbosity=0, debug=True)
    (i4,w4) = tf.get_results()

    eps_i = maxdiff((i3*w3)[:,:nt_tot],(i4*w4)[:,:nt_tot])
    eps_w = maxdiff(w3[:,:nt_tot], w4[:,:nt_tot])

    assert eps_i < 1.0e-5
    assert eps_w < 1.0e-5
    assert np.all(w3[:,nt_tot:] == 0.0)
    assert np.all(w4[:,nt_tot:] == 0.0)

    
####################################################################################################


niter = 100

for iter in xrange(100):
    if iter % 10 == 0:
        print 'test-wi-sub-pipeline: iteration %d/%d' % (iter, niter)
    run_test()

print 'test-wi-sub-pipeline: pass'
