#!/usr/bin/env python
#
# This unit test builds up a random pipeline (make_random_*()), and compares with a reference 
# python implementation (emulate_pipeline()).
#
# In some cases, for example intensity_clipper, the reference python implementation calls the
# same C++ kernel which is called by the real pipeline.  In these cases, we are not testing 
# correctness of the C++ kernel, we are just testing consistency of the pipeline ring-buffering/
# chunking logic.  (There are other unit tests which test correctness of the kernels!)
#
# Note that the "real" pipeline is run first, then the reference pipeline.  (This detail matters
# in a few cases, such as testing mask_(de)serializer, see below.)


import os
import sys
import glob
import json
import tempfile
import numpy as np
import numpy.random as rand

import rf_pipelines

have_hdf5 = False

try:
    import h5py
    have_hdf5 = True
except ImportError:
    print "test-almost-everything.py: 'import h5py' failed; pipeline features which use hdf5 will not be tested"


###########################################  helpers  ##############################################


def xdiv(m, n):
    assert n > 0
    assert m >= 0
    return m // n


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


def temporary_hdf5_filename():
    (f, filename) = tempfile.mkstemp(dir='.', suffix='.h5')
    return filename
    

####################################################################################################


class pipeline_params:
    """
    Represents high-level pipeline parameters (number of frequency channels, time samples. etc.)

    FIXME: we define a timestamp range using attributes named 'fpga_counts_per_sample' 
    and  'initial_fpga_count'.  This is currently expected by some transforms, for example
    mask_(de)serializer.  Need to settle on a systematic general convention for timestamps!
    """

    def __init__(self, nfreq_ds=None, nds=1, nt_tot=None, initial_fpga_count=None, fpga_counts_per_sample=384):
        """If any constructor arguments are None, they will be replaced by randomly generated values."""

        if nfreq_ds is None:
            nfreq_ds = 2**rand.randint(0, 10)
            nfreq_ds *= rand.randint(128//nfreq_ds + 1, 2048//nfreq_ds + 1)

        # We choose endpoint timestamps which are multiples of 256,
        # since this is currently assumed by mask_(de)serializer.

        if nt_tot is None:
            nt_tot = 256 * rand.randint(40, 100)

        if initial_fpga_count is None:
            initial_fpga_count = 256 * rand.randint(1024*1024, 2*1024*1024) * fpga_counts_per_sample

        self.nfreq_ds = nfreq_ds
        self.nds = nds
        self.nt_tot = nt_tot
        self.initial_fpga_count = initial_fpga_count
        self.fpga_counts_per_sample = fpga_counts_per_sample

    
    def downsample(self, Df, Dt):
        assert (Df > 0) and (Dt > 0) and (self.nfreq_ds % Df == 0)
        return pipeline_params(self.nfreq_ds // Df, self.nds * Dt, self.nt_tot, self.initial_fpga_count, self.fpga_counts_per_sample)


####################################################################################################
#
# make_random_*
#
# These functions take a 'pipeline_params' argument, and return json-seralized random pipeline_objects 
# (not the pipeline_objects themselves)


def make_random_polynomial_detrender(p):
    assert p.nfreq_ds > 0
    assert p.nds > 0

    if (rand.randint(0,2)) and (p.nfreq_ds >= 16):
        axis = 'AXIS_FREQ'
        maxdeg = (p.nfreq_ds // 8) - 2
        maxdeg = min(maxdeg, 4)
        polydeg = rand.randint(0, maxdeg+1)
        nt_chunk = 8 * p.nds * max(rand.randint(-5,5), 0)
    else:
        axis = 'AXIS_TIME'
        polydeg = rand.randint(0,5)
        nt_chunk = 8 * p.nds * rand.randint(polydeg+2, 2*polydeg+4)

    return {
        'class_name': 'polynomial_detrender',
        'nt_chunk': nt_chunk,
        'axis': axis,
        'polydeg': polydeg,
        'epsilon': rand.uniform(1.0e-2, 2.0e-2)
    }


def make_random_spline_detrender(p):
    assert p.nfreq_ds >= 32
    max_nbins = p.nfreq_ds // 32

    return {
        'class_name': 'spline_detrender',
        'nt_chunk': 8 * p.nds * max(rand.randint(-5,10),0),
        'axis': 'AXIS_FREQ',
        'nbins': rand.randint(1, min(max_nbins+1,5)),
        'epsilon': rand.uniform(3.0e-4, 6.0e-4)
    }


def make_random_intensity_clipper(p):
    assert p.nfreq_ds >= 32

    a = rand.randint(0,3)
    Dt = 2**max(rand.randint(-3,5),0)

    if (a == 0) and (p.nfreq_ds >= 32):
        axis = 'AXIS_FREQ'
        m = fact2(p.nfreq_ds, kmin=32)
        Df = 2**max(rand.randint(-3,m+1),0)
        nt_chunk = 8 * p.nds * Dt * rand.randint(1, 8)
    elif a == 1:
        axis = 'AXIS_TIME'
        m = fact2(p.nfreq_ds)
        Df = 2**max(rand.randint(-3,m+1),0)
        nt_chunk = 8 * p.nds * Dt * rand.randint(8, 16)
    else:
        axis = 'AXIS_NONE'
        m = fact2(p.nfreq_ds)
        Df = 2**max(rand.randint(-3,m+1),0)
        n = (8 * Df) // p.nfreq_ds
        nt_chunk = 8 * p.nds * Dt * rand.randint(n+1, 10)

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


def make_random_std_dev_clipper(p):
    assert p.nfreq_ds >= 32
    assert p.nfreq_ds % 8 == 0

    axis = 'AXIS_FREQ' if rand.randint(0,2) else 'AXIS_TIME'
    m = fact2(p.nfreq_ds, kmin=4)
    Df = 2**max(rand.randint(-3,m-2),0)
    Dt = 2**max(rand.randint(-3,5),0)
    nt_chunk = 8 * p.nds * Dt * rand.randint(8,16)

    return {
        'class_name': 'std_dev_clipper',
        'axis': axis,
        'Df': Df,
        'Dt': Dt,
        'nt_chunk': nt_chunk,
        'sigma': rand.uniform(1.7, 2.0),
        'two_pass': True if rand.randint(0,2) else False
    }


def make_random_mask_serializer(p):

    # We test the mask_serializer as follows:
    #
    #   - When the "real" pipeline runs, it serializes the mask
    #     to a temporary file.
    #
    #   - When the reference pipeline runs, it reads the mask file
    #     and compares it to the mask in the reference pipeline.
    #     It then cleans up by deleting the temporary file.

    assert (p.nds == 1) and have_hdf5

    return {
        'class_name': 'mask_serializer',
        'hdf5_filename': temporary_hdf5_filename()
    }


def make_random_mask_deserializer(p):

    # We test the mask_deserializer as follows:
    #
    #   - When the random pipeline is generated, we create a random 
    #     bitmask and write it to a temporary file.  When both the
    #     "real" and reference pipelines run, they will read this
    #     file and apply the bitmask.
    #
    #   - When the reference pipeline runs (after the "real" pipeline), 
    #     it also cleans up by deleting the temporary file.

    assert (p.nds == 1) and have_hdf5

    # These restrictions are convenient, but could be relaxed if necessary.
    assert (p.nt_tot % 256) == 0
    assert (p.initial_fpga_count % (256 * p.fpga_counts_per_sample) == 0)

    nf = p.fpga_counts_per_sample
    it0 = xdiv(p.initial_fpga_count, nf)
    it1 = it0 + p.nt_tot

    # Expand mask
    it0 -= 256 * max(rand.randint(-10,10), 0)
    it1 += 256 * max(rand.randint(-10,10), 0)

    filename = temporary_hdf5_filename()
    shape = (p.nfreq_ds, xdiv(it1-it0,8))

    # For file format, see mask_serializer.cpp.
    f = h5py.File(filename, 'w')
    f.attrs['nfreq'] = p.nfreq_ds
    f.attrs['initial_fpga_count'] = it0 * nf
    f.attrs['fpga_counts_per_sample'] = nf

    ds = f.create_dataset('bitmask', shape=shape, dtype=np.uint8)
    ds[:,:] = rand.randint(0, 256, size=shape)

    return {
        'class_name': 'mask_serializer',
        'hdf5_filename': filename
    }


def make_random_transform_list(p, nelements):
    assert nelements > 0

    ret = [ ]

    while (rand.uniform() < 0.5) and (nelements > 1):
        n = rand.randint(1, nelements+1)
        ret.append(make_random_pipeline(p, n))
        nelements -= (n-1)

    while len(ret) < nelements:
        # If more transforms are added to the list below, 
        # don't forget to increment the randint() upper limit!

        r = rand.randint(0,6)

        if (r == 0):
            ret.append(make_random_polynomial_detrender(p))
        elif (r == 1) and (p.nfreq_ds >= 32):
            ret.append(make_random_spline_detrender(p))
        elif (r == 2) and (p.nfreq_ds >= 32):
            ret.append(make_random_intensity_clipper(p))
        elif (r == 3) and (p.nfreq_ds % 8 == 0) and (p.nfreq_ds >= 32):
            ret.append(make_random_std_dev_clipper(p))
        elif (r == 4) and (rand.uniform() < 0.1) and (p.nds == 1) and have_hdf5:
            ret.append(make_random_mask_serializer(p))
        elif (r == 5) and (rand.uniform() < 0.1) and (p.nds == 1) and have_hdf5:
            ret.append(make_random_mask_deserializer(p))

    return ret



def make_random_pipeline(p, nelements, allow_downsampling=True):
    m = fact2(p.nfreq_ds)
    n = 5 - fact2(p.nds)
    Df = 2**max(rand.randint(-5,m+1),0) if allow_downsampling else 1
    Dt = 2**max(rand.randint(-5,n+1),0) if allow_downsampling else 1

    p = p.downsample(Df, Dt)
    tl = make_random_transform_list(p, nelements)

    if (len(tl) > 1) or (rand.uniform() < 0.1):
        ret = { 'class_name': 'pipeline', 'name': 'pipeline', 'elements': tl }
    else:
        ret = tl[0]
    
    if (Df == 1) and (Dt == 1) and (rand.random() < 0.9):
        return ret

    ret = { 
        'class_name': 'wi_sub_pipeline',
        'sub_pipeline': ret,
        'w_cutoff': rand.uniform(0.0, 0.01),
        'nfreq_out': p.nfreq_ds,
        'nds_out': p.nds,
        'Df': Df,
        'Dt': Dt
    }

    if rand.uniform() < 0.33:
        ret['nfreq_out'] = 0
    elif rand.uniform() < 0.5:
        ret['Df'] = 0

    if rand.uniform() < 0.33:
        ret['nds_out'] = 0
    elif rand.uniform() < 0.5:
        ret['Dt'] = 0

    return ret


####################################################################################################
#
# emulate_pipeline


def emulate_pipeline(pipeline_json, intensity, weights, params):
    """
    Returns new (intensity, weights) pair.

    Note that emulate_pipeline() runs second, after the real pipeline run.  The consequence
    of this is that transforms which write files (such as 'mask_serializer') can be unit-tested
    by reading back the file and checking its contents.
    """

    assert intensity.ndim == 2
    assert intensity.shape == weights.shape

    (nfreq, nt_ds) = intensity.shape
    nds = params.nds

    if pipeline_json['class_name'] == 'pipeline':
        for q in pipeline_json['elements']:
            (intensity, weights) = emulate_pipeline(q, intensity, weights, params)
        return (intensity, weights)

    if pipeline_json['class_name'] == 'wi_sub_pipeline':
        sub_pipeline = pipeline_json['sub_pipeline']
        w_cutoff = pipeline_json['w_cutoff']
        nfreq_out = pipeline_json['nfreq_out']
        nds_out = pipeline_json['nds_out']
        Df = pipeline_json['Df']
        Dt = pipeline_json['Dt']

        if (nfreq_out > 0) and (Df > 0):
            assert nfreq == nfreq_out * Df
        elif (nfreq_out == 0) and (Df > 0):
            assert nfreq % Df == 0
            nfreq_out = nfreq // Df
        elif (nfreq_out > 0) and (Df == 0):
            assert nfreq % nfreq_out == 0
            Df = nfreq // nfreq_out

        if (nds_out > 0) and (Dt > 0):
            assert nds_out == nds * Dt
        elif (nds_out == 0) and (Dt > 0):
            nds_out = nds * Dt
        elif (nds_out > 0) and (Dt == 0):
            assert nds_out % nds == 0
            Dt = nds_out // nds

        (i_copy, w_copy) = wi_copy(intensity, weights, 8*Dt)
        (i_ds, w_ds) = rf_pipelines.wi_downsample(i_copy, w_copy, Df, Dt)
        n = i_ds.shape[1]
        
        (i_ds, w_ds) = emulate_pipeline(sub_pipeline, i_ds, w_ds, params.downsample(Df,Dt))
        rf_pipelines.weight_upsample(w_copy, w_ds[:,:n], w_cutoff)

        return (i_copy, w_copy)


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
            n = xdiv(nt_chunk, nds)
            (i_copy, w_copy) = wi_copy(intensity, weights, n)
            for it in xrange(0, nt_ds, n):
                rf_pipelines.apply_polynomial_detrender(i_copy[:,(it):(it+n)], w_copy[:,(it):(it+n)], axis, polydeg, epsilon)
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

        n = xdiv(nt_chunk, nds)
        (i_copy, w_copy) = wi_copy(intensity, weights, n)

        for it in xrange(0, nt_ds, n):
            rf_pipelines.apply_intensity_clipper(i_copy[:,(it):(it+n)], w_copy[:,(it):(it+n)], axis, sigma, niter, iter_sigma, Df, Dt, two_pass)

        return (i_copy, w_copy)


    if pipeline_json['class_name'] == 'std_dev_clipper':
        axis = pipeline_json['axis']
        Df = pipeline_json['Df']
        Dt = pipeline_json['Dt']
        sigma = pipeline_json['sigma']
        nt_chunk = pipeline_json['nt_chunk']
        two_pass = pipeline_json['two_pass']

        n = xdiv(nt_chunk, nds)
        (i_copy, w_copy) = wi_copy(intensity, weights, n)

        for it in xrange(0, nt_ds, n):
            rf_pipelines.apply_std_dev_clipper(i_copy[:,(it):(it+n)], w_copy[:,(it):(it+n)], axis, sigma, Df, Dt, two_pass)
    
        return (i_copy, w_copy)


    if pipeline_json['class_name'] == 'mask_serializer':
        assert have_hdf5
        assert params.nds == 1

        hdf5_filename = pipeline_json['hdf5_filename']        

        f = h5py.File(hdf5_filename, 'r')
        assert f.attrs['nfreq'] == params.nfreq_ds
        assert f.attrs['initial_fpga_count'] == params.initial_fpga_count
        assert f.attrs['fpga_counts_per_sample'] == params.fpga_counts_per_sample

        m8 = f['bitmask'][:,:]
        f.close()

        # Temporary hdf5 file gets deleted here (see comment at the top of
        # make_random_mask_serializer())
        os.remove(hdf5_filename)

        assert m8.ndim == 2
        assert m8.shape[0] == weights.shape[0]
        assert m8.dtype == np.uint8

        mbool1 = np.zeros((m8.shape[0], m8.shape[1]*8), dtype=np.bool)
        for i in xrange(8):
            mbool1[:,i::8] = (m8[:,:] & (np.uint8(1) << i))

        mbool2 = (weights > 0.0)

        # Compare mbool1 and mbool2, but first discard trailing zeros.

        for n1 in xrange(mbool1.shape[1], -1, -1):
            if (n1 > 0) and np.any(mbool1[:,n1-1]):
                break

        for n2 in xrange(mbool2.shape[1], -1, -1):
            if (n2 > 0) and np.any(mbool2[:,n2-1]):
                break

        assert n1 == n2
        assert np.array_equal(mbool1[:,:n1], mbool2[:,:n2])
        # print 'mask_serializer debug:', n1, np.count_nonzero(mbool1) / float(mbool1.size)

        return (intensity, weights)


    if pipeline_json['class_name'] == 'mask_deserializer':
        assert have_hdf5
        assert params.nds == 1
        assert params.nt_tot % 8 == 0
        assert weights.shape[0] == params.nfreq_ds
        assert weights.shape[1] >= params.nt_tot

        if weights.shape[1] > params.nt_tot:
            assert np.all(weights[:,params.nt_tot:] == 0.0)

        hdf5_filename = pipeline_json['hdf5_filename']

        f = h5py.File(hdf5_filename, 'r')
        assert f.attrs['fpga_counts_per_sample'] == params.fpga_counts_per_sample

        m8 = f['bitmask'][:,:]
        assert m8.ndim == 2
        assert m8.shape[0] == params.nfreq_ds
        assert m8.dtype == np.uint8

        file_i0 = xdiv(f.attrs['initial_fpga_count'], params.fpga_counts_per_sample)
        file_i1 = i0 + (8 * m8.shape[1])
        f.close()

        # Temporary hdf5 file gets deleted here (see comment at the top of
        # make_random_mask_deserializer())
        os.remove(hdf5_filename)

        data_i0 = xdiv(params.initial_fpga_count, params.fpga_counts_per_sample)

        assert data_i0 % 8 == 0
        assert file_i0 % 8 == 0
        assert file_i0 <= data_i0
        assert file_i1 >= data_i0 + params.nt_tot

        # Shrink m8 to shape (nfreq, nt_tot/8).  Data type is still uint8.
        i0 = xdiv(data_i0 - file_i0, 8)
        i1 = i0 + xdiv(params.nt_tot, 8)
        m8 = m8[:,i0:i1]

        mbool = np.zeros((params.nfreq_ds, params.nt_tot), dtype=np.bool)
        for i in xrange(8):
            mbool[:,i::8] = (m8[:,:] & (np.uint8(1) << i))
        
        wout = weights.copy()
        wout[:,:params.nt_tot] *= mbool

        return (intensity, wcopy)


    raise RuntimeError('emulate_pipeline: unsupported class_name "%s"' % pipeline_json['class_name'])


####################################################################################################


class initial_stream(rf_pipelines.wi_stream):
    def __init__(self, p):
        assert isinstance(p, pipeline_params)
        assert p.nds == 1

        rf_pipelines.wi_stream.__init__(self, 'initial_stream')
        shape = (p.nfreq_ds, p.nt_tot)

        self.params = p
        self.nfreq = p.nfreq_ds
        self.nds = 1
        self.intensity = rand.standard_normal(size=shape)
        self.weights = rand.uniform(0.5, 1.0, size=shape)
        self.weights *= rand.choice([0.0,1.0], p=[0.1,0.9], size=shape)   # randomly mask 10%
        self.nt_chunk = rand.randint(20, 50)


    def _start_pipeline(self, json_attrs):
        # Catering to mask_serializer.
        json_attrs['initial_fpga_count'] = self.params.initial_fpga_count
        json_attrs['fpga_counts_per_sample'] = self.params.fpga_counts_per_sample


    def _fill_chunk(self, intensity, weights, pos):
        nt_tot = self.params.nt_tot
        intensity[:,:] = 0.
        weights[:,:] = 0.

        if pos >= nt_tot:
            return False

        n = min(nt_tot - pos, self.nt_chunk)
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
    # Make random pipeline params
    params = pipeline_params()
    nfreq = params.nfreq_ds
    nt_tot = params.nt_tot

    s = initial_stream(params)
    tj = make_random_transform_list(params, nelements=rand.randint(10,20))
    t = [ rf_pipelines.pipeline_object.from_json(j) for j in tj ]
    u = final_transform()

    p = rf_pipelines.pipeline([s] + t + [u])
    p.bind(outdir=None, verbosity=0, debug=True)

    # Check jsonization (test is slightly stronger if bind() comes first)
    tj2 = [ x.jsonize() for x in t ]
    rf_pipelines.utils.json_assert_equal(tj, tj2, name1='reference_json', name2='pipeline_json')

    # First run
    p.run(outdir=None, verbosity=0, debug=True)
    (i0,w0) = u.get_results()

    nt0 = i0.shape[1]
    assert nt0 >= nt_tot
    assert i0.shape == w0.shape == (nfreq, nt0)
    assert np.all(w0[:,nt_tot:] == 0.0)
    
    # Second run
    pj = { 'class_name': 'pipeline', 'elements': tj }
    (i1,w1) = emulate_pipeline(pj, s.intensity, s.weights, params)
    
    nt1 = i1.shape[1]
    assert nt1 >= nt_tot
    assert i1.shape == w1.shape == (nfreq, nt1)
    assert np.all(w1[:,nt_tot:] == 0.0)

    # Compare
    eps_i = maxdiff((i0*w0)[:,:nt_tot], (i1*w1)[:,:nt_tot])
    eps_w = maxdiff(w0[:,:nt_tot], w1[:,:nt_tot])

    assert eps_i < 1.0e-5
    assert eps_w < 1.0e-5


####################################################################################################

rand.seed(23)
niter = 100

def delete_temporary_files():
    for f in glob.glob('tmp*.h5'):
        print 'deleting leftover temporary file %s' % f
        os.remove(f)

for iter in xrange(niter):
    delete_temporary_files()
    print 'iteration %d/%d' % (iter, niter)
    run_test()

delete_temporary_files()
print 'test-almost-everything: pass'
