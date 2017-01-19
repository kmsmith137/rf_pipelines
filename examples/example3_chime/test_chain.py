#!/usr/bin/env python
import os
import sys
import glob
import rf_pipelines

if not os.path.exists('bonsai_config.hdf5'):
    print "Before running this script, you need to create the file 'bonsai_config.hdf5', using this command:"
    print "  bonsai-mkweight bonsai_config.txt bonsai_config.hdf5"
    sys.exit(1)

filename_list = [ '00000327.h5', '00000344.h5' ]
filename_list = [ os.path.join('/data/pathfinder/16-09-19-incoherent-without-noise-source',f) for f in filename_list ]

# --------------
nt_chunk = 1024
nfreq = 1024
# --------------
detrend_nt = 1024
clip_nt = 1024
# --------------
rms_cut = 1e10
mask_cut = 0.0
max_niter = 3
test = True
# --------------

py_clippers = [ 'rf_pipelines.clip_fx(intensity, weights, thr=3, n_internal=12, axis=None, dsample_nfreq=%d/2, dsample_nt=%d/16)' % (nfreq, clip_nt),
                'rf_pipelines.clip_fx(intensity, weights, thr=3, n_internal=1, axis=0)',
                'rf_pipelines.clip_fx(intensity, weights, thr=3, n_internal=1, axis=1)',
                'rf_pipelines.filter_stdv(intensity, weights, thr=3, axis=1, dsample_nt=%d/16)' % clip_nt ]

imitate_cpp_clippers = [ 'rf_pipelines.clip_fx(intensity, weights, thr=3, n_internal=12, axis=None, dsample_nfreq=%d/2, dsample_nt=%d/16, imitate_cpp=True)' % (nfreq, clip_nt),
                         'rf_pipelines.clip_fx(intensity, weights, thr=3, n_internal=1, axis=0, imitate_cpp=True)',
                         'rf_pipelines.clip_fx(intensity, weights, thr=3, n_internal=1, axis=1, imitate_cpp=True)',
                         'rf_pipelines.filter_stdv(intensity, weights, thr=3, axis=1, dsample_nt=%d/16, imitate_cpp=True)' % clip_nt ]

cpp_clippers = [ 'rf_pipelines_c.apply_intensity_clipper(intensity, weights, sigma=3, niter=12, iter_sigma=3, axis=None, Df=2, Dt=16)',
                 'rf_pipelines_c.apply_intensity_clipper(intensity, weights, sigma=3, niter=1, iter_sigma=3, axis=0, Df=1, Dt=1)',
                 'rf_pipelines_c.apply_intensity_clipper(intensity, weights, sigma=3, niter=1, iter_sigma=3, axis=1, Df=1, Dt=1)',
                 'rf_pipelines_c.apply_std_dev_clipper(intensity, weights, sigma=3, axis=1, Dt=16)' ]

test_fdict = {'py' : py_clippers,
              'imitate_cpp' : imitate_cpp_clippers,
              'cpp' : cpp_clippers
             }

s = rf_pipelines.chime_stream_from_filename_list(filename_list, nt_chunk=nt_chunk, noise_source_align=detrend_nt)

transform_chain = [ rf_pipelines.badchannel_mask('/data/pathfinder/rfi_masks/rfi_20160705.dat', nt_chunk=clip_nt),
                    rf_pipelines.polynomial_detrender_cpp(polydeg=0, axis=1, nt_chunk=detrend_nt),
                    rf_pipelines.master_transform(nt_chunk=nt_chunk, fdict=test_fdict, rms_cut=rms_cut, mask_cut=mask_cut, max_niter=max_niter, test=test) ]

s.run(transform_chain)

print "test_chain.py completed successfully.\n ----->>>> test_results[key] = [unmasked(weights), np.mean(weights), np.std(weights)] <<<<-----"
