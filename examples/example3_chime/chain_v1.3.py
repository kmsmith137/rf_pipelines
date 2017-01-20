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

clip_nt = 1024
detrend_nt = 2048

clipper_niter = 3
detrender_niter = 2

s = rf_pipelines.chime_stream_from_filename_list(filename_list, nt_chunk=1024, noise_source_align=detrend_nt)

def clipper_chain(ix):
    return [ rf_pipelines.intensity_clipper_cpp(sigma=3, niter=12, iter_sigma=3, axis=None, nt_chunk=clip_nt, Df=2, Dt=16),
             rf_pipelines.intensity_clipper_cpp(sigma=3, niter=1, iter_sigma=3, axis=0, nt_chunk=clip_nt, Df=1, Dt=1),
             rf_pipelines.intensity_clipper_cpp(sigma=3, niter=1, iter_sigma=3, axis=1, nt_chunk=clip_nt, Df=1, Dt=1),
             rf_pipelines.std_dev_clipper_cpp(sigma=3, axis=1, Dt=16)
           ]

def detrender_chain(jx):
    return [ rf_pipelines.polynomial_detrender_cpp(polydeg=4, axis=1, nt_chunk=detrend_nt),
             rf_pipelines.polynomial_detrender_cpp(polydeg=8, axis=0, nt_chunk=detrend_nt),
             rf_pipelines.plotter_transform('dc_out%d' % jx , img_nfreq=512, img_nt=1200, downsample_nt=16),  
           ]

transform_chain = [ rf_pipelines.plotter_transform('raw', img_nfreq=512, img_nt=1200, downsample_nt=16),
                    rf_pipelines.badchannel_mask('/data/pathfinder/rfi_masks/rfi_20160705.dat', nt_chunk=clip_nt)
                  ]

for ix in xrange(detrender_niter):
    for jx in xrange(clipper_niter):
        transform_chain += clipper_chain(jx)
    transform_chain += detrender_chain(ix)

transform_chain += [ rf_pipelines.bonsai_dedisperser('bonsai_config.hdf5', 'triggers.hdf5', nt_per_file=16*1200) ]

s.run(transform_chain)

print "chain_v1.3.py completed successfully"
print "You can plot the bonsai triggers with 'bonsai-plot-triggers.py triggers.hdf5'"
