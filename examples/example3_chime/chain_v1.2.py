#!/usr/bin/env python
import os
import sys
import glob
import rf_pipelines

if not os.path.exists('bonsai_config.hdf5'):
    print "Before running this script, you need to create the file 'bonsai_config.hdf5', using this command:"
    print "  bonsai-mkweight bonsai_config.txt bonsai_config.hdf5"
    sys.exit(1)

# Note: The utility 'ch-show-intensity-file' may be useful for quickly inspecting a CHIME hdf5 file.
filename_list = [ '00000327.h5', '00000344.h5' ]
filename_list = [ os.path.join('/data/pathfinder/16-09-19-incoherent-without-noise-source',f) for f in filename_list ]

dsample_nt = 1024 # standard downsampling nt for 'master_clipper'
rms_cut = 0.002 # if the rms of intensity is above this threshold, then don't apply clippers AND set weights to 0  
detrend_nt = 2048 # standard nt for 'legendre_detrender'
niter = 2 # number of chain iterations (not counting internal iterations in 'master_clipper')

s = rf_pipelines.chime_stream_from_filename_list(filename_list, nt_chunk=1024, noise_source_align=detrend_nt)

transform_chain = [ rf_pipelines.plotter_transform('raw', img_nfreq=512, img_nt=1200, downsample_nt=16),
                    rf_pipelines.badchannel_mask('/data/pathfinder/rfi_masks/rfi_20160705.dat', nt_chunk=dsample_nt) ] 

for ix in xrange(niter):
    transform_chain += [ rf_pipelines.master_clipper(nt_chunk=1024, dsample_nt=dsample_nt, rms_cut=rms_cut, max_niter=4),
                         rf_pipelines.legendre_detrender(deg=4, axis=1, nt_chunk=detrend_nt), 
                         rf_pipelines.legendre_detrender(deg=8, axis=0, nt_chunk=detrend_nt), 
                         rf_pipelines.bonsai_dedisperser('bonsai_config.hdf5', 'triggers%d.hdf5' % ix, nt_per_file=16*1200),
                         rf_pipelines.plotter_transform('dc_out%d' % ix, img_nfreq=512, img_nt=1200, downsample_nt=16) ]

s.run(transform_chain)

print "chain_v1.2.py completed successfully"
print "You can plot the bonsai triggers with 'bonsai-plot-triggers.py triggers1.hdf5'"
