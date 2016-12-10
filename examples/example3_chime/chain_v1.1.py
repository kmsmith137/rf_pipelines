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
#filename_list = sorted(glob.glob('/data/pathfinder/16-09-19-incoherent-without-noise-source/*.h5'))

detrend_nt = 2048
dsample_nt = 1024
nc = 2
nd = 2

s = rf_pipelines.chime_stream_from_filename_list(filename_list, nt_chunk=1024, noise_source_align=detrend_nt)
#frb = rf_pipelines.frb_injector_transform(snr=1000, undispersed_arrival_time=1045.0, sample_rms=0.05, dm=500.)

transform_chain = [ frb, rf_pipelines.plotter_transform('raw', img_nfreq=512, img_nt=1200, downsample_nt=16),
                    rf_pipelines.badchannel_mask('/data/pathfinder/rfi_masks/rfi_20160705.dat', nt_chunk=clipper_nt) ] 

for jx in xrange(nd):
    for ix in xrange(nc):
        transform_chain += master_clipper(dsample_nt=dsample_nt)
    transform_chain += [ rf_pipelines.legendre_detrender(deg=4, axis=1, nt_chunk=detrend_nt), 
                         rf_pipelines.legendre_detrender(deg=8, axis=0, nt_chunk=detrend_nt), 
                         rf_pipelines.bonsai_dedisperser('bonsai_config.hdf5', 'triggers%d.hdf5' % jx, nt_per_file=16*1200),
                         rf_pipelines.plotter_transform('dc_out%d' % jx, img_nfreq=512, img_nt=1200, downsample_nt=16)]

s.run(transform_chain)

print "chain_v1.1.py completed successfully"
print "You can plot the bonsai triggers with 'bonsai-plot-triggers.py triggers1.hdf5'"
