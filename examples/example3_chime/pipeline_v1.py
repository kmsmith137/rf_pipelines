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
#       and the utility 'ch-show-intensity-file' makes a quick waterfall plot.

#filename_list = sorted(glob.glob('/data/pathfinder/16-09-19-incoherent-without-noise-source/*.h5'))[0:13]
filename_list = [ '00000327.h5', '00000344.h5', '00000360.h5', '00000376.h5', '00000393.h5', '00000409.h5', '00000426.h5' ]
filename_list = [ os.path.join('/data/pathfinder/16-09-19-incoherent-without-noise-source',f) for f in filename_list ]

detrend_nt = 2048
clipper_nt = 4096
niterations = 6

s = rf_pipelines.chime_stream_from_filename_list(filename_list, nt_chunk=1024, noise_source_align=detrend_nt)
#frb = rf_pipelines.frb_injector_transform(snr=30, undispersed_arrival_time=2105.0, sample_rms=0.005, dm=500.)

def make_dc_chain(ix):
    return [ rf_pipelines.legendre_detrender(deg=2, axis=1, nt_chunk=detrend_nt),
             rf_pipelines.legendre_detrender(deg=2, axis=0, nt_chunk=detrend_nt),
             rf_pipelines.clipper_transform(thr=3, axis=0, nt_chunk=clipper_nt, dsample_nfreq=1024/2, dsample_nt=clipper_nt/128),
             rf_pipelines.clipper_transform(thr=3, axis=0, nt_chunk=clipper_nt, dsample_nfreq=1024, dsample_nt=clipper_nt),
             rf_pipelines.clipper_transform(thr=3, axis=1, nt_chunk=clipper_nt, dsample_nfreq=1024, dsample_nt=clipper_nt),
             rf_pipelines.clipper_transform(thr=3, nt_chunk=clipper_nt, dsample_nfreq=1024/2, dsample_nt=clipper_nt/64),
             rf_pipelines.clipper_transform(thr=3, axis=1, nt_chunk=clipper_nt, dsample_nfreq=1024/128, dsample_nt=clipper_nt/4),
             rf_pipelines.mask_expander(thr=0.3, nt_chunk=clipper_nt/2**10),
             rf_pipelines.mask_expander(thr=0.3, nt_chunk=clipper_nt/2**8),
             rf_pipelines.mask_expander(thr=0.3, nt_chunk=clipper_nt/2**6),
             rf_pipelines.mask_expander(thr=0.3, nt_chunk=clipper_nt/2**4),
             rf_pipelines.legendre_detrender(deg=10, axis=0, nt_chunk=detrend_nt),
             rf_pipelines.plotter_transform('clipper_output%d' % ix, img_nfreq=512, img_nt=2400, downsample_nt=16) ]

transform_chain = [ rf_pipelines.plotter_transform('raw', img_nfreq=512, img_nt=2400, downsample_nt=16),
                    rf_pipelines.badchannel_mask('/data/pathfinder/rfi_masks/rfi_20160705.dat', nt_chunk=512) ]

for ix in xrange(niterations):
    transform_chain += make_dc_chain(ix)

transform_chain += [ rf_pipelines.bonsai_dedisperser('bonsai_config.hdf5', 'triggers.hdf5', nt_per_file=16*2400) ]

s.run(transform_chain)

print "kendrick-experiment.py completed successfully"
print "You can plot the bonsai triggers with 'bonsai-plot-triggers.py triggers.hdf5'"
