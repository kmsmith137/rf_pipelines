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
rms_cut = 0.0
mask_cut = 0.0
detrend_nt = 2048

s = rf_pipelines.chime_stream_from_filename_list(filename_list, nt_chunk=1024, noise_source_align=detrend_nt)

transform_chain = [ rf_pipelines.plotter_transform('raw', img_nfreq=512, img_nt=1200, downsample_nt=16),
                    rf_pipelines.badchannel_mask('/data/pathfinder/rfi_masks/rfi_20160705.dat', nt_chunk=clip_nt),
                    
                    rf_pipelines.intensity_clipper(thr=3, n_internal=12, axis=None, nt_chunk=clip_nt, dsample_nfreq=512, dsample_nt=clip_nt/16, imitate_cpp=True),
                    rf_pipelines.intensity_clipper(thr=3, n_internal=1, axis=0, nt_chunk=clip_nt, dsample_nfreq=1024, dsample_nt=clip_nt, imitate_cpp=True),
                    rf_pipelines.intensity_clipper(thr=3, n_internal=1, axis=1, nt_chunk=clip_nt, dsample_nfreq=1024, dsample_nt=clip_nt, imitate_cpp=True),
                    rf_pipelines.std_dev_clipper(thr=3, axis=1, dsample_nt=clip_nt/16, imitate_cpp=True),

                    rf_pipelines.intensity_clipper(thr=3, n_internal=12, axis=None, nt_chunk=clip_nt, dsample_nfreq=512, dsample_nt=clip_nt/16, imitate_cpp=True),
                    rf_pipelines.intensity_clipper(thr=3, n_internal=1, axis=0, nt_chunk=clip_nt, dsample_nfreq=1024, dsample_nt=clip_nt, imitate_cpp=True),
                    rf_pipelines.intensity_clipper(thr=3, n_internal=1, axis=1, nt_chunk=clip_nt, dsample_nfreq=1024, dsample_nt=clip_nt, imitate_cpp=True),
                    rf_pipelines.std_dev_clipper(thr=3, axis=1, dsample_nt=clip_nt/16, imitate_cpp=True),

                    rf_pipelines.intensity_clipper(thr=3, n_internal=12, axis=None, nt_chunk=clip_nt, dsample_nfreq=512, dsample_nt=clip_nt/16, imitate_cpp=True),
                    rf_pipelines.intensity_clipper(thr=3, n_internal=1, axis=0, nt_chunk=clip_nt, dsample_nfreq=1024, dsample_nt=clip_nt, imitate_cpp=True),
                    rf_pipelines.intensity_clipper(thr=3, n_internal=1, axis=1, nt_chunk=clip_nt, dsample_nfreq=1024, dsample_nt=clip_nt, imitate_cpp=True),
                    rf_pipelines.std_dev_clipper(thr=3, axis=1, dsample_nt=clip_nt/16, imitate_cpp=True),

                    rf_pipelines.polynomial_detrender(deg=4, axis=1, nt_chunk=detrend_nt),
                    rf_pipelines.polynomial_detrender(deg=8, axis=0, nt_chunk=detrend_nt),
                    rf_pipelines.plotter_transform('dc_out%d' % 1, img_nfreq=512, img_nt=1200, downsample_nt=16),
                    
                    rf_pipelines.intensity_clipper(thr=3, n_internal=12, axis=None, nt_chunk=clip_nt, dsample_nfreq=512, dsample_nt=clip_nt/16, imitate_cpp=True),
                    rf_pipelines.intensity_clipper(thr=3, n_internal=1, axis=0, nt_chunk=clip_nt, dsample_nfreq=1024, dsample_nt=clip_nt, imitate_cpp=True),
                    rf_pipelines.intensity_clipper(thr=3, n_internal=1, axis=1, nt_chunk=clip_nt, dsample_nfreq=1024, dsample_nt=clip_nt, imitate_cpp=True),
                    rf_pipelines.std_dev_clipper(thr=3, axis=1, dsample_nt=clip_nt/16, imitate_cpp=True),

                    rf_pipelines.intensity_clipper(thr=3, n_internal=12, axis=None, nt_chunk=clip_nt, dsample_nfreq=512, dsample_nt=clip_nt/16, imitate_cpp=True),
                    rf_pipelines.intensity_clipper(thr=3, n_internal=1, axis=0, nt_chunk=clip_nt, dsample_nfreq=1024, dsample_nt=clip_nt, imitate_cpp=True),
                    rf_pipelines.intensity_clipper(thr=3, n_internal=1, axis=1, nt_chunk=clip_nt, dsample_nfreq=1024, dsample_nt=clip_nt, imitate_cpp=True),
                    rf_pipelines.std_dev_clipper(thr=3, axis=1, dsample_nt=clip_nt/16, imitate_cpp=True),

                    rf_pipelines.intensity_clipper(thr=3, n_internal=12, axis=None, nt_chunk=clip_nt, dsample_nfreq=512, dsample_nt=clip_nt/16, imitate_cpp=True),
                    rf_pipelines.intensity_clipper(thr=3, n_internal=1, axis=0, nt_chunk=clip_nt, dsample_nfreq=1024, dsample_nt=clip_nt, imitate_cpp=True),
                    rf_pipelines.intensity_clipper(thr=3, n_internal=1, axis=1, nt_chunk=clip_nt, dsample_nfreq=1024, dsample_nt=clip_nt, imitate_cpp=True),
                    rf_pipelines.std_dev_clipper(thr=3, axis=1, dsample_nt=clip_nt/16, imitate_cpp=True),
                    
                    rf_pipelines.polynomial_detrender(deg=4, axis=1, nt_chunk=detrend_nt),
                    rf_pipelines.polynomial_detrender(deg=8, axis=0, nt_chunk=detrend_nt),
                    rf_pipelines.bonsai_dedisperser('bonsai_config.hdf5', 'triggers.hdf5', nt_per_file=16*1200),
                    rf_pipelines.plotter_transform('dc_out%d' % 2, img_nfreq=512, img_nt=1200, downsample_nt=16) ] 

s.run(transform_chain)

print "chain_v1.2.py completed successfully"
print "You can plot the bonsai triggers with 'bonsai-plot-triggers.py triggers.hdf5'"
