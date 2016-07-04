#!/usr/bin/env python
#
# This example is a little more streamlined than example1-toy.py, which includes more comments!

import rf_pipelines

# An example GBNCC file on chimer
# Note: 8192 frequency channels, sample length 8.192e-5 sec, total size is 120.1 sec (1466368 samples)
s = rf_pipelines.psrfits_stream('/data/pathfinder/gbncc_example/guppi_56990_GBNCC123399_0165_0001.fits')

# Detrending timescale is 8192 samples = 0.67 sec
t1 = rf_pipelines.simple_detrender(8192)

# We downsample by a factor 4 on the frequency axis, and a factor 1024 on the time axis.
# Thus each pixel on the time axis corresponds to 1024 samples = 0.084 sec.
#t2 = rf_pipelines.plotter_transform(img_prefix='detrended_gbncc', img_nfreq=512, img_nt=2048, downsample_nt=1024)

# Before running the dedisperser, you'll need to run the command 'bonsai-mkweight bonsai_config.txt bonsai_config.hdf5'
# which creates the config hdf5 file from the config text file.  The output file triggers.hdf5 contains coarse-grained
# triggers which can be plotted with 'bonsai-plot-triggers.py'.
# t3 = rf_pipelines.bonsai_dedisperser('bonsai_config.hdf5', 'triggers.hdf5')

s.run([t1])
