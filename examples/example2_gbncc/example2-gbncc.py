#!/usr/bin/env python
#
# See the README file in this directory for instructions for running this script.
#
# For more documentation of the rf_transform API and builtin stream/transforms, see
# the python docstrings.

import os
import sys
import rf_pipelines

if not os.path.exists('bonsai_config.hdf5'):
    print "Before running this script, you need to create the file 'bonsai_config.hdf5', using this command:"
    print "  bonsai-mkweight bonsai_config.txt bonsai_config.hdf5"
    sys.exit(1)

# An example GBNCC file on chimer
# Note: 8192 frequency channels, sample length 8.192e-5 sec, total size is 120.1 sec (1466368 samples)
s = rf_pipelines.psrfits_stream('/data/pathfinder/gbncc_example/guppi_56990_GBNCC123399_0165_0001.fits')

# Detrending timescale is 8192 samples = 0.67 sec
t1 = rf_pipelines.simple_detrender(8192)

# We put the plotter_transform after the detrender, so that it incrementally generates
# plots of detrended data.
#
# We downsample by a factor 4 on the frequency axis, and a factor 1024 on the time axis.
# Thus each pixel on the time axis corresponds to 1024 samples = 0.084 sec.
t2 = rf_pipelines.plotter_transform(img_prefix='detrended_gbncc', img_nfreq=512, img_nt=2048, downsample_nt=1024)

# Before running the dedisperser, you'll need to run the command 'bonsai-mkweight bonsai_config.txt bonsai_config.hdf5'
# which creates the config hdf5 file from the config text file.  The output file triggers.hdf5 contains coarse-grained
# triggers which can be plotted with 'bonsai-plot-triggers.py'.
#
# Note: the optional arugment 'nt_per_file' to rf_pipelines.bonsai_dedisperser() can be used to break the bonsai output
# into multiple files.  This can make the bonsai output files easier to interpret by putting them in 1-1 correspondence
# with the intensity waterfall plots.  See example3_chime for an explicit example.
t3 = rf_pipelines.bonsai_dedisperser('bonsai_config.hdf5', 'triggers.hdf5')

s.run([t1,t2,t3])

print "example2.py completed successfully"
print "You can plot the bonsai triggers with 'bonsai-plot-triggers.py triggers.hdf5'"
