#!/usr/bin/env python

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

# Analyze four arbitrarily chosen files from the run.  We now use the 16-07-08 run,
# which supersedes the earlier buggy runs.  The filenames below assume you're running 
# on chimer.physics.mcgill.ca.
#
# Note: The utility 'ch-show-intensity-file' may be useful for quickly inspecting a CHIME hdf5 file.
#       and the utility 'ch-show-intensity-file' makes a quick waterfall plot.

filename_list = [ '00000770.h5', '00000786.h5', '00000802.h5', '00000819.h5' ]
filename_list = [ os.path.join('/data/pathfinder/16-07-08',f) for f in filename_list ]

# Construct CHIME stream object.
s = rf_pipelines.chime_stream_from_filename_list(filename_list)

# This plotter_transform is before the detrender, so it generates "raw" (non-detrended)
# plots.  Downsampling by a factor 16 in time, and using 1200 coarse-grained times per
# waterfall plot, we end up with 4 plots (filenames raw_chime_0.png, raw_chime_1.png, ...)

t1 = rf_pipelines.plotter_transform('raw_chime', img_nfreq=512, img_nt=1200, downsample_nt=16)

# Mask out bad channels (i.e., weights[bad] = 0)
t2 = rf_pipelines.badchannel_mask('/data/pathfinder/rfi_masks/rfi_20160705.dat', nt_chunk=512)

# A very simple detrender (should implement something better soon!)
t3 = rf_pipelines.simple_detrender(1024)

# This plotter_transform is after the detrender, so it generates detrended plots.
t4 = rf_pipelines.plotter_transform('detrended_chime', img_nfreq=512, img_nt=1200, downsample_nt=16)

# Dedisperse and write coarse-grained triggers to the file 'triggers.hdf5'.
#
# The bonsai_config.hdf5 input file can be made with 'bonsai-mkweight'.
#
# The triggers.hdf5 file can be plotted with 'bonsai-plot-triggers.py' (warning: this script needs 
# improvement, in particular if run on a large stream it will make a "monster" plot with a huge
# number of pixels).

t5 = rf_pipelines.bonsai_dedisperser('bonsai_config.hdf5', 'triggers.hdf5')

s.run([t1,t2,t3,t4,t5])

print "example3.py completed successfully"
print "You can plot the bonsai triggers with 'bonsai-plot-triggers.py triggers.hdf5'"

