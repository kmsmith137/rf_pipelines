#!/usr/bin/env python

import os
import rf_pipelines

# Analyze four arbitrarily chosen files from the June 21 run.
# The utility 'ch-show-intensity-file' may be useful for quickly inspecting a CHIME hdf5 file.
filename_list = [ '00003014.h5', '00003035.h5', '00003056.h5', '00003078.h5' ]
filename_list = [ os.path.join('/data/pathfinder/16-06-21',f) for f in filename_list ]

# Construct CHIME stream object.
s = rf_pipelines.chime_stream_from_filename_list(filename_list)

# Currently 'rf_pipelines' contains a simple detrending transform, and no RFI-removing transforms!
t1 = rf_pipelines.simple_detrender(1024)

# We put the plotter_transform after the detrender, so that it incrementally generates
# plots of detrended data.  Downsampling by a factor 16 in time, and using 1200 coarse-grained
# times per waterfall plot, we get 4 waterfall plots in this stream.
t2 = rf_pipelines.plotter_transform('detrended_chime', img_nfreq=512, img_nt=1200, downsample_nt=16)

# Dedisperse and write coarse-grained triggers to the file 'triggers.hdf5'.
#
# The bonsai_config.hdf5 input file can be made with 'bonsai-mkweight'.
#
# The triggers.hdf5 file can be plotted with 'bonsai-plot-triggers.py' (warning: this script needs 
# improvement, in particular if run on a large stream it will make a "monster" plot with a huge
# number of pixels).

t3 = rf_pipelines.bonsai_dedisperser('bonsai_config.hdf5', 'triggers.hdf5')

s.run([t1,t2,t3])
