#!/usr/bin/env python

import os
import rf_pipelines

# Analyze four arbitrarily chosen files
filename_list = [ '00003014.h5', '00003035.h5', '00003056.h5', '00003078.h5' ]
filename_list = [ os.path.join('/data/pathfinder/16-06-21',f) for f in filename_list ]

# Construct stream object.  Note: 'rf_pipelines' also contains stream constructors
# 'make_chime_stream_from_filename', 'make_chime_stream_from_acqdir'.  The last one
# should be used with caution, e.g. if pointed to /data/pathfinder/16-06-21 then it
# will analyze all 43GB of pathfinder data as a single stream!
# s = rf_pipelines.make_chime_stream_from_filename_list(filename_list)
s = rf_pipelines.make_chime_stream_from_acqdir('/data/pathfinder/16-06-21')

# Currently 'rf_pipelines' contains a simple detrending transform, and no RFI-removing transforms!
t1 = rf_pipelines.make_simple_detrender(8192)

# We put the plotter_transform after the detrender, so that it incrementally generates
# plots of detrended data.
t2 = rf_pipelines.plotter_transform('intensity', img_nfreq=512, img_nt=1024, downsample_nt=8)

# Dedisperse and write coarse-grained triggers to the file 'triggers.hdf5'.
#
# The bonsai_chime.hdf5 input file can be made with 'bonsai-mkweight'.
#
# The triggers.hdf5 file can be plotted with 'bonsai-plot-triggers.py' (warning: this script needs 
# improvement, in particular if run on a large stream it will make a "monster" plot with a huge
# number of pixels).

t3 = rf_pipelines.make_bonsai_dedisperser('bonsai_chime.hdf5', 'triggers.hdf5')

s.run([t1,t2,t3])