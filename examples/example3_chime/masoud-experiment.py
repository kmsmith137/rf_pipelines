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

#
# Construct CHIME stream object.  
#
# We use the noise_source_align optional arg, which ensures that the noise source is aligned 
# with the detrender chunks, by discarding a small amount of initial data if necessary.  Note 
# that the value of 'noise_source_align' should be equal to the DETRENDER chunk size, not the
# chime_stream's nt_chunk.  (In this example the two happen to be equal.)
#
s = rf_pipelines.chime_stream_from_filename_list(filename_list, nt_chunk=1024, noise_source_align=1024)

frb = rf_pipelines.frb_injector_transform(snr=300, undispersed_arrival_time=10,\
        sample_rms=1, dm=50)

# This plotter_transform is before the detrender, so it generates "raw" (non-detrended)
# plots.  Downsampling by a factor 16 in time, and using 1200 coarse-grained times per
# waterfall plot, we end up with 4 plots (filenames raw_chime_0.png, raw_chime_1.png, ...)

t1 = rf_pipelines.plotter_transform('raw_chime', img_nfreq=512, img_nt=1200, downsample_nt=16)

# Mask out bad channels (i.e., weights[bad] = 0)
t2 = rf_pipelines.badchannel_mask('/data/pathfinder/rfi_masks/rfi_20160705.dat', nt_chunk=512)

# A very simple detrender (should implement something better soon!)
# The argument to the simple_detrender constructor is the detrender chunk size.
# The value of 'noise_source_align' above should be chosen to equal this.
t3 = rf_pipelines.simple_detrender(1024)

# -------------------------------------------------------------
# ------------------ Clipper Transforms -----------------------
# -------------------------------------------------------------
# Localize Horizontal RFIs: dsample along freq; clip along time
# Localize Vertical RFIs: dsample along time; clip along freq
# -------------------------------------------------------------
# Pre-fit clippers
c1 = rf_pipelines.clipper_transform(thr=3, axis=0, nt_chunk=1024,\
        dsample_nfreq=1024, dsample_nt=1024/256, test=False)

c2 = rf_pipelines.clipper_transform(thr=3, axis=2, nt_chunk=1024,\
        dsample_nfreq=1024, dsample_nt=1024/256)

c3 = rf_pipelines.clipper_transform(thr=3, axis=1, nt_chunk=1024,\
        dsample_nfreq=1024/128, dsample_nt=1024)

# Post-fit clippers
c4 = rf_pipelines.clipper_transform(thr=3, axis=1, nt_chunk=2048,\
        dsample_nfreq=1024/64, dsample_nt=1024)

c5 = rf_pipelines.clipper_transform(thr=3, axis=2, nt_chunk=1024,\
        dsample_nfreq=1024/4, dsample_nt=1024/256)

c6 = rf_pipelines.clipper_transform(thr=3, axis=0, nt_chunk=1024,\
        dsample_nfreq=1024, dsample_nt=4)

# Legendre detrenders
l1 = rf_pipelines.legendre_detrender(deg=4, axis=1, nt_chunk=1024, test=False)

l2 = rf_pipelines.legendre_detrender(deg=4, axis=1, nt_chunk=1024)
l3 = rf_pipelines.legendre_detrender(deg=10, axis=0, nt_chunk=1024)

l4 = rf_pipelines.legendre_detrender(deg=10, axis=0, nt_chunk=1024)

# This plotter_transform is after the detrender, so it generates detrended plots.
t4 = rf_pipelines.plotter_transform('detrended_chime', img_nfreq=512, img_nt=1200, downsample_nt=16)

# Dedisperse and write coarse-grained triggers to the file 'triggers.hdf5'.
#
# We run bonsai in multifile mode, with nt_per_file=16*1200 so that the bonsai output files will
# be in 1-1 correspondence with the plotter_transforms above.  This makes it easier to interpret
# the outputs!
#
# The bonsai_config.hdf5 input file can be made with 'bonsai-mkweight'.
# The triggers.hdf5 file can be plotted with 'bonsai-plot-triggers.py'.

t5 = rf_pipelines.bonsai_dedisperser('bonsai_config.hdf5', 'triggers.hdf5', nt_per_file=16*1200)

s.run([t1,t2,\
    l1,\
    c1,c2,c3,\
    l2,l3,\
    c4,c5,c6,l4,\
    t4,t5])

print "example3.py completed successfully"
print "You can plot the bonsai triggers with 'bonsai-plot-triggers.py triggers.hdf5'"

