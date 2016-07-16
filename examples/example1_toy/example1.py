#!/usr/bin/env python

import os
import sys
import numpy as np
import rf_pipelines


if not os.path.exists('bonsai_config.hdf5'):
    print "Before running this script, you need to create the file 'bonsai_config.hdf5', using this command:"
    print "  bonsai-mkweight bonsai_config.txt bonsai_config.hdf5"
    sys.exit(1)


# To illustrate writing a transform in python, here is a toy transform
# which injects a sinusoidal signal, and also masks sporadic blocks of data.
#
# This code will be much clearer if read in conjunction with the rf_pipelines docstrings,
# especially the docstring for 'class py_wi_transform' in rf_pipelines/__init__.py !
#
# See the README file in this directory for instructions for running this script.
#
# For more documentation of the rf_transform API and builtin stream/transforms, see
# the python docstrings.


class toy_transform(rf_pipelines.py_wi_transform):
    def __init__(self, nt_chunk):
        # self.nt_chunk defines the chunk size for passing data into the transform.
        # As explained in the py_wi_transform docstring, each transform can choose 
        # this independently of the other transforms.
        self.nt_chunk = nt_chunk
        self.nt_prepad = 0
        self.nt_postpad = 0


    def set_stream(self, stream):
        # As explained in the py_wi_transform docstring, the four members
        #   self.nfreq, self.nt_chunk, self.nt_prepad, self.nt_postpad
        # must all be initialized either in the transform constructor or in set_stream().
        # 
        # In our toy transform, we initialize nt_chunk, nt_prepad, nt_postpad in
        # the constructor, but defer initialization of nfreq to set_stream(), where
        # it can be conveniently initialized from the stream.

        self.nfreq = stream.nfreq


    def process_chunk(self, t0, t1, intensity, weights, pp_intensity, pp_weights):
        # The API for process_chunk() is in the docstring for 'class py_wi_transform'
        # in rf_pipelines/__init__.py.  One thing to note here!  The frequency channels
        # are always ordered from highest frequency to lowest.  This is the same as
        # the 'bonsai' ordering.

        assert intensity.shape == weights.shape == (self.nfreq, self.nt_chunk)

        # These two lines add a fake intensity signal which is sinusodial in both
        # frequency and time (the details here are completely arbitrary!)
        (x, y) = np.meshgrid(np.linspace(0,1,self.nfreq), np.linspace(0,1,self.nt_chunk), indexing='ij')
        intensity[:,:] += 3. * (np.sin(np.pi*x) * np.sin(np.pi*y))**2
        
        # These three lines mask an arbitrarily chosen block of (frequency, time) pixels.
        # The masked block is repeated every time a length-nt_chunk block of data is
        # processed by process_chunk(), so you'll see a repeating periodic mask in the plots.
        (ifreq, jfreq) = (int(0.2*self.nfreq), int(0.3*self.nfreq))
        (it, jt) = (int(0.4*self.nfreq), int(0.5*self.nfreq))
        weights[ifreq:jfreq,it:jt] = 0.
        
# To run an rf_pipeline, we need an input stream and a sequence of transforms.
# For the stream,  we use the 'gaussian_noise_stream', which just simulates
# an independent random Gaussian for each (freq channel, time sample) pair.

s = rf_pipelines.gaussian_noise_stream(nfreq = 512,           # Number of frequency channels
                                       nt_tot = 12000,        # Length of stream in samples
                                       freq_lo_MHz = 400.0,   # Lower limit of band
                                       freq_hi_MHz = 800.0,   # Upper limit of band
                                       dt_sample= 1.0e-3)     # Sample duration in sec

# We define a sequence of transforms as follows.
# First, the 'toy_transform', which injects a fake sinusoidal signal as well as some masks
t1 = toy_transform(nt_chunk=512)

# Now inject a simulated FRB!
# We use signal-to-noise 100 so that it will be prominent by eye in plots.
# The undispersed_arrival_time is in seconds.  For this 
t2 = rf_pipelines.frb_injector_transform(snr=100, undispersed_arrival_time=5, dm=50)

# The plotter_transform can be put anywhere in the transform chain and is very convenient 
# for inspecting the pipeline.  (You can include multiple plotter_transforms, just use different 
# filenames of course!).  You'll be able to see the simulated FRB in 'waterfall_5.png'

t3 = rf_pipelines.plotter_transform(img_prefix = 'waterfall',   # filenames will be waterfall_0.png, ...
                                    img_nfreq = 512,     # no downsampling in freq axis, since img_nfreq = stream_nfreq
                                    img_nt = 1024,       # output images will be (512, 1024)
                                    downsample_nt = 1)   # no downsampling in time axis either

# The last "transform" doesn't actually modify the data, it just runs incoherent dedispersion.
# You'll need to create the config file with: 'bonsai-mkweight bonsai_config.txt'
# Coarse-grained triggers are written to 'triggers.hdf5', which can be plotted with 'bonsai-plot-triggers.py triggers.hdf5'.
# You'll see some wide diagonal stripes from the toy_transform output, plus a sharply peaked bowtie from the FRB.
t4 = rf_pipelines.bonsai_dedisperser('bonsai_config.hdf5', 'triggers.hdf5')

# Run the rf_pipeline!
s.run([t1,t2,t3,t4])

print "example1.py completed successfully"
print "You can plot the bonsai triggers with 'bonsai-plot-triggers.py triggers.hdf5'"
