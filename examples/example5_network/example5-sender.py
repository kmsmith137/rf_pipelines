#!/usr/bin/env python
#
# This toy script simulates some data, using the same arbitrary procedure as in example1,
# makes waterfall plots, and sends the data over the network.


import os
import sys
import numpy as np
import rf_pipelines


# This toy transform is cut-and-paste from example1.  It adds a sinusoidal signal and
# masks some rectangular subarrays.  (No particular motivation for these steps, just
# results in a visually interesting plot.)
class toy_transform(rf_pipelines.py_wi_transform):
    def __init__(self, nt_chunk):
        self.name = 'toy_transform'
        self.nt_chunk = nt_chunk
        self.nt_prepad = 0
        self.nt_postpad = 0

    def set_stream(self, stream):
        self.nfreq = stream.nfreq

    def process_chunk(self, t0, t1, intensity, weights, pp_intensity, pp_weights):
        assert intensity.shape == weights.shape == (self.nfreq, self.nt_chunk)

        (x, y) = np.meshgrid(np.linspace(0,1,self.nfreq), np.linspace(0,1,self.nt_chunk), indexing='ij')
        intensity[:,:] += 3. * (np.sin(np.pi*x) * np.sin(np.pi*y))**2

        (ifreq, jfreq) = (int(0.2*self.nfreq), int(0.3*self.nfreq))
        (it, jt) = (int(0.4*self.nfreq), int(0.5*self.nfreq))
        weights[ifreq:jfreq,it:jt] = 0.


#
# The first step in the pipeline is a gaussian_noise stream, which just outputs Gaussian noise.
# NOTE two artificial restrictions of the current packet format:
#   - number of frequency channels must be a multiple of 1024
#   - dt_sample must be an integer multiple of 2.56e-6
#
s = rf_pipelines.gaussian_noise_stream(nfreq = 1024,                  # Number of frequency channels
                                       nt_tot = 12000,                # Length of stream in samples
                                       freq_lo_MHz = 400.0,           # Lower limit of band
                                       freq_hi_MHz = 800.0,           # Upper limit of band
                                       dt_sample = 384 * 2.56e-6)     # Sample duration in sec

# The toy_transform injects a fake sinusoidal signal as well as some masks
t1 = toy_transform(nt_chunk=512)

# The frb_injector_transform adds a simulated FRB.
t2 = rf_pipelines.frb_injector_transform(snr=100, undispersed_arrival_time=5, dm=50)

# This plotter_transform plots the data before sending over the network
t3 = rf_pipelines.plotter_transform(img_prefix = 'waterfall',   # filenames will be sender_waterfall_0.png, ...
                                    img_nfreq = 512,            # downsample by factor 2 in freq axis
                                    img_nt = 1024,              # output images will be (512, 1024)
                                    downsample_nt = 1)          # no downsampling in time axis either

# The packetizing_transform sends the data over the network.
# See chime_packetizer docstring for an explanation of all the arguments.
t4 = rf_pipelines.chime_packetizer(dstname = '127.0.0.1',
                                   nfreq_per_packet = 32,
                                   nt_per_chunk = 512,
                                   nt_per_packet = 32, 
                                   wt_cutoff = 0.5, 
                                   target_gbps = 0.1)

# Run the rf_pipeline.  We specify an outdir to avoid filename collisions with the receiving pipeline.
s.run([t1,t2,t3,t4], outdir='sending_pipeline')
