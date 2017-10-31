#!/usr/bin/env python

import os
import sys
import numpy as np
import rf_pipelines

# This stream will receive data over the network (source IP address doesn't need 
# to be specified; the stream will accept packets from any source).
s = rf_pipelines.chime_network_stream()

# This transform makes waterfall plots with the same parameters as the sending pipeline, for comparison afterwards.
t = rf_pipelines.plotter_transform(img_prefix='waterfall', img_nfreq=512, img_nt=1024, downsample_nt=1)

# Run the rf_pipeline.  We specify an outdir to avoid filename collisions with the sending pipeline.
s.run([t], outdir='receiving_pipeline')
