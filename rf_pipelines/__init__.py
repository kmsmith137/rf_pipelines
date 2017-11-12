"""
rf_pipelines: plugin-based radio astronomy pipelines.

Currently under construction and documentation is sparse!


Abstract base classes
---------------------
  pipeline_object
  chunked_pipeline_object
  wi_stream
  wi_transform

Container classes
-----------------
  pipeline
  wi_sub_pipeline

Streams
-------
  chime_stream_from_acqdir
  chime_stream_from_filename
  chime_stream_from_filename_list
  chime_frb_stream_from_filename
  chime_frb_stream_from_filename_list
  chime_frb_stream_from_glob
  chime_network_stream
  gaussian_noise_stream

Detrenders
----------
  spline_detrender
  polynomial_detrender

Clippers
--------
  intensity_clipper
  std_dev_clipper
  mask_expander

CHIME-specific
--------------
  chime_file_writer
  chime_packetizer
  chime_16k_spike_mask
  chime_16k_derippler
  chime_16k_stripe_analyzer

Miscellaneous transforms
------------------------
  adversarial_masker (*)
  badchannel_mask (*)
  bonsai_dedisperser (*)
  chime_file_writer (*)
  chime_packetizer (*)
  frb_injector_transform (*)
  mask_filler (*) 
  noise_filler (*)
  plotter_transform (*)
  variance_estimator (*)

Utilities
---------
  apply_polynomial_detrender
  apply_intensity_clipper
  apply_std_dev_clipper
  weighted_mean_and_rms
  wi_downsample

(*) = python-only
"""


import sys
import json
import numpy as np

# Not sure why, but this import has to be in the toplevel module,
# or png writing will sometimes segfault!
try:
    import PIL.Image
except:
    print >>sys.stderr, "rf_pipelines: import PIL.Image failed; many things will work but plotting will fail"

from .rf_pipelines_c import *

from .streams.chime_streams import chime_stream_from_times
from .transforms.adversarial_masker import adversarial_masker
#from .transforms.chime_packetizer import chime_packetizer
#from .transforms.chime_transforms import chime_file_writer
from .transforms.plotter_transform import plotter_transform
from .transforms.bonsai_dedisperser import bonsai_dedisperser
from .transforms.frb_injector_transform import frb_injector_transform
#from .transforms.badchannel_mask import badchannel_mask
#from .transforms.intensity_clipper import intensity_clipper
#from .transforms.polynomial_detrender import polynomial_detrender
#from .transforms.mask_expander import mask_expander
#from .transforms.kurtosis_filter import kurtosis_filter
#from .transforms.std_dev_clipper import std_dev_clipper
#from .transforms.thermal_noise_weight import thermal_noise_weight
#from .transforms.RC_detrender import RC_detrender
from .transforms.variance_estimator import variance_estimator
from .transforms.mask_filler import mask_filler
from .transforms.noise_filler import noise_filler
#from .transforms.online_mask_filler import online_mask_filler

# Helper routines for implementing new transforms in python.

from .utils import write_png, upsample, json_show, json_str, json_assert_equal, json_read, json_write, Variance_Estimates, run_for_web_viewer, run_anonymously

# Grouping code to handle output of bonsai_dedisperser
from L1b import L1Grouper
