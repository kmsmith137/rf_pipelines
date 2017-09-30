"""
rf_pipelines.retirement_home

This is where python pipeline_objects go to retire, after
they are replaced by their younger and more efficient C++
counterparts.

The python versions are still useful for reference,
or in unit tests which verify that fast C++ transforms 
and python reference transforms are equivalent.
"""

from .intensity_clipper import intensity_clipper, clip_fx
