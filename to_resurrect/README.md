At the end of September, the core rf_pipelines API was changed substantially!  ("rf_pipelines2")

This required revisiting every transform, and making minor changes.  This directory
contains transforms which have not yet been ported to the new API.

```
reverter.hpp
bitmask_maker.cpp
online_mask_filler.cpp
psrfits_stream.cpp
pulse_adder.cpp
reverter.cpp
test-mask-expander.py
rf_pipelines/transforms/kurtosis_filter.py
rf_pipelines/transforms/mask_expander.py
rf_pipelines/transforms/online_mask_filler.py
rf_pipelines/transforms/RC_detrender.py
rf_pipelines/transforms/thermal_noise_weight.py
```
