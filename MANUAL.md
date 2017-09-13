Some day, there will be systematic documentation here!
In the meantime, here are some notes on things that have changed recently.

The "sept7" branches contain some new features, which we want to
merge to master branches soon.  In the meantime, you should either
use the "master branch everywhere setup":
```
simd_helpers  master
rf_kernels     [ not needed ]
bonsai        master
rf_pipelines  master
ch_frb_rfi    master
```
or the "sept7 branch everywhere" setup:
```
simd_helpers  sept7
rf_kernels    sept7
bonsai        sept7
rf_pipelines  sept7
ch_frb_rfi    sept7
```
Note that when switching between these two setups, you'll
want to rebuild every package from scratch (make clean; make -j all install)
**in the order shown above**.

Here are some recent changes to rf_pipelines transforms:

 - `spline_detrender(nt_chunk, axis, nbins, epsilon = 3.0e-4)`

   This is an experimental transform designed to address the high
   computational cost of the polynomial_detrender.
   
   A spline_detrender with N bins should be roughly equivalent to a
   polynomial_detrender with degree (2N+1), but its computational cost
   should be independent of N.
    
   The 'epsilon' parameter regulates the spline fit by penalizing large time derivatives.
   If epsilon is too small, then overfitting may occur in regions with sparse weights.
   If epsilon is too large, then the fitter may have difficulty "keeping up" with rapid
   variations in the data.  I think that 3.0e-4 (the default) is a reasonable choice of
   epsilon, but I haven't experimented systematically.
    
   FIXME: currently, the only allowed axis type is AXIS_FREQ (=0).
   I propose that we experiment with replacing the polynomial_detrender in the
   AXIS_FREQ case, and if this looks good, then I'll implement AXIS_TIME and AXIS_NONE.


 - `bonsai_dedisperser`

   The bonsai_dedisperser now has a boolean `fill_rfi_mask` flag.
   If set to true, then bonsai will do "online" variance estimation, and replace
   RFI-masked data with randomly simulated data (as we plan to do in the real-time
   CHIME search).  If you use the fill_rfi_mask flag, then you'll want to note
   the following!

   There are some new bonsai config_params, including one required
   parameter `variance_timescale`.  The parameter `reweighting_timescale`
   may also be worth playing with.  You'll want to read the section
   "Config_params which are used by bonsai's 'fill_rfi_mask' feature"
   in bonsai/MANUAL.md
   ([Link](https://github.com/CHIMEFRB/bonsai/blob/sept7/MANUAL.md#user-content-config-fill-rfi-mask),
   note that this is the "sept7" bonsai branch)

   Because there is a new required parameter, you'll need to update the
   bonsai config files.  You'll also want to re-generate the hdf5 files
   (using the `bonsai-mkweight` utility).  I did this for the "7-tree, 1024-freq"
   config, see these files on frb1:
   ```
   /data/bonsai_configs/bonsai_nfreq1024_7tree_v3.txt
   /data/bonsai_configs/bonsai_nfreq1024_7tree_v3.hdf5
   ```
   Note that `ch_frb_rfi.bonsai.nfreq1K_7tree(params, v=3)` will use the
   new config files.

 - I made some minor changes to ch_frb_rfi to support these new transforms
   (`transform_parameters.spline`, `transform_parameters.bonsai_fill_rfi_mask`).

   The new script `ch_frb_rfi/scripts/s1-kms.py` is an example which
   shows how to use the spline_detrender and the fill_rfi_mask features.
   By comparing with `ch_frb_rfi/scripts/s1.py`, you can see what needs
   to be changed.
