- Version 6:

     - From Masoud: New transform legendre_detrender (polynomial fitting along either time or frequency axis)

     - From Masoud: New transform clipper_transform (clips outlier intensities, relative to a variance
       estimate which can be computed either in 2D, in 1D along the time axis, or in 1D along the frequency
       axis)

- Version 5:

     - The bonsai trigger hdf5 file output can now be generated in multifile mode.  Previously
       all triggers were written to one "monster file" of unbounded size.  Multifile mode is
       needed for long-running network streams, and also makes the rf_pipelines trigger plots
       easier to interpret, since they can be generated in 1-1 correspondence with intensity
       waterfall plots.  (See example3_chime for an example)

- Version 4

  - plotter_transform: fix bug which sometimes caused outlier intensities to be superfluously masked
  
  - plotter_transform: determine color mapping using outlier-clipped mean/variance (without this, the
    plots were sometimes misleading!)

- Version 3

  - From Masoud: New transform badchannel_mask, which reads a text file containing a list of bad 
    frequency ranges, and masks channels which overlap (currently python-only)

- Version 2

  - new transform chime_file_writer, which writes a stream to disk in CHIME hdf5 format.

  - fix end-of-stream bug in plotter transform
