- version 4

  - plotter_transform: fix bug which sometimes caused outlier intensities to be superfluously masked
  
  - plotter_transform: determine color mapping using outlier-clipped mean/variance (without this, the
    plots were sometimes misleading!)

- version 3

  - From Masoud: New transform badchannel_mask, which reads a text file containing a list of bad 
    frequency ranges, and masks channels which overlap (currently python-only)

- version 2

  - new transform chime_file_writer, which writes a stream to disk in CHIME hdf5 format.

  - fix end-of-stream bug in plotter transform
