import sys
import numpy as np
import rf_pipelines


class plotter_transform(rf_pipelines.py_wi_transform):
    """
    This is a pseudo-transform (meaning that it does not actually modify its input)
    which makes waterfall plots.  It's primitive so please feel free to improve it!

    Note that the y-axis of the waterfall plots corresponds to frequency channel, with 
    lowest frequency on the bottom and highest frequency on the top.

    Constructor syntax:
    
      t = plotter_transform(img_prefix, img_nfreq, img_nt, downsample_nt=1, nt_chunk=0)
      
      The 'img_prefix' arg determines the output filenames: ${img_prefix}_0.png,
      ${img_prefix}_1.png ...

      The 'img_nfreq' arg determines the number of y-pixels in the plot.
      If the number of instrumental frequency channels is larger it will be downsampled.
    
      If the 'downsample_nt' optional arg is specfied, then each x-pixel in the plot will
      correspond to multiple time samples.

      The 'img_nt' arg determines the number of x-pixels before a waterfall plot gets
      written, and a new waterfall plot is started.
      
      By default, the color scheme is assigned by computing the mean and rms after clipping
      3-sigma outliers using three masking iterations.  The 'clip_niter' and 'sigma_clip' 
      arguments can be used to override these defaults.

      The 'nt_chunk' argument is the number of samples of data which are processed in
      each call to process_chunk().
    """

    def __init__(self, img_prefix, img_nfreq, img_nt, downsample_nt=1, nt_chunk=0, clip_niter=3, sigma_clip=3.0):
        # Call base class constructor
        rf_pipelines.py_wi_transform.__init__(self, 'plotter_transform')

        assert img_nt > 0
        assert img_nfreq > 0
        assert downsample_nt > 0
        assert clip_niter >= 1
        assert sigma_clip >= 2.0    # following assert in rf_pipelines.utils.write_png()

        if nt_chunk == 0:
            nt_chunk = min(img_nt, 1024//downsample_nt+1) * downsample_nt    # default value

        assert nt_chunk > 0

        if nt_chunk % downsample_nt != 0:
            raise RuntimeError("plotter_transform: current implementation requires 'nt_chunk' to be a multiple of 'downsample_nt'")

        # As explained in the rf_pipelines.py_wi_transform docstring, the following members
        # of the py_wi_transform base class must be initialized in the subclass.  (We also
        # need to initialize 'nfreq', but that will be done in 'set_stream', see below.)
        self.nt_chunk = nt_chunk
        self.nt_prepad = 0
        self.nt_postpad = 0

        self.img_prefix = img_prefix
        self.img_nfreq = img_nfreq
        self.img_nt = img_nt
        self.downsample_nt = downsample_nt
        self.nt_chunk_ds = nt_chunk // downsample_nt
        self.clip_niter = clip_niter
        self.sigma_clip = sigma_clip

        # Create plot_group 0.  See docstrings for py_wi_transform, py_wi_transform.add_plot_group(),
        # and py_wi_transform.add_plot() for more details.
        self.add_plot_group("waterfall", nt_per_pix=downsample_nt, ny=img_nfreq)


    def set_stream(self, s):
        # As explained in the rf_pipelines.py_wi_transform docstring, this function is
        # called once per pipeline run, when the stream (the 's' argument) has been specified.

        # Now that the stream has been specified, we can initialize self.nfreq as required.
        self.nfreq = s.nfreq

        if s.nfreq % self.img_nfreq != 0:
            raise RuntimeError("plotter_transform: current implementation requires 'img_nfreq' to be a divisor of stream nfreq")


    def start_substream(self, isubstream, t0):
        # Called once per substream (a stream can be split into multiple substreams).
        # We initialize per-substream state: a buffer for the downsampled data which
        # will be used to construct the plot.

        self.intensity_buf = np.zeros((self.img_nfreq,self.img_nt), dtype=np.float32)
        self.weight_buf = np.zeros((self.img_nfreq,self.img_nt), dtype=np.float32)
        self.isubstream = isubstream
        self.ifile = 0    # keeps track of which png file we're accumulating
        self.ipos = 0     # keeps track of how many (downsampled) time samples have been accumulated into file so far


    def process_chunk(self, t0, t1, intensity, weights, pp_intensity, pp_weights):
        # This is the main computational routine defining the transform, which is called
        # once per incoming "block" of data.  For documentation on the interface see
        # rf_pipelines.py_wi_transform docstring.

        # Invariant: When this routine is called, downsampled arrays of intensity and
        # weights data have been partially accumulated, in self.intensity_buf and self.weight_buf.
        # When the arrays are complete, an image will be written.  The number of files written
        # so far is 'self.ifile', and the number of (downsampled) time samples accumulated into
        # the current file is 'self.ipos'

        # Downsample data block being processed, before accumulating into arrays.
        (intensity, weights) = rf_pipelines.wi_downsample(intensity, weights, self.img_nfreq, self.nt_chunk_ds)

        # Keeps track of how much downsampled data has been moved from input chunk to the plots.
        ichunk = 0

        while ichunk < self.nt_chunk_ds:
            # Move to end of chunk or end of current plot, whichever comes first.
            n = min(self.nt_chunk_ds - ichunk, self.img_nt - self.ipos)
            assert n > 0
            
            self.intensity_buf[:,self.ipos:(self.ipos+n)] = intensity[:,ichunk:(ichunk+n)]
            self.weight_buf[:,self.ipos:(self.ipos+n)] = weights[:,ichunk:(ichunk+n)]
            self.ipos += n
            ichunk += n

            if self.ipos == self.img_nt:
                self._write_file()


    def end_substream(self):
        if self.ipos > 0:
            self._write_file()


    def _write_file(self):
        # When we reach end-of-stream, the buffer might be partially full (i.e. self.ipos < self.img_nt).
        # In this case, the plotting convention which I like best cosmetically is to truncate the image if there
        # is only one file in the output (i.e. self.ifile==0), otherwise pad with black.

        intensity = self.intensity_buf if (self.ifile > 0) else self.intensity_buf[:,:self.ipos]
        weights = self.weight_buf if (self.ifile > 0) else self.weight_buf[:,:self.ipos]

        basename = self.img_prefix
        if self.isubstream > 0:
            basename += str(isubstream+1)
        basename += ('_%s.png' % self.ifile)

        # The add_plot() method adds the plot to the JSON output, and returns the filename that should be written.
        # Note that a transform which writes multiple plot_groups would need to specify a group_id in add_plot().
        # (By default the group_id is zero.)

        filename = self.add_plot(basename, 
                                 it0 = self.ifile * self.img_nt * self.downsample_nt,
                                 nt = self.img_nt * self.downsample_nt,
                                 nx = intensity.shape[1], 
                                 ny = intensity.shape[0])

        rf_pipelines.write_png(filename, intensity, weights=weights, transpose=True, ytop_to_bottom=True, 
                               clip_niter=self.clip_niter, sigma_clip=self.sigma_clip)

        self.intensity_buf[:,:] = 0.
        self.weight_buf[:,:] = 0.
        self.ifile += 1
        self.ipos = 0
