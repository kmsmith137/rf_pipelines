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
    
      t = plotter_transform(img_prefix, img_nfreq, img_nt, downsample_nt=1, n_zoom=4,  nt_chunk=0)
      
      The 'img_prefix' arg determines the output filenames: ${img_prefix}_0.png,
      ${img_prefix}_1.png ...

      The 'img_nfreq' arg determines the number of y-pixels in the plot.
      If the number of instrumental frequency channels is larger it will be downsampled.
    
      The 'img_nt' arg determines the number of x-pixels before a waterfall plot gets
      written, and a new waterfall plot is started.

      If the 'downsample_nt' optional arg is specfied, then each x-pixel in the plot will
      correspond to multiple time samples.

      The n_zoom arg determines how many zoom levels will be produced.

      The 'nt_chunk' argument is the number of samples of data which are processed in
      each call to process_chunk().

      By default, the color scheme is assigned by computing the mean and rms after clipping
      3-sigma outliers using three masking iterations.  The 'clip_niter' and 'sigma_clip' 
      arguments can be used to override these defaults.
    """

    def __init__(self, img_prefix, img_nfreq, img_nt, downsample_nt=1, n_zoom = 1, nt_chunk=0, clip_niter=3, sigma_clip=3.0):
        # Build up name string, showing arguments which differ from non-default values
        name = "plotter_transform('%s', img_nfreq=%d, img_nt=%d, downsample_nt=%d" % (img_prefix, img_nfreq, img_nt, downsample_nt)
        if nt_chunk > 0:
            name += ", nt_chunk=%d" % nt_chunk
        if clip_niter != 3:
            name += ", clip_niter=%d" % clip_niter
        if sigma_clip != 3.0:
            name += ", sigma_clip=%s" % sigma_clip
        name += ')'

        # Call base class constructor
        rf_pipelines.py_wi_transform.__init__(self, name)

        assert img_nt > 0
        assert img_nfreq > 0
        assert downsample_nt > 0
        assert clip_niter >= 1
        assert sigma_clip >= 2.0    # following assert in rf_pipelines.utils.write_png()
        assert n_zoom > 0
        
        max_downsample = downsample_nt * 2**(n_zoom-1)

        if nt_chunk == 0:
            # Choose default nt_chunk to be close to 1024, with added constraint of being evenly divisible by max_downsample.
            nt_chunk = (1024//max_downsample + 1) * max_downsample
        elif nt_chunk % max_downsample != 0:
            raise RuntimeError("plotter_transform: specified nt_chunk(=%d) must be a multiple of downsampling factor at max zoom level (=%d)" % (nt_chunk,max_downsample))

        self.nt_chunk = nt_chunk
        self.nt_prepad = 0
        self.nt_postpad = 0
        self.img_nfreq = img_nfreq
        self.img_nt = img_nt
        self.clip_niter = clip_niter
        self.sigma_clip = sigma_clip
        self.n_zoom = n_zoom
        
        # Parameters implemented as lists that change for each zoom level
        self.downsample_nt = [downsample_nt]
        self.nt_chunk_ds = [nt_chunk // downsample_nt]     # number of times/chunk divided by number of times per pixel -> pixels/chunk
        self.add_plot_group("waterfall", nt_per_pix=downsample_nt, ny=img_nfreq)
        self.img_prefix = [img_prefix + "_zoom0"]

        if self.n_zoom > 1:
            for zoom_level in xrange(self.n_zoom - 1):
                self.downsample_nt += [self.downsample_nt[zoom_level] * 2]   # zoom_level = previous element's index because of the original value added
                self.nt_chunk_ds += [nt_chunk // self.downsample_nt[zoom_level + 1]]
                self.add_plot_group("waterfall", nt_per_pix=self.downsample_nt[zoom_level + 1], ny=img_nfreq)
                self.img_prefix += [img_prefix + "_zoom" + str(zoom_level+1)]

 
    def set_stream(self, s):
        # As explained in the rf_pipelines.py_wi_transform docstring, this function is
        # called once per pipeline run, when the stream (the 's' argument) has been specified.
        self.nfreq = s.nfreq

        if s.nfreq % self.img_nfreq != 0:
            raise RuntimeError("plotter_transform: current implementation requires 'img_nfreq' to be a divisor of stream nfreq")


    def start_substream(self, isubstream, t0):
        # Called once per substream (a stream can be split into multiple substreams).
        # We initialize per-substream state: a buffer for the downsampled data which
        # will be used to construct the plot.

        # The buffers store intensity/weights data until a plot can be written out
        # This means that the buffers will have dimensions n_zoom x n_xpix x n_ypix
        self.intensity_buf = np.zeros((self.n_zoom, self.img_nfreq, self.img_nt), dtype=np.float32)
        self.weight_buf = np.zeros((self.n_zoom, self.img_nfreq, self.img_nt), dtype=np.float32)
        self.isubstream = isubstream
        self.ifile = np.zeros((self.n_zoom), int)    # keeps track of which png file we're accumulating
        self.ipos = np.zeros((self.n_zoom), int)     # keeps track of how many (downsampled) time samples have been accumulated into file so far


    def process_chunk(self, t0, t1, intensity, weights, pp_intensity, pp_weights):
        # This is the main computational routine defining the transform, which is called
        # once per incoming "block" of data.  For documentation on the interface see
        # rf_pipelines.py_wi_transform docstring.

        # Invariant: When this routine is called, downsampled arrays of intensity and
        # weights data have been partially accumulated, in self.intensity_buf and self.weight_buf.
        # When the arrays are complete, an image will be written.  The number of files written
        # so far is 'self.ifile', and the number of (downsampled) time samples accumulated into
        # the current file is 'self.ipos'

        preserved_intensity = intensity
        preserved_weights = weights

        for zoom_level in xrange(self.n_zoom):
            intensity = preserved_intensity.copy()
            weights = preserved_weights.copy()

            (intensity, weights) = rf_pipelines.wi_downsample(intensity, weights, self.img_nfreq, self.nt_chunk_ds[zoom_level])
            
            # Keeps track of how much downsampled data has been moved from input chunk to the plots.
            ichunk = 0

            while ichunk < self.nt_chunk_ds[zoom_level]:
                # Move to end of chunk or end of current plot, whichever comes first.
                n = min(self.nt_chunk_ds[zoom_level] - ichunk, self.img_nt - self.ipos[zoom_level])
                assert n > 0
            
                self.intensity_buf[zoom_level, :, self.ipos[zoom_level]:(self.ipos[zoom_level]+n)] = intensity[:,ichunk:(ichunk+n)]
                self.weight_buf[zoom_level, :, self.ipos[zoom_level]:(self.ipos[zoom_level]+n)] = weights[:,ichunk:(ichunk+n)]
                self.ipos[zoom_level] += n
                ichunk += n

                if self.ipos[zoom_level] == self.img_nt:
                    self._write_file(zoom_level)


    def end_substream(self):
        for zoom_level in xrange(self.n_zoom):
            if self.ipos[zoom_level] > 0:
                self._write_file(zoom_level)


    def _write_file(self, zoom_level):
        # When we reach end-of-stream, the buffer might be partially full (i.e. self.ipos < self.img_nt).
        # In this case, pad with black
        intensity = self.intensity_buf[zoom_level, :, :]
        weights = self.weight_buf[zoom_level, :, :]

        basename = self.img_prefix[zoom_level]
        if self.isubstream > 0:
            basename += str(isubstream+1)
        basename += ('_%s.png' % self.ifile[zoom_level])

        # The add_plot() method adds the plot to the JSON output, and returns the filename that should be written.
        filename = self.add_plot(basename, 
                                 it0 = int(self.ifile[zoom_level] * self.img_nt * self.downsample_nt[zoom_level]),
                                 nt = self.img_nt * self.downsample_nt[zoom_level],
                                 nx = intensity.shape[1], 
                                 ny = intensity.shape[0], 
                                 group_id = zoom_level)


        rf_pipelines.write_png(filename, intensity, weights=weights, transpose=True, ytop_to_bottom=True, 
                               clip_niter=self.clip_niter, sigma_clip=self.sigma_clip)

        self.intensity_buf[zoom_level, :, :] = 0.
        self.weight_buf[zoom_level, :, :] = 0.
        self.ifile[zoom_level] += 1
        self.ipos[zoom_level] = 0
