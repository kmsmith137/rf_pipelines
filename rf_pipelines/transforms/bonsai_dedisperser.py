import sys
import numpy as np
import rf_pipelines


class bonsai_dedisperser(rf_pipelines.py_wi_transform):
    """
    Returns a "transform" which doesn't actually modify the data, it just runs the bonsai dedisperser.


    Constructor arguments
    ---------------------

       - config_filename:  The configuration file used to initialize the dedisperser.  

       - img_prefix: Determines output filenames, using a similar convention to the plotter_transform.
           If img_prefix=None (the default), then no plots are generated.

       - img_ndm: Number of y-pixels in each output plot (same as img_nfreq in the plotter_transform).
           Note that the y-axis of the bonsai plots corresponds to dispersion measure.

       - img_nt: Number of x-pixels in each output plot.

       - downsample_nt: If > 1, then each x-pixel in the plot will correspond to multiple time samples.

       - n_zoom: Number of zoom levels in plots.

       - track_global_max: If True, then the following json output will be written:
             frb_global_max_trigger, frb_global_max_trigger_dm, frb_global_max_trigger_tfinal

       - dm_min, dm_max: Only meaningful if track_global_max is True.

       - hdf5_output_filename: If specified, HDF5 file(s) containing coarse-grained triggers will be written.

       - nt_per_hdf5_file: Only meaningful if hdf5_output_filename=True.  Zero means "one big file".

       - deallocate_between_substreams: infrequently-used option, used in frb_olympics

       - use_analytic_normalization: if True, then the dedisperser will use the exact trigger
           normalization, assuming a toy model in which each input (frequency, time) sample
           is an uncorrelated Gaussian.  Not suitable for real data!
    """

    def __init__(self, config_filename, img_prefix=None, img_ndm=256, img_nt=256, downsample_nt=1, n_zoom=1, 
                 track_global_max=False, dm_min=None, dm_max=None, hdf5_output_filename=None, nt_per_hdf5_file=0,
                 deallocate_between_substreams=False, use_analytic_normalization=False, dynamic_plotter=False,
                 plot_threshold1=6, plot_threshold2=10):

        # We import the bonsai module here, rather than at the top of the file, so that bonsai isn't
        # required to import rf_pipelines (but is required when you try to construct a bonsai_dedisperser).
        try:
            import bonsai
        except ImportError:
            raise RuntimeError("rf_pipelines: couldn't import the 'bonsai' module.  You may need to clone https://github.com/CHIMEFRB/bonsai and install.")
        
        if img_prefix is None:
            self.make_plot = False
        else:
            self.make_plot = True

        name = "bonsai_dedisperser('%s')" % config_filename
        rf_pipelines.py_wi_transform.__init__(self, name)

        self.config_filename = config_filename
        self.deallocate_between_substreams = deallocate_between_substreams
        
        initially_allocated = not deallocate_between_substreams
        self.dedisperser = bonsai.Dedisperser(config_filename, allocate=initially_allocated, use_analytic_normalization=use_analytic_normalization)
        self.global_max_tracker = None

        if track_global_max:
            self.global_max_tracker = bonsai.global_max_tracker(dm_min, dm_max)
            self.dedisperser.add_processor(self.global_max_tracker)

        if hdf5_output_filename:
            t = bonsai.trigger_hdf5_file_writer(hdf5_output_filename, nt_per_hdf5_file)
            self.dedisperser.add_processor(t)

        # Note that 'nfreq' is determined by the config file.  If the stream's 'nfreq' differs,
        # then an exception will be thrown.  The 'nt_chunk' parameter is also determined by the
        # config file, not a constructor argument.
        self.nfreq = self.dedisperser.nfreq
        self.nt_chunk = self.dedisperser.nt_chunk
        self.nt_prepad = 0
        self.nt_postpad = 0

        # Set plotting parameters
        if self.make_plot:
            self.dynamic_plotter = dynamic_plotter
            self.plot_threshold1 = plot_threshold1
            self.plot_threshold2 = plot_threshold2
            self.n_zoom = n_zoom
            self.img_ndm = img_ndm
            self.img_nt = img_nt
            self.downsample_nt = [downsample_nt]
            self.nt_chunk_ds = [self.nt_chunk // self.downsample_nt[0]]
            self.img_prefix = [str(img_prefix) + "_zoom0"]
            if self.n_zoom > 1:
                for zoom_level in xrange(self.n_zoom - 1):
                    self.downsample_nt += [self.downsample_nt[zoom_level] * 2]   # zoom_level = previous element's index because of the original value added
                    self.nt_chunk_ds += [self.nt_chunk // self.downsample_nt[zoom_level + 1]]
                    self.img_prefix += [img_prefix + "_zoom" + str(zoom_level+1)] 

            if self.nt_chunk % self.downsample_nt[-1] != 0:
                raise RuntimeError("bonsai plotter transform: specified nt_chunk(=%d) must be a multiple of downsampling factor at max zoom level (=%d)" 
                                   % (self.nt_chunk, self.downsample_nt[-1]))

            # Set incoming triger dimension paramaters for assertions
            # self.trigger_dim stores dimensions as a list of tuples containing dimensions for each tree
            # This is where our hard-coded tree logic begins, though it's easy enough to modify as desired
            self.trigger_dim = [(ndm, nt) for (ndm, nt) in zip(self.dedisperser.ndm_coarse, self.deidsperser.nt_coarse_per_chunk)]
            self.ntrees = len(self.trigger_dim)
            assert self.ntrees > 1
            assert self.trigger_dim[0][0]% self.img_ndm  == 0 or self.img_ndm % self.trigger_dim[0][0] == 0  # downsample or upsample dm for plot0
            assert all([tup[0] % self.img_ndm/self.ntrees  == 0 or self.img_ndm/self.ntrees % tup[0] == 0 for tup in self.trigger_dim])  # for plot1 
            assert all(tup[1] % (self.nt_chunk_ds[-1]) == 0 or self.nt_chunk_ds[0] % tup[1] == 0 for tup in self.trigger_dim)  # downsample or upsample t

            # Add a cheat-y assert, but I think this should be fine
            # Ensures that each chunk will evenly divide the plot, so no fancy buffering is required
            assert self.img_nt % self.nt_chunk_ds[-1] == 0

            # Convention for plot groups: 1st half will be the tree0 plot and 2nd half will be the plot for the rest of the trees
            for i in xrange(2*self.n_zoom):
                self.add_plot_group("waterfall", nt_per_pix=downsample_nt[i % (self.ntrees-1)], ny=img_ndm)


    def set_stream(self, stream):
        if stream.nfreq != self.nfreq:
            raise RuntimeError("rf_pipelines: number of frequencies in stream (nfreq=%d) does not match bonsai config file '%s' (nfreq=%d)" % (stream.nfreq, self.config_filename, self.nfreq))


    def start_substream(self, isubstream, t0):
        if self.deallocate_between_substreams:
            self.dedisperser.allocate()

        if self.make_plot:
            self.isubstream = isubstream
            plot0 = Plotter(self, ntrees=1, nt_chunk_ds=self.nt_chunk_ds[0], downsample_nt=self.downsample_nt[0], plot=0)
            plot1 = Plotter(self, ntrees=self.ntrees-1, nt_chunk_ds=self.nt_chunk_ds[1:], downsample_nt=self.downsample_nt[1:], plot=1)


    def process_chunk(self, t0, t1, intensity, weights, pp_intensity, pp_weights):
        # FIXME some day I'd like to remove this extra copy required by current cython implementation
        intensity = np.array(intensity, dtype=np.float32, order='C')
        weights = np.array(weights, dtype=np.float32, order='C')

        # Send the inputs (intensity, weights) to the dedisperser.
        self.dedisperser.run(intensity, weights)

        # Retrieve the outputs (trigger arrays) from the dedisperser.
        #
        # self.dedisperser.get_triggers() returns a list of 4D numpy arrays.
        #
        # Each array in the list corresponds to one dedispersion tree.  (The bonsai dedisperser
        # generally uses multiple trees internally, to dedisperse the data in different parts
        # of the (DM, pulse_width) parameter space.)
        #
        # Each 4D array is indexed by (DM_index, SM_index, beta_index, time_index).
        if self.make_plot:
            triggers = self.dedisperser.get_triggers()

            # Checking for invalid values in the incoming array
            if not all(np.all(np.isfinite(t)) for t in triggers):
                # Was the problem in the input arrays... ?
                if not np.all(np.isfinite(intensity)) or not np.all(np.isfinite(weights)):
                    raise RuntimeError('bonsai_dedisperser: input intensity/weights arrays contained Inf or NaN!')
                # If not, then maybe it's a 16-bit overflow... ?
                raise RuntimeError('bonsai returned Inf or NaN triggers!  Try reducing the normalization of the intensity and weights arrays, or setting nbits=32 in the bonsai config.')

            # Max along the beta and SM indices, then add the new triggers to the plots! 
            flat_triggers = [np.amax(np.amax(tree, axis=1), axis=1) for tree in triggers]
            plot0.process(flat_triggers[0])
            plot1.process(flat_triggers[1:])


    def end_substream(self):
        if self.make_plot:
            for zoom_level in xrange(self.n_zoom):
                if self.ipos[zoom_level] > 0:
                    self._write_file(zoom_level)

        self.dedisperser.end_dedispersion()
        
        if self.global_max_tracker is not None:
            # Add global max trigger data to pipeline json output
            self.json_per_substream["frb_global_max_trigger"] = self.global_max_tracker.global_max_trigger
            self.json_per_substream["frb_global_max_trigger_dm"] = self.global_max_tracker.global_max_trigger_dm
            self.json_per_substream["frb_global_max_trigger_tfinal"] = self.global_max_tracker.global_max_trigger_arrival_time

        if self.deallocate_between_substreams:
            self.dedisperser.deallocate()


    def _write_file(self, arr, zoom_level, downsample_nt, ifile, iplot):
        # When we reach end-of-stream, the buffer might be partially full. In this case, pad with black.
        basename = self.img_prefix[zoom_level] + '_plot' + str(iplot)
        if self.isubstream > 0:
            basename += str(isubstream+1)
        basename += ('_%s.png' % ifile)

        # The add_plot() method adds the plot to the JSON output, and returns the filename that should be written.                                                                                         
        filename = self.add_plot(basename,
                                 it0 = int(ifile[zoom_level] * self.img_nt * downsample_nt),
                                 nt = self.img_nt * downsample_nt, 
                                 nx = arr.shape[1],
                                 ny = arr.shape[0], 
                                 group_id = zoom_level*(iplot+1))
        
        if self.dynamic_plotter:
            rf_pipelines.write_png(filename, arr, transpose=True)
        else:
            rf_pipelines.utils.triggers_png(filename, arr, transpose=True, 
                                            threshold1=self.plot_threshold1, threshold2=self.plot_threshold2)


class Plotter():
    """A plotter object holds all desired zoom levels for a plot"""
    def __init__(self, transform, ntrees, nt_chunk_ds, downsample_nt, iplot):
        self.transform = transform       # To use _write_file()
        self.iplot = iplot               # For determining outfile name and which plot group a plot belongs to
        self.nzoom = transform.n_zoom    # Number of zoom levels to be produced
        self.nx = transform.img_nt       # Number of x pixels per plot
        self.ny = transform.img_ndm      # Number of y pixels per plot
        self.ntrees = ntrees             # Number of trees being __plotted together__ (different from transform.ntrees)
        self.nt_chunk_ds = nt_chunk_ds   # Number of x pixels that should be written per chunk (2D array)
        self.ndm = self.nx / self.ntrees # Number of y pixels per tree
        self.ix = np.zeros(self.nzoom)        # Keep track of what x position to add chunks to
        self.isubstream = transform.isubstream

        assert len(self.nt_chunk_ds) == self.ntrees
        assert all([len(zoom) == self.nzoom for zoom in self.nt_chunk_ds])
        assert np.all(np.sum(self.nt_chunk_ds, axis=1) == self.ny)

        self.plots = np.zeros((self.nzoom, self.ny, self.nx))        

        # Extra info for transform.write_file
        self.img_prefix = transform.img_prefix
        self.downsample_nt = downsample_nt  # In constructor so only relevant tree info passed (the same is true of nt_chunk_ds)
        self.ifile = np.zeros(nzoom)


        def process(self, arrs):
            # Check that the arrays passed to process contain the expected number of trees
            assert len(arrs) == self.ntrees
            
            # Zooming only happens in the time axis, so we can reshape the dm axis outside of the loop
            for i in xrange(self.ntrees):
                arr_shape = arrs[i].shape
                if arr_shape[0] > self.ndm:
                    arrs = self.max_ds(arrs[i], self.ndm, arr_shape[1])
                elif arr_shape[0] < self.ndm:
                    arrs = rf_pipelines.upsample(arrs[i], self.ndm, arr_shape[1])

            iy = 0  # Keep track of where to add each tree vertically
            for itree in xrange(self.ntrees):
                for zoom_level in xrange(self.nzoom):
                    dm_t = arrs[itree].copy()
                    dm_t_shape = dm_t.shape
                    # In the x (time) axis, we need to transform dm_t_shape[1] to self.nt_chunk_ds[zoom_level] - may need to up/downsample
                    if dm_t_shape[1] > self.nt_chunk_ds[zoom_level]:
                        dm_t = self.max_ds(dm_t, dm_t_shape[0], self.nt_chunk_ds[zoom_level])
                    elif dm_t_shape[1] < self.nt_chunk_ds[zoom_level]:
                        dm_t = rf_pipelines.upsample(dm_t, dm_t_shape[0], self.nt_chunk_ds[zoom_level])
                    # Now the array will be scaled properly to stick into the plot accumulator arrays
                    self.plot[zoom_level, self.iy:self.iy+self.ndm, self.ix[zoom_level]:self.ix[zoom_level]+self.nt_chunk_ds[zoom_level]] = dm_t
                    self.ix[zoom_level] += self.nt_chunk_ds[zoom_level]
                iy += self.ndm

            # After adding, check whether we need to write any files
            for zoom_level in xrange(self.nzoom):
                if self.ix[zoom_level] >= self.nx:
                    self.transform._write_file(self.plot[zoom_level, :, :], zoom_level, self.downsample_nt[zoom_level], self.ifile[zoom_level], self.iplot)
                    self.ifile[zoom_level] += 1
                    self.ix[zoom_level] = 0


    def max_ds(self, arr, new_dm, new_t):
        """Takes maxima along axes to downsample"""
        assert arr.ndim == 2
        assert new_dm > 0
        assert new_t > 0
        (ndm, nt) = arr.shape
        assert ndm % new_dm == 0
        assert nt % new_t == 0
        arr = np.reshape(arr, (new_dm, ndm//new_dm, new_t, nt//new_t))
        arr = np.amax(arr, axis=3)
        arr = np.amax(arr, axis=1)
        return arr
