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

       - dynamic_plotter, plot_threshold1, plot_threshold2: if False, plot_threshold1 and plot_threshold2 define
            the points of colour disontinuities in rf_pipelines.utils.triggers_png. If True, write_png is used
            and the plot threshold parameters are meaningless.

       - nplot_groups: The number of rows of plots the bonsai dedisperser should produce as displayed in the 
            web viewer. If 1, only tree 0 will be plotted. If >1, tree 0 will be plotted in the bottom row, and
            the remaining trees will be divided between the remaining plots. 

    """

    def __init__(self, config_filename, img_prefix=None, img_ndm=256, img_nt=256, downsample_nt=1, n_zoom=1, 
                 track_global_max=False, dm_min=None, dm_max=None, hdf5_output_filename=None, nt_per_hdf5_file=0,
                 deallocate_between_substreams=False, use_analytic_normalization=False, dynamic_plotter=False,
                 plot_threshold1=6, plot_threshold2=10, nplot_groups=1):

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

            self.trigger_dim = [(ndm, nt) for (ndm, nt) in zip(self.dedisperser.ndm_coarse, self.dedisperser.nt_coarse_per_chunk)]
            self.ntrees = len(self.trigger_dim)
            self.nplot_groups = nplot_groups

            assert self.ntrees > 0
            assert self.nplot_groups > 0
            assert self.trigger_dim[0][0]% self.img_ndm  == 0 or self.img_ndm % self.trigger_dim[0][0] == 0  # Downsample or upsample dm for plot0
            if self.ntrees > 1:
                assert all([tup[0] % self.img_ndm/(self.ntrees -1)  == 0 or self.img_ndm/(self.ntrees - 1) % tup[0] == 0 for tup in self.trigger_dim[1:]])  # For plot1 +
            assert all(tup[1] % (self.nt_chunk_ds[-1]) == 0 or self.nt_chunk_ds[0] % tup[1] == 0 for tup in self.trigger_dim)  # Downsample or upsample t
            assert self.img_nt % self.nt_chunk_ds[-1] == 0  # Each chunk evenly divides the plots
            assert (self.ntrees - 1) % (self.nplot_groups - 1) == 0  # Trees are evenly divisible between plots

            # Convention for plot groups: 1st half will be the tree0 plot and 2nd half will be the plot for the rest of the trees
            for i in self.nplot_groups*range(self.n_zoom):
                self.add_plot_group("waterfall", nt_per_pix=self.downsample_nt[i], ny=img_ndm)


    def set_stream(self, stream):
        if stream.nfreq != self.nfreq:
            raise RuntimeError("rf_pipelines: number of frequencies in stream (nfreq=%d) does not match bonsai config file '%s' (nfreq=%d)" 
                               % (stream.nfreq, self.config_filename, self.nfreq))


    def start_substream(self, isubstream, t0):
        if self.deallocate_between_substreams:
            self.dedisperser.allocate()

        if self.make_plot:
            self.isubstream = isubstream
            self.plot_groups = [Plotter(self, ntrees=(1 if i==0 else (self.ntrees - 1) / (self.nplot_groups - 1)), iplot=i) for i in xrange(self.nplot_groups)]
            self.json_per_substream["n_plot_groups"] = self.nplot_groups  # Helpful parameter for the web viewer


    def process_chunk(self, t0, t1, intensity, weights, pp_intensity, pp_weights):
        # FIXME some day I'd like to remove this extra copy required by current cython implementation
        intensity = np.array(intensity, dtype=np.float32, order='C')
        weights = np.array(weights, dtype=np.float32, order='C')

        # Send the inputs (intensity, weights) to the dedisperser.
        self.dedisperser.run(intensity, weights, t0)

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
                raise RuntimeError('bonsai returned Inf or NaN triggers! ' + 
                                   'Try reducing the normalization of the intensity and weights arrays, or setting nbits=32 in the bonsai config.')

            # Max along the beta and SM indices, then add the new triggers to the plots! 
            flat_triggers = [np.amax(np.amax(tree, axis=1), axis=1) for tree in triggers]
            i = 0
            while i < self.nplot_groups:
                if i == 0:
                    self.plot_groups[i].process([flat_triggers[0]])
                else:
                    self.plot_groups[i].process(flat_triggers[i : i + (self.ntrees - 1) / (self.nplot_groups - 1)])
                i += (1 if i == 0 else (self.ntrees - 1) / (self.nplot_groups - 1))


    def end_substream(self):
        if self.make_plot:
            for zoom_level in xrange(self.n_zoom):
                for i in xrange(self.nplot_groups):
                    self._write_file(self.plot_groups[i].plots[zoom_level], 
                                     zoom_level,
                                     self.plot_groups[i].downsample_nt[zoom_level],
                                     self.plot_groups[i].ifile[zoom_level],
                                     self.plot_groups[i].iplot,
                                     self.n_zoom)

        self.dedisperser.end_dedispersion()
        
        if self.global_max_tracker is not None:
            # Add global max trigger data to pipeline json output
            self.json_per_substream["frb_global_max_trigger"] = self.global_max_tracker.global_max_trigger
            self.json_per_substream["frb_global_max_trigger_dm"] = self.global_max_tracker.global_max_trigger_dm
            self.json_per_substream["frb_global_max_trigger_tfinal"] = self.global_max_tracker.global_max_trigger_arrival_time

        if self.deallocate_between_substreams:
            self.dedisperser.deallocate()


    def _write_file(self, arr, zoom_level, downsample_nt, ifile, iplot, nzoom):
        # When we reach end-of-stream, the buffer might be partially full. In this case, pad with black.
        basename = self.img_prefix[zoom_level] + '_plot' + str(iplot)
        if self.isubstream > 0:
            basename += str(isubstream+1)
        basename += ('_%s.png' % ifile)

        group_id = iplot * nzoom + zoom_level

        # The add_plot() method adds the plot to the JSON output, and returns the filename that should be written.
        filename = self.add_plot(basename,
                                 it0 = int(ifile * self.img_nt * downsample_nt),
                                 nt = self.img_nt * downsample_nt, 
                                 nx = arr.shape[1],
                                 ny = arr.shape[0], 
                                 group_id = group_id)
        
        if self.dynamic_plotter:
            rf_pipelines.write_png(filename, arr, transpose=True)
        else:
            rf_pipelines.utils.triggers_png(filename, arr, transpose=True, 
                                            threshold1=self.plot_threshold1, threshold2=self.plot_threshold2)


class Plotter():
    """A plotter object holds all desired zoom levels for a plot"""
    def __init__(self, transform, ntrees, iplot):
        self.transform = transform       # To use _write_file()
        self.iplot = iplot               # Helpful for establishing plot group
        self.nzoom = transform.n_zoom    # Number of zoom levels to be produced
        self.nx = transform.img_nt       # Number of x pixels per plot
        self.ny = transform.img_ndm      # Number of y pixels per plot
        self.ntrees = ntrees             # Number of trees being __plotted together__ (different from transform.ntrees)
        self.nt_chunk_ds = transform.nt_chunk_ds   # Number of x pixels that should be written per chunk (2D array)
        self.ndm = self.nx / self.ntrees # Number of y pixels per tree
        self.ix = np.zeros(self.nzoom)   # Keep track of what x position to add chunks to
        self.isubstream = transform.isubstream
        self.img_prefix = transform.img_prefix
        self.downsample_nt = transform.downsample_nt
        self.ifile = np.zeros(self.nzoom)
        self.plots = np.zeros((self.nzoom, self.ny, self.nx))        

        assert self.nzoom == len(self.nt_chunk_ds)


    def process(self, arrs):
        # Check that the arrays passed to process contain the expected number of trees
        assert len(arrs) == self.ntrees
        shape = [triggers.shape for triggers in arrs]

        # Zooming only happens in the time axis, so we can reshape the dm axis outside of the loop
        for i in xrange(self.ntrees):
            if shape[i][0] > self.ndm:
                arrs[i] = self.max_ds(arrs[i], self.ndm, shape[i][1])
            elif shape[i][0] < self.ndm:
                arrs[i] = rf_pipelines.upsample(arrs[i], self.ndm, shape[i][1])

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
                self.plots[zoom_level, iy:iy+self.ndm, self.ix[zoom_level]:self.ix[zoom_level]+self.nt_chunk_ds[zoom_level]] = dm_t
            iy += self.ndm
        self.ix += self.nt_chunk_ds  # We can do this here because of the assertion that each chunks evenly divide into plots

        # After adding, check whether we need to write any files
        for zoom_level in xrange(self.nzoom):
            if self.ix[zoom_level] >= self.nx:
                self.transform._write_file(self.plots[zoom_level, :, :], zoom_level, self.downsample_nt[zoom_level], self.ifile[zoom_level], self.iplot, self.nzoom)
                self.ifile[zoom_level] += 1
                self.ix[zoom_level] = 0
                self.plots[zoom_level, :, :] = 0
                

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
