import sys
import numpy as np

from rf_pipelines.rf_pipelines_c import pipeline_object, wi_transform
from rf_pipelines.utils import write_png, triggers_png, upsample
from rf_pipelines.L1b import L1Grouper


class bonsai_dedisperser(wi_transform):
    """
    Returns a "transform" which doesn't actually modify the data, it just runs the bonsai dedisperser.


    Constructor arguments
    ---------------------

       - config_filename:  The configuration file used to initialize the dedisperser.  

       - fill_rfi_mask: Determines whether the online_mask_filler will be run.

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

       - dyanmic_plotter, plot_threshold1, plot_threshold2: if the first is true, the old red-blue plotter
           will be used, else the new fixed-scale plotter will be used (the second two arguments are arguments
           for the new plotter, representing the sigma value for colour transitions) 

       - event_outfile: None for running without peak finding or a full path for running with
           peak finding (i.e. L1Grouper) and writing the results in a text file

       - (L1Grouper_thr, L1Grouper_beam, L1Grouper_addr) are L1Grouper parameters (see its docstring!)

       - plot_all_trees: if False, only tree 0 will be plotted; if True, all trees will be plotted.

    """

    def __init__(self, config_filename, fill_rfi_mask, img_prefix=None, img_ndm=256, img_nt=256, downsample_nt=1,
                 n_zoom=1, track_global_max=False, dm_min=None, dm_max=None, hdf5_output_filename=None, nt_per_hdf5_file=0,
                 deallocate_between_substreams=False, use_analytic_normalization=False, dynamic_plotter=False,
                 plot_threshold1=6, plot_threshold2=10, event_outfile=None, L1Grouper_thr=7, L1Grouper_beam=0, 
                 L1Grouper_addr=None, plot_all_trees=False):

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
        wi_transform.__init__(self, name)

        # Need to save all constructor arguments, for later use in jsonize()
        self.config_filename = config_filename
        self.fill_rfi_mask = fill_rfi_mask
        self.hdf5_output_filename = hdf5_output_filename
        self.use_analytic_normalization = use_analytic_normalization
        self.deallocate_between_substreams = deallocate_between_substreams
        self._img_prefix = img_prefix  # special case: self.img_prefix has another meaning
        
        initially_allocated = not deallocate_between_substreams
        self.dedisperser = bonsai.Dedisperser(config_filename, fill_rfi_mask=fill_rfi_mask, allocate=initially_allocated, use_analytic_normalization=use_analytic_normalization)
        self.global_max_tracker = None

        if track_global_max:
            self.global_max_tracker_dm_min = dm_min
            self.global_max_tracker_dm_max = dm_max
            self.global_max_tracker = bonsai.global_max_tracker(dm_min, dm_max)
            self.dedisperser.add_processor(self.global_max_tracker)

        if hdf5_output_filename:
            t = bonsai.trigger_hdf5_file_writer(hdf5_output_filename, nt_per_hdf5_file)
            self.dedisperser.add_processor(t)
            self.nt_per_hdf5_file = nt_per_hdf5_file

        # For grouper code
        self.event_outfile = event_outfile
        if self.event_outfile is not None:
            self.L1Grouper_thr = L1Grouper_thr
            self.L1Grouper_beam = L1Grouper_beam
            self.L1Grouper_addr = L1Grouper_addr
            self.grouper = L1Grouper(self.dedisperser, L1Grouper_thr, L1Grouper_beam, L1Grouper_addr)
            self.detected_events = []

        # The 'nt_chunk' parameter is also determined by the config file, not a constructor argument.
        # The same is true for 'nfreq', but rather than initializing self.nfreq here, we let bind()
        # initialize self.nfreq to the "stream nfreq", and check in _bind_transform() that it agrees
        # with the "dedisperser nfreq".  (We do it this way so that the error message will be more helpful.)
        self.nt_chunk = self.dedisperser.nt_chunk

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
            self.plot_all_trees = plot_all_trees
            assert self.ntrees > 0
            assert self.trigger_dim[0][0] % self.img_ndm  == 0 or self.img_ndm % self.trigger_dim[0][0] == 0  # Downsample or upsample dm for plot0
            if self.plot_all_trees:
                assert all([tup[0] % (self.img_ndm / 2)  == 0 
                            or (self.img_ndm / 2) % tup[0] == 0 for tup in self.trigger_dim[1:]])
            assert all(tup[1] % (self.nt_chunk_ds[-1]) == 0 or self.nt_chunk_ds[0] % tup[1] == 0 for tup in self.trigger_dim)  # Downsample or upsample t


    def _bind_transform(self, json_data):
        if self.nfreq != self.dedisperser.nfreq:
            raise RuntimeError("rf_pipelines: number of frequencies in stream (nfreq=%d) does not match bonsai config file '%s' (nfreq=%d)" 
                               % (self.nfreq, self.config_filename, self.dedisperser.nfreq))


    def _allocate(self):
        if self.deallocate_between_substreams:
            self.dedisperser.allocate()

        if self.make_plot:
            self.plot_groups = [Plotter(self, ny=self.img_ndm/(1 if i==0 else 2), iplot=i) for i in xrange(1 + self.plot_all_trees * (self.ntrees-1))]


    def _start_pipeline(self, json_attrs):
        # Calls to add_plot_group() must go here.

        if not self.make_plot:
            return

        for i in xrange(self.n_zoom):
            self.add_plot_group("waterfall", nt_per_pix=self.downsample_nt[i], ny=self.img_ndm)

        if self.plot_all_trees:
            for i in (self.ntrees-1) * range(self.n_zoom):
                self.add_plot_group("waterfall", nt_per_pix=self.downsample_nt[i], ny=self.img_ndm//2)


    def _process_chunk(self, intensity, weights, pos):
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
                raise RuntimeError('bonsai returned Inf or NaN triggers! ' + 
                                   'Try reducing the normalization of the intensity and weights arrays, or setting nbits=32 in the bonsai config.')

            # Max along the beta and SM indices, then add the new triggers to the plots! 
            flat_triggers = [np.amax(np.amax(tree, axis=1), axis=1) for tree in triggers]

            if self.plot_all_trees:
                for i in xrange(1 + self.plot_all_trees * (self.ntrees-1)):
                    self.plot_groups[i].process(flat_triggers[i])
            else:
                self.plot_groups[0].process(flat_triggers[0])
           
        # For grouper code
        if self.event_outfile is not None:
            events = self.grouper.process_triggers()
            if events is not None:
                self.detected_events.append(events) 
 

    def _end_pipeline(self, json_output):
        if self.make_plot:
            # Helpful parameter for the web viewer
            json_output["n_plot_groups"] = 1 + self.plot_all_trees * (self.ntrees-1) 

            # Write whatever may be left in the plots
            for zoom_level in xrange(self.n_zoom):
                for i in xrange(1 + self.plot_all_trees * (self.ntrees-1)):
                    self._write_file(self.plot_groups[i].plots[zoom_level], 
                                     zoom_level,
                                     self.plot_groups[i].ifile[zoom_level],
                                     self.plot_groups[i].iplot)

        self.dedisperser.end_dedispersion()
        
        if self.global_max_tracker is not None:
            # Add global max trigger data to pipeline json output
            json_output["frb_global_max_trigger"] = self.global_max_tracker.global_max_trigger
            json_output["frb_global_max_trigger_dm"] = self.global_max_tracker.global_max_trigger_dm
            json_output["frb_global_max_trigger_tfinal"] = self.global_max_tracker.global_max_trigger_arrival_time

        if self.deallocate_between_substreams:
            self.dedisperser.deallocate()

        # For grouper code
        if self.event_outfile is not None and self.detected_events:
            self.detected_events = np.hstack(self.detected_events)
            np.save(self.event_outfile, self.detected_events)
            with open(self.event_outfile, 'w') as f:
                for i, event in enumerate(self.detected_events):
                    print ("--- Recovered Pulse --- # %d, DM: %.2f,  SNR: %.2f, Arrival Time: %.4f sec"
                           % (i+1, event['dm'], event['snr'], 
                              event['time'].astype(float)/1e6))
                    if i == 0:
                        f.write("# \t DM \t SNR \t Arrival Time (sec)\n")
                    f.write("%d \t %.2f \t %.2f \t %.4f\n" % (i+1, event['dm'], event['snr'], event['time'].astype(float)/1e6))


    def _write_file(self, arr, zoom_level, ifile, iplot):
        basename = self.img_prefix[zoom_level] + '_tree' + str(iplot)
        basename += ('_%s.png' % ifile)

        group_id = iplot * self.n_zoom + zoom_level

        # The add_plot() method adds the plot to the JSON output, and returns the filename that should be written.
        filename = self.add_plot(basename,
                                 it0 = int(ifile * self.img_nt * self.downsample_nt[zoom_level]),
                                 nt = self.img_nt * self.downsample_nt[zoom_level],
                                 nx = arr.shape[1],
                                 ny = arr.shape[0], 
                                 group_id = group_id)
        
        if self.dynamic_plotter:
            write_png(filename, arr, transpose=True)
        else:
            triggers_png(filename, arr, transpose=True, threshold1=self.plot_threshold1, threshold2=self.plot_threshold2)


    def jsonize(self):
        return { 'class_name': 'bonsai_dedisperser',
                 'config_filename': self.config_filename,
                 'fill_rfi_mask': self.fill_rfi_mask,
                 'img_prefix': self._img_prefix if self.make_plot else None,
                 'img_ndm': self.img_ndm if self.make_plot else 0,
                 'img_nt': self.img_nt if self.make_plot else 0,
                 'downsample_nt': self.downsample_nt[0] if self.make_plot else 0,
                 'n_zoom': self.n_zoom if self.make_plot else 0,
                 'track_global_max': (self.global_max_tracker is not None),
                 'global_max_tracker_dm_min': self.global_max_tracker_dm_min if (self.global_max_tracker is not None) else None,
                 'global_max_tracker_dm_max': self.global_max_tracker_dm_min if (self.global_max_tracker is not None) else None,
                 'hdf5_output_filename': self.hdf5_output_filename,
                 'nt_per_hdf5_file': self.nt_per_hdf5_file if self.hdf5_output_filename else 0,
                 'deallocate_between_substreams': self.deallocate_between_substreams,
                 'use_analytic_normalization': self.use_analytic_normalization,
                 'dynamic_plotter': self.dynamic_plotter if self.make_plot else False,
                 'plot_threshold1': self.plot_threshold1 if self.make_plot else 6,
                 'plot_threshold2': self.plot_threshold2 if self.make_plot else 10,
                 'event_outfile': self.event_outfile,
                 'L1Grouper_thr': self.L1Grouper_thr if (self.event_outfile is not None) else 7,
                 'L1Grouper_beam': self.L1Grouper_beam if (self.event_outfile is not None) else 0,
                 'L1Grouper_addr': self.L1Grouper_addr if (self.event_outfile is not None) else None,
                 'plot_all_trees': self.plot_all_trees if self.make_plot else False }


    @staticmethod
    def from_json(j):
        return bonsai_dedisperser(config_filename = j['config_filename'],
                                  fill_rfi_mask = j['fill_rfi_mask'],
                                  img_prefix = j['img_prefix'],
                                  img_ndm = j['img_ndm'],
                                  img_nt = j['img_nt'],
                                  downsample_nt = j['downsample_nt'],
                                  n_zoom = j['n_zoom'],
                                  track_global_max = j['track_global_max'],
                                  dm_min = j['global_max_tracker_dm_min'],
                                  dm_max = j['global_max_tracker_dm_max'],
                                  hdf5_output_filename = j['hdf5_output_filename'],
                                  nt_per_hdf5_file = j['nt_per_hdf5_file'],
                                  deallocate_between_substreams = j['deallocate_between_substreams'],
                                  use_analytic_normalization = j['use_analytic_normalization'],
                                  dynamic_plotter = j['dynamic_plotter'],
                                  plot_threshold1 = j['plot_threshold1'],
                                  plot_threshold2 = j['plot_threshold2'],
                                  event_outfile = j['event_outfile'],
                                  L1Grouper_thr = j['L1Grouper_thr'],
                                  L1Grouper_beam = j['L1Grouper_beam'],
                                  L1Grouper_addr = j['L1Grouper_addr'],
                                  plot_all_trees = j['plot_all_trees'])


class Plotter():
    """A plotter object holds all desired zoom levels for a plot."""

    def __init__(self, transform, ny, iplot):
        self.transform = transform       # Access transform parameters
        self.iplot = iplot               # Helpful for establishing plot group
        self.ny = ny                     # Number of y pixels that will be written
        self.ix = np.zeros(self.transform.n_zoom, dtype=np.int)   # Keep track of what x position to add chunks to
        self.ifile = np.zeros(self.transform.n_zoom, dtype=np.int)
        self.plots = np.zeros((self.transform.n_zoom, self.ny, self.transform.img_nt))

        assert transform.n_zoom == len(self.transform.nt_chunk_ds)


    def process(self, arr):
        # Check that the arrays passed to process contain the expected number of trees
        shape = arr.shape

        # Zooming only happens in the time axis, so we can reshape the dm axis outside of the loop
        if shape[0] > self.ny:
            arr = self.max_ds(arr, self.ny, shape[1])
        elif shape[0] < self.ny:
            arr = upsample(arr, self.ny, shape[1])

        for zoom_level in xrange(self.transform.n_zoom):
            dm_t = arr.copy()
            dm_t_shape = dm_t.shape
            # In the x (time) axis, we need to transform dm_t_shape[1] to self.nt_chunk_ds[zoom_level] - may need to up/downsample
            if dm_t_shape[1] > self.transform.nt_chunk_ds[zoom_level]:
                dm_t = self.max_ds(dm_t, dm_t_shape[0], self.transform.nt_chunk_ds[zoom_level])
            elif dm_t_shape[1] < self.transform.nt_chunk_ds[zoom_level]:
                dm_t = upsample(dm_t, dm_t_shape[0], self.transform.nt_chunk_ds[zoom_level])

            ichunk = 0
            while ichunk < self.transform.nt_chunk_ds[zoom_level]:
                # Move to end of chunk or end of current plot, whichever comes first.
                n = min(self.transform.nt_chunk_ds[zoom_level] - ichunk, self.transform.img_nt - self.ix[zoom_level])
                assert n > 0
                self.plots[zoom_level, :, int(self.ix[zoom_level]):(int(self.ix[zoom_level])+n)] = dm_t[:, ichunk:(ichunk+n)]

                self.ix[zoom_level] += n
                ichunk += n

                # Check whether to write out
                if self.ix[zoom_level] >= self.transform.img_nt:
                    self.transform._write_file(self.plots[zoom_level, :, :], 
                                               zoom_level, 
                                               self.ifile[zoom_level], 
                                               self.iplot)
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


pipeline_object.register_json_deserializer('bonsai_dedisperser', bonsai_dedisperser.from_json)
