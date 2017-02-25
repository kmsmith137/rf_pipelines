import sys
import numpy as np
import rf_pipelines
import rf_pipelines.rf_pipelines_c as rf_pipelines_c

class bonsai_dedisperser(rf_pipelines.py_wi_transform):
    """
    Returns a "transform" which doesn't actually modify the data, it just runs the bonsai dedisperser.


    Constructor arguments
    ---------------------

       - config_filename:  The configuration file used to initialize the dedisperser.  

       - img_prefix: Determines output filenames, using a similar convention to the plotter_transform.
           There should be an option to create a bonsai_dedisperser without trigger plots, so let's
           say that if img_prefix=None, then no plots are generated.

       - img_ndm: Number of y-pixels in each output plot (same as img_nfreq in the plotter_transform).
           Note that the y-axis of the bonsai plots corresponds to dispersion measure.

       - img_nt: Number of x-pixels in each output plot.

       - downsample_nt: If > 1, then each x-pixel in the plot will correspond to multiple time samples.

       - n_zoom: Number of zoom levels in plots.

       - track_global_max: If True, then the following json output will be written:
             frb_global_max_trigger, frb_global_max_trigger_dm, frb_global_max_trigger_tfinal

       - deallocate_between_substreams: if True, then the dedisperser will be deallocated initialially,
           allocate buffers in start_substream(), and deallocate in end_substream().

    Currently, we use "analytic weights", i.e. the trigger normalization for a toy noise model in
    which every (frequency_channel, time) sample is a *unit* Gaussian.  This is a placeholder for
    implementing something better soon!

    If 'config_filename' is an hdf5 file with precomputed analytic weights, they will be read from
    disk.  Otherwise they will be computed in the transform constructor.
    
    If the deallocate_between_substreams=False, then the dedisperser will be deallocated initially,
    allocate buffers in start_substream(), and deallocate in end_substream().
    """

    def __init__(self, config_filename, img_prefix="triggers", img_ndm=256, img_nt=256, downsample_nt=1, n_zoom=1, track_global_max=False, deallocate_between_substreams=False):
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
        self.dedisperser = bonsai.Dedisperser(config_filename, init_analytic_variance=True, allocate=initially_allocated)
        self.global_max_tracker = None

        if track_global_max:
            self.global_max_tracker = bonsai.global_max_tracker()
            self.dedisperser.add_processor(self.global_max_tracker)
        
        # Note that 'nfreq' is determined by the config file.  If the stream's 'nfreq' differs,
        # then an exception will be thrown.  The 'nt_chunk' parameter is also determined by the
        # config file, not a constructor argument.
        self.nfreq = self.dedisperser.nfreq
        self.nt_chunk = self.dedisperser.nt_chunk
        self.nt_prepad = 0
        self.nt_postpad = 0

        # Set plotting parameters
        if self.make_plot:
            self.n_zoom = n_zoom
            self.img_ndm = img_ndm
            self.img_nt = img_nt
            self.downsample_nt = [downsample_nt]
            self.nt_chunk_ds = [self.nt_chunk // self.downsample_nt[0]]
            self.img_prefix = [str(img_prefix) + "_zoom0"]
            self.add_plot_group("waterfall", nt_per_pix=downsample_nt, ny=img_ndm)
            if self.n_zoom > 1:
                for zoom_level in xrange(self.n_zoom - 1):
                    self.downsample_nt += [self.downsample_nt[zoom_level] * 2]   # zoom_level = previous element's index because of the original value added
                    self.nt_chunk_ds += [self.nt_chunk // self.downsample_nt[zoom_level + 1]]
                    self.add_plot_group("waterfall", nt_per_pix=self.downsample_nt[zoom_level + 1], ny=img_ndm)
                    self.img_prefix += [img_prefix + "_zoom" + str(zoom_level+1)] 
            self.dimensions_init = False  # temporary! 

            if self.nt_chunk % self.downsample_nt[-1] != 0:
                raise RuntimeError("bonsai plotter transform: specified nt_chunk(=%d) must be a multiple of downsampling factor at max zoom level (=%d)" 
                                   % (self.nt_chunk, self.downsample_nt[-1]))
        

    def set_stream(self, stream):
        if stream.nfreq != self.nfreq:
            raise RuntimeError("rf_pipelines: number of frequencies in stream (nfreq=%d) does not match bonsai config file '%s' (nfreq=%d)" % (stream.nfreq, self.config_filename, self.nfreq))


    def start_substream(self, isubstream, t0):
        if self.deallocate_between_substreams:
            self.dedisperser.allocate()

        if self.make_plot:
            self.buf = np.zeros((self.n_zoom, self.img_ndm, self.img_nt), dtype=np.float32)
            self.isubstream = isubstream
            self.ifile = np.zeros((self.n_zoom))    # keeps track of which png file we're accumulating 
            self.ipos = np.zeros((self.n_zoom))     # keeps track of how many (downsampled) time samples have been accumulated into file so far


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

            if not all(np.all(np.isfinite(t)) for t in triggers):
                # Was the problem in the input arrays... ?
                if not np.all(np.isfinite(intensity)) or not np.all(np.isfinite(weights)):
                    raise RuntimeError('bonsai_dedisperser: input intensity/weights arrays contained Inf or NaN!')

                # If not, then maybe it's a 16-bit overflow... ?
                raise RuntimeError('bonsai returned Inf or NaN triggers!  Try reducing the normalization of the intensity and weights arrays, or setting nbits=32 in the bonsai config.')

            # First, let's flatten the SM_index and beta_index axes by taking max values to get an array indexed only by dm and time
            preserved_dm_t = np.amax(np.amax(triggers[0], axis=1), axis=1)
    
            # Here we check that some parameters are okay - temp. until trigger dimensions can be accessed from __init__! 
            if self.dimensions_init == False:
                self.trigger_dim = preserved_dm_t.shape
                assert self.trigger_dim[0] % self.img_ndm  == 0 or self.img_ndm % self.trigger_dim[0] == 0   # downsample or upsample dm
                assert self.trigger_dim[1] % (self.nt_chunk_ds[-1]) == 0 or self.nt_chunk_ds[0] % self.trigger_dim[1] == 0   # downsample or upsample t
                self.dimensions_init = True

            # Because "zooming" only happens in the time axis, we can reshape the dm axis outside of the loop
            # In the y (dm) axis, we need to transform self.trigger_dim[0] to self.img_ndm - may need to downsample or upsample
            if self.trigger_dim[0] > self.img_ndm:
                preserved_dm_t = self._max_downsample(preserved_dm_t, self.img_ndm, preserved_dm_t.shape[1])
            elif self.trigger_dim[0] < self.img_ndm:
                preserved_dm_t = rf_pipelines.upsample(preserved_dm_t, self.img_ndm, preserved_dm_t.shape[1])
                    
            for zoom_level in xrange(self.n_zoom): 
                dm_t = preserved_dm_t.copy()

                # In the x (time) axis, we need to transform self.trigger_dim[1] to self.nt_chunk / self.downsample_nt - may need to downsample or upsample
                if dm_t.shape[1] > self.nt_chunk_ds[zoom_level]:
                    dm_t = self._max_downsample(dm_t, dm_t.shape[0], self.nt_chunk_ds[zoom_level])
                elif dm_t.shape[1] < self.nt_chunk_ds[zoom_level]:
                    dm_t = rf_pipelines.upsample(dm_t, dm_t.shape[0], self.nt_chunk_ds[zoom_level])
 
                # Now the array will be scaled properly to stick into the plot accumulator array
                ichunk = 0
                while ichunk < self.nt_chunk_ds[zoom_level]:
                    # Move to end of chunk or end of current plot, whichever comes first.                                                                                        
                    n = min(self.nt_chunk_ds[zoom_level] - ichunk, self.img_nt - self.ipos[zoom_level])
                    assert n > 0
                    self.buf[zoom_level, :, self.ipos[zoom_level]:(self.ipos[zoom_level]+n)] = dm_t[:, ichunk:(ichunk+n)]
                    self.ipos[zoom_level] += n
                    ichunk += n
            
                    if self.ipos[zoom_level] >= self.img_nt:
                        self._write_file(zoom_level)
 

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


    def _max_downsample(self, arr, new_dm, new_t):
        """Takes maxima along axes"""
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


    def _write_file(self, zoom_level):
        # When we reach end-of-stream, the buffer might be partially full (i.e. self.ipos < self.img_nt).                                                                                           
        # In this case, pad with black                                                                                                  
        basename = self.img_prefix[zoom_level]
        if self.isubstream > 0:
            basename += str(isubstream+1)
        basename += ('_%s.png' % self.ifile[zoom_level])

        # The add_plot() method adds the plot to the JSON output, and returns the filename that should be written.                                                                                         
        filename = self.add_plot(basename,
                                 it0 = int(self.ifile[zoom_level] * self.img_nt * self.downsample_nt[zoom_level]),
                                 nt = self.img_nt * self.downsample_nt[zoom_level],
                                 nx = self.buf[zoom_level, :, :].shape[1],
                                 ny = self.buf[zoom_level, :, :].shape[0], 
                                 group_id = zoom_level)

        rf_pipelines.write_png(filename, self.buf[zoom_level, :, :], transpose=True)

        self.buf[zoom_level, :, :] = 0.
        self.ifile[zoom_level] += 1
        self.ipos[zoom_level] = 0
