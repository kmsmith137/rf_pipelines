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

        # Set plotting parameters ---- these things all need to change lol
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
            self.add_plot_group("waterfall", nt_per_pix=downsample_nt, ny=img_ndm)
            if self.n_zoom > 1:
                for zoom_level in xrange(self.n_zoom - 1):
                    self.downsample_nt += [self.downsample_nt[zoom_level] * 2]   # zoom_level = previous element's index because of the original value added
                    self.nt_chunk_ds += [self.nt_chunk // self.downsample_nt[zoom_level + 1]]
                    self.add_plot_group("waterfall", nt_per_pix=self.downsample_nt[zoom_level + 1], ny=img_ndm)
                    self.img_prefix += [img_prefix + "_zoom" + str(zoom_level+1)] 

            if self.nt_chunk % self.downsample_nt[-1] != 0:
                raise RuntimeError("bonsai plotter transform: specified nt_chunk(=%d) must be a multiple of downsampling factor at max zoom level (=%d)" 
                                   % (self.nt_chunk, self.downsample_nt[-1]))

            # Set incoming triger dimension paramaters ---- as do these things
            self.trigger_dim = self.dedisperser.ndm_coarse[0], self.dedisperser.nt_coarse_per_chunk[0]
            assert self.trigger_dim[0] % self.img_ndm  == 0 or self.img_ndm % self.trigger_dim[0] == 0   # downsample or upsample dm
            assert self.trigger_dim[1] % (self.nt_chunk_ds[-1]) == 0 or self.nt_chunk_ds[0] % self.trigger_dim[1] == 0   # downsample or upsample t      
            # Add as cheat-y assert, but I think this should be fine. Saves a lot of work
            assert self.img_nt % self.nt_chunk_ds[-1] == 0
            

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



# First, we establish the number of zoom levels
# Then, for each tree, we identify the amount of downsampling we need to do for each zoom level
# Because we only zoom in the time axis, we only need a 1D list of downsampling amounts for the dm axis
# i.e. [dm_ds_tree0, dm_ds_tree1, ..., dm_ds_treen]

# For the time axis, we need a downsampling parameter for each zoom level for each tree, so this calls for
# a 2D numpy array as follows of size (ntrees, nzoom)
# [[ds_t_z0_t0, ds_t_z1_t0, ...,  ds_t_zn_t0], [[ds_t_z0_t1, ds_t_z1_t1, ...,  ds_t_zn_t1], ..., [[ds_t_z0_tn, ds_t_z1_tn, ...,  ds_t_zn_tn]]

# Then, in process_chunk, we need to get the triggers
# Then, for each tree, we should down/upsample in the dm axis as required. Here it is fine to modify the actual trigger arrays

# Next comes the trickier part... doing the time downsampling and constructing the plot objects
# I think the easiest thing to do will be to construct a plot class such that the class contains a plot object
# that actually contains all the zoom levels for that plot... so actually I guess it would be a tree object, not a plot object

# Then, when we pass a dm-down/upsampled array to this tree object with an add() function, it will down/upsample in the time
# access as required by each of the zoom levels and build up the array itself
# This means that each tree object needs to know by how much to down/upsample each chunk in time. This can be accomplished by 
# the nt_chunk_ds parameter, which defines how many x pixels we have per chunk

# In an ideal world, this would be everything to worry about, but no! We must also worry about stacking in the dm axis if
# multiple trees are being plotted in the same plot
# For this purpose, I propose keeping two parameters it and idm, similar to the old ipos counter. 
# The tree object will take as an argument the number. Wait I just realized it's not a tree. It is in the trivial case of 
# just plotting one tree... uhm I guess it is a plot object, but with all the zooms built in. Back to calling it a plot 
# object I guess. Anyway, it will take the number of trees being passed to it as an argument. The add() method will then
# expect a length ntrees list with the arrays for all the trees in it and use a for loop to do the things. Okay
# I think this is good enough to start coding. 

class Plotter():
    """A plotter object holds all desired zoom levels for a plot"""
    def __init__(self, dim, nzoom, ntrees, nt_chunk_ds, img_prefix, downsample_nt, substream=0):
        # dim - (ny, nx) number of y and x pixels per plot
        # nzoom - number of zoom levels
        # ntrees - number of trees that will be used to construct the plot
        # nt_chunk_ds - 2D array containing the number of x pixels per chunk for each (ntrees, nzoom)
        # img_prefix, downsample_nt, substream - all copied from the transform class for write_file
        self.plots = np.zeros((nzoom, ny, nx))
        self.nzoom = nzoom
        self.ntrees = ntrees
        self.nt_chunk_ds = nt_chunk_ds  # number of x pixels per chunk (nzoom-dim list)
        self.ix = np.zeros(nzoom)  # keep track of what x position to add chunks to

        # Blarg stuff for write_file...
        self.img_prefix = img_prefix  # for write_file
        self.ifile = np.zeros(nzoom)  # for write_file
        self.isubstream = substream  # for write_file
        self.downsample_nt = downsample_nt  # for write_file

        assert len(self.nt_chunk_ds) == self.ntrees
        assert all([len(zoom) == self.nzoom for zoom in self.nt_chunk_ds])
        assert np.all(np.sum(self.nt_chunk_ds, axis=1) == ny)


        def add(self, arrs):
            # First, we check that the list of arrays that have been passed to this method is equal to the number of 
            # trees being plotted
            assert len(arrs) == self.ntrees
            # Also note that one thing we are assuming here is that the arr will fit evenly into each plot and that we will
            # never need to buffer 
            # Now, we need to iterate over each tree
            iy = 0 # keep track of where to add each tree vertically
            for i, tree in enumerate(arrs):
                # Now, we need to resize whatever darn number of times arrrrrgh this is going to be sad
                # Iterate over each zoom level I guess :'(
                for zoom_level in range(self.nzoom):
                    # Oh right! I already wrote this! Hooray for copying and pasting old code! 
                    dm_t = tree.copy()
                    dm_t_shape = dm_t.shape
                    # In the x (time) axis, we need to transform self.trigger_dim[1] to self.nt_chunk / self.downsample_nt - may need to downsample or upsample
                    if dm_t_shape[1] > self.nt_chunk_ds[zoom_level]:
                        dm_t = self.max_ds(dm_t, dm_t.shape[0], self.nt_chunk_ds[zoom_level])
                    elif dm_t_shape[1] < self.nt_chunk_ds[zoom_level]:
                        dm_t = rf_pipelines.upsample(dm_t, dm_t.shape[0], self.nt_chunk_ds[zoom_level])
                    # Now the array will be scaled properly to stick into the plot accumulator arrays
                    # Because of the nice cheat-y assert, we don't need to worry about like anything! :D
                    self.plot[zoom_level, self.iy, self.ix[zoom_level]] = dm_t
                    self.ix[zoom_level] += dm_t_shape[1]
                
                iy += dm_t_shape[1]
            for zoom_level in range(self.nzoom):
                if self.ix[zoom_level] >= self.img_nt:
                    self._write_file(zoom_level)
 

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


    def write_file(self, zoom_level):
        # This method is bugging me because it forces us to copy lots of stuff from the transfrom class
        # into the plotter class, but I don't know of a better way other than making write_file a method
        # in rf_pipelines/utils.py, which could be used by both plotter_transform.py and 
        # bonsai_dedisperser.py.

        # When we reach end-of-stream, the buffer might be partially full (i.e. self.ipos < self.img_nt).                                                                                           
        # In this case, pad with black     

        # Also, I need to eventually establish a convention for this group_id thing
        basename = self.img_prefix[zoom_level]
        if self.isubstream > 0:
            basename += str(isubstream+1)
        basename += ('_%s.png' % self.ifile[zoom_level])

        # The add_plot() method adds the plot to the JSON output, and returns the filename that should be written.                                                                                         
        filename = self.add_plot(basename,
                                 it0 = int(self.ifile[zoom_level] * self.nx * self.downsample_nt[zoom_level]),
                                 nt = self.nx * self.downsample_nt[zoom_level],
                                 nx = self.buf[zoom_level, :, :].shape[1],
                                 ny = self.buf[zoom_level, :, :].shape[0], 
                                 group_id = zoom_level)
        
        if self.dynamic_plotter:
            rf_pipelines.write_png(filename, self.buf[zoom_level, :, :], transpose=True)
        else:
            rf_pipelines.utils.triggers_png(filename, self.buf[zoom_level, :, :], transpose=True, 
                                            threshold1=self.plot_threshold1, threshold2=self.plot_threshold2)

        self.buf[zoom_level, :, :] = 0.
        self.ifile[zoom_level] += 1
        self.ipos[zoom_level] = 0
