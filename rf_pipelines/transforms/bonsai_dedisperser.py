import sys
import numpy as np
import rf_pipelines


class bonsai_dedisperser(rf_pipelines.py_wi_transform):
    """
    Returns a "transform" which doesn't actually modify the data, it just runs the bonsai dedisperser.


    Constructor arguments
    ---------------------

       - config_hdf5_filename:  The configuration file used to initialize the dedisperser.  
           This must be produced with the program 'bonsai-mkweight' in the bonsai github repo.

       - img_prefix: Determines output filenames, using a similar convention to the plotter_transform.
           There should be an option to create a bonsai_dedisperser without trigger plots, so let's
           say that if img_prefix=None, then no plots are generated.

       - img_ndm: Number of y-pixels in each output plot (same as img_nfreq in the plotter_transform).
           Note that the y-axis of the bonsai plots corresponds to dispersion measure.

       - img_nt: Number of x-pixels in each output plot.

       - downsample_nt: If > 1, then each x-pixel in the plot will correspond to multiple time samples.

       - n_zoom: Number of zoom levels in plots.

       - trigger_hdf5_filename: If specified, coarse-grained triggers will be written to an HDF5 file.


    FIXME: Currently the dedisperser must be initialized from a config hdf5 file (rather than
    the simpler config text file) since we use analytic weights to normalize the triggers.
    Since the analytic weights are only correct for unit-variance noise, the trigger normalization
    will be wrong for a real experiment, and the triggers won't be meaningfully normalized to
    "sigmas".  All of this is just a placeholder until Monte Carlo trigger variance estimation
    is implemented in bonsai.
    """

    def __init__(self, config_hdf5_filename, img_prefix=None, img_ndm=None, img_nt=None, downsample_nt=1, n_zoom=1, trigger_hdf5_filename=None):
        # We import the bonsai module here, rather than at the top of the file, so that bonsai isn't
        # required to import rf_pipelines (but is required when you try to construct a bonsai_dedisperser).
        try:
            import bonsai
        except ImportError:
            raise RuntimeError("rf_pipelines: couldn't import the 'bonsai' module.  You may need to clone https://github.com/CHIMEFRB/bonsai and install.")

        name = "bonsai_dedisperser('%s')" % config_hdf5_filename
        rf_pipelines.py_wi_transform.__init__(self, name)

        config_params = bonsai.ConfigParams(config_hdf5_filename, True)

        self.config_hdf5_filename = config_hdf5_filename
        self.dedisperser = bonsai.Dedisperser(config_params, True)
        
        # Note that 'nfreq' is determined by the config file.  If the stream's 'nfreq' differs,
        # then an exception will be thrown.  The 'nt_chunk' parameter is also determined by the
        # config file, not a constructor argument.

        self.nfreq = self.dedisperser.nfreq
        self.nt_chunk = self.dedisperser.nt_data
        self.nt_prepad = 0
        self.nt_postpad = 0

        if trigger_hdf5_filename:
            self.dedisperser.start_trigger_file(trigger_hdf5_filename, nt_per_file=0)

        # Activates some fields in self.dedisperser which are used in the frb_olympics.
        self.dedisperser.global_max_trigger_active = True


    def set_stream(self, stream):
        if stream.nfreq != self.nfreq:
            raise RuntimeError("rf_pipelines: number of frequencies in stream (nfreq=%d) does not match bonsai config file '%s' (nfreq=%d)" % (stream.nfreq, self.config_hdf5_filename, self.nfreq))


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

        triggers = self.dedisperser.get_triggers()

        # Plotter logic can go here!


    def end_substream(self):
        self.dedisperser.end_dedispersion()


####################################################################################################


def old_bonsai_dedisperser(config_hdf5_filename, trigger_hdf5_filename=None, trigger_plot_stem=None, nt_per_file=0, ibeam=0):
    """
    This is the old C++ bonsai_dedisperser, which we're trying to phase out, in favor of the
    python implementation which has been partially implemented above!

    The plotting behavior of the C++ dedisperser is not what we want:

       - one plot_group per tree
       - number of y-pixels in plots is determined by bonsai config file (not selectable)
       - time downsampling factor in plots is determined by bonasi config file (not selectable)
    
    When the new plotting behavior is implemented in the python bonsai_dedisperser, then the
    old_bonsai_dedisperser can be removed.
    """

    if trigger_hdf5_filename is None:
        trigger_hdf5_filename = ''

    if trigger_plot_stem is None:
        trigger_plot_stem = ''

    # Note: 'ibeam' argument ignored, as of bonsai v7_devel.
    return rf_pipelines_c.make_bonsai_dedisperser(config_hdf5_filename, trigger_hdf5_filename, trigger_plot_stem, nt_per_file)
