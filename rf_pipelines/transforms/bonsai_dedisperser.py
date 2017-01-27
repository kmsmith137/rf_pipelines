import sys
import numpy as np
import rf_pipelines


class bonsai_dedisperser(rf_pipelines.py_wi_transform):
    """
    Returns a "transform" which doesn't actually modify the data, it just runs the bonsai dedisperser.  
    The dedisperser must be initialized from a config hdf5 file produced with the program 
    'bonsai-mkweight' in the bonsai github repo.

    If 'trigger_hdf5_filename' is a nonempty string, then triggers will be written to one
    or more HDF5 output files.  If 'nt_per_file' is zero, then all triggers will be written
    to a single "monster file".  Otherwise multiple files will be written.  Note that nt_per_file
    is the number of input time samples (the number of coarse-grained triggers is usually
    much smaller).

    If 'trigger_plot_stem' is a nonempty string, then realtime trigger plots will be written.  
    In this case, the nt_per_file arg must be positive.  Filenames are of the form
       ${trigger_plot_stem}_${plot_number}_tree${tree_index}.png.

    FIXME: Currently the dedisperser must be initialized from a config hdf5 file (rather than
    the simpler config text file) since we use analytic weights to normalize the triggers.
    Since the analytic weights are only correct for unit-variance noise, the trigger normalization
    will be wrong for a real experiment, and the triggers won't be meaningfully normalized to
    "sigmas".  All of this is just a placeholder until Monte Carlo trigger variance estimation
    is implemented in bonsai.
    """

    def __init__(self, config_hdf5_filename, trigger_hdf5_filename=None, trigger_plot_stem=None, nt_per_file=0):
        try:
            import bonsai
        except ImportError:
            raise RuntimeError("rf_pipelines: couldn't import the 'bonsai' module.  You may need to clone https://github.com/CHIMEFRB/bonsai and install.")

        name = "bonsai_dedisperser('%s')" % config_hdf5_filename
        rf_pipelines.py_wi_transform.__init__(self, name)

        if (trigger_hdf5_filename is not None) or (trigger_plot_stem is not None) or (nt_per_file > 0):
            print >>sys.stderr, 'XXX bonsai_dedisperser warning: trigger_hdf5_filename, trigger_plot_stem, nt_per_file currently ignored'

        config_params = bonsai.ConfigParams(config_hdf5_filename, True)

        self.config_hdf5_filename = config_hdf5_filename
        self.dedisperser = bonsai.Dedisperser(config_params, True)
        self.dedisperser.global_max_trigger_active = True

        self.nfreq = self.dedisperser.nfreq
        self.nt_chunk = self.dedisperser.nt_data
        self.nt_prepad = 0
        self.nt_postpad = 0


    def set_stream(self, stream):
        if stream.nfreq != self.nfreq:
            raise RuntimeError("rf_pipelines: number of frequencies in stream (nfreq=%d) does not match bonsai config file '%s' (nfreq=%d)" % (stream.nfreq, self.config_hdf5_filename, self.nfreq))


    def process_chunk(self, t0, t1, intensity, weights, pp_intensity, pp_weights):
        # FIXME remove this extra copy required by current cython implementation
        intensity = np.array(intensity, dtype=np.float32, order='C')
        weights = np.array(intensity, dtype=np.float32, order='C')

        self.dedisperser.run(intensity, weights)
        self.process_triggers(t0, t1, self.dedisperser.get_triggers())


    # Subclass may wish to override this.
    def process_triggers(self, t0, t1, triggers):
        pass

        
    # Subclass may wish to override this.
    def end_substream(self):
        pass
