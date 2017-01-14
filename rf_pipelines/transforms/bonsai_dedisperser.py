"""
Bonsai dedisperser.  This is a thin wrapper around a C++ implementation.
See bonsai_dedisperser.cpp, and python linkage in rf_pipelines_c.cpp.
"""

from rf_pipelines import rf_pipelines_c

def bonsai_dedisperser(config_hdf5_filename, trigger_hdf5_filename=None, trigger_plot_stem=None, nt_per_file=0, ibeam=0):
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

    The 'ibeam' argument determines the assignment of threads to cores and can probably
    be zero except in special situations.

    FIXME: Currently the dedisperser must be initialized from a config hdf5 file (rather than
    the simpler config text file) since we use analytic weights to normalize the triggers.
    Since the analytic weights are only correct for unit-variance noise, the trigger normalization
    will be wrong for a real experiment, and the triggers won't be meaningfully normalized to
    "sigmas".  All of this is just a placeholder until Monte Carlo trigger variance estimation
    is implemented in bonsai.
    """

    if trigger_hdf5_filename is None:
        trigger_hdf5_filename = ''

    if trigger_plot_stem is None:
        trigger_plot_stem = ''

    return rf_pipelines_c.make_bonsai_dedisperser(config_hdf5_filename, trigger_hdf5_filename, trigger_plot_stem, nt_per_file, ibeam)
