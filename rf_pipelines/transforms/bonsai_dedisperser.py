"""
Bonsai dedisperser.  This is a thin wrapper around a C++ implementation.
See simple_detrender.cpp, and python linkage in rf_pipelines_c.cpp.
"""

from rf_pipelines import rf_pipelines_c


def bonsai_dedisperser(config_hdf5_filename, output_hdf5_filename, nt_per_file=0, ibeam=0):
    """
    Returns a "transform" (object of class wi_transform) which doesn't actually modify the data,
    it just runs the bonsai dedisperser.  The output is a stream of coarse-grained triggers
    which are written to an output hdf5 file.  The dedisperser must be initialized from a
    config hdf5 file produced with the program 'bonsai-mkweight' in the bonsai github repo.

    Note that the program 'bonsai-plot-triggers.py' in the bonsai github repo may be useful
    for quick visual inspection of the bonsai output.

    If 'nt_per_file' is zero, then all triggers will be written to a single "monster file".  
    Otherwise multiple files will be written.  Note that nt_per_file is the number of input 
    time samples (the number of coarse-grained triggers is usually much smaller).

    The 'ibeam' argument determines the assignment of threads to cores and can probably
    be zero except in special situations.

    FIXME 1: Currently the only trigger "processing" which can be done is writing the triggers
    to an hdf5 file for later analysis.  It would be better if we could use bonsai's python
    interface, which allows the dedisperser process_triggers() callback to be written in
    python.  Right now, it's not possible to use bonsai's python interface with rf_pipelines!

    FIXME 2: Currently the dedisperser must be initialized from a config hdf5 file (rather than
    the simpler config text file) since we use analytic weights to normalize the triggers.
    Since the analytic weights are only correct for unit-variance noise, the trigger normalization
    will be wrong for a real experiment, and the triggers won't be meaningfully normalized to
    "sigmas".  All of this is just a placeholder until Monte Carlo trigger variance estimation
    is implemented in bonsai.
    """
    
    return rf_pipelines_c.make_bonsai_dedisperser(config_hdf5_filename, output_hdf5_filename, nt_per_file, ibeam)
