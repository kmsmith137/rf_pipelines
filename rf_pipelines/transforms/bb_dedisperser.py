import rf_pipelines.rf_pipelines_c

# This "wrapper around a wrapper" is kind of silly, but has the convenient side effect
# that it can be called with a flexible python syntax (i.e. arguments can be reordered
# or assigned defaults)

def bb_dedisperser(dm_start, dm_end, dm_tol, pulse_width_ms):
    """
    bb_dedisperser(dm_start, dm_end, dm_tol, pulse_width_ms)
    Returns a wrapper around the Ben Barsdell 'dedisp' GPU code.
    """

    return rf_pipelines.rf_pipelines_c.make_bb_dedisperser(dm_start, dm_end, dm_tol, pulse_width_ms)
