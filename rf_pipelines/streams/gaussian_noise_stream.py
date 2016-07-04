"""
Gaussian noise stream, useful in conjunction with frb_injector_transform for simulating FRB data.

This is a thin wrapper around a C++ implmentation.
See gaussian_noise_stream.cpp, and python linkage in rf_pipelines_c.cpp.
"""

from rf_pipelines import rf_pipelines_c


def gaussian_noise_stream(nfreq, nt_tot, freq_lo_MHz, freq_hi_MHz, dt_sample, sample_rms=1.0, nt_chunk=0):
    """
    Returns a weighted intensity stream (wi_stream) which simulates Gaussian random noise for each frequency channel and time sample.
    
    The 'nt_tot' arg is the total length of the stream, in samples.
    The 'dt_sample' arg is the length of a sample in seconds.
    The 'sample_rms' arg is the Gaussian RMS of a single (freq_channel, time_sample) pair.

    The 'nt_chunk' arg is the chunk size used internally when moving data into the rf_pipelines buffer.
    If unspecified or zero, it will default to a reasonable value.
    """

    return rf_pipelines_c.make_gaussian_noise_stream(nfreq, nt_tot, freq_lo_MHz, freq_hi_MHz, dt_sample, sample_rms, nt_chunk)
