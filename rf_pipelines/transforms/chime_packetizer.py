"""
CHIME packetizer.  This is a thin wrapper around a C++ implementation (rf_pipelines/chime_packetizer.cpp)
which in turn in a wrapper around code in libch_frb_io.
"""

from rf_pipelines import rf_pipelines_c


def chime_packetizer(dstname, nfreq_per_packet, nt_per_chunk, nt_per_packet, wt_cutoff, target_gbps):
    """
    Returns a pseudo-transform which packetizes the stream without modifying it.  The chime_packetizer
    can be inserted anywhere in the pipeline, and will packetize the data as processed up to that point.

    The 'dstname' argument is a string of the form HOSTNAME:PORT.  For example 'localhost:13178' or
    'chimer.physics.ubc.ca:13178'.  (Be careful sending packets over the internet since the bandwidth
    can be very high!)

    The 'wt_cutoff' argument is used to convert the rf_pipelines 'weights' array to a boolean mask.
    This conversion is necessary because the CHIME L0_L1 packet format doesn't support a floating-point
    weight array.  Samples with weight below the cutoff will be masked.

    If the 'target_gbps' argument is nonzero, then output will be "throttled" to the target bandwidth, specified
    in Gbps.  If target_gbps=0, then packets will be sent as quickly as possible.
    """

    return rf_pipelines_c.make_chime_packetizer(dstname, nfreq_per_packet, nt_per_chunk, nt_per_packet, wt_cutoff, target_gbps)

