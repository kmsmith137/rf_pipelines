import sys
import numpy as np

from rf_pipelines.rf_pipelines_c import pipeline_object, wi_transform


####################################################################################################


import_successful = False

try:
    import bz_fdmt
    import_successful = True
except ImportError:
    pass


####################################################################################################


def is_power_of_two(n):
    return (n > 0) and (n & (n-1)) == 0

def round_up_to_power_of_two(n):
    assert n > 0
    return 2**(int(np.log2(n-0.5)) + 1)


class bz_fdmt_dedisperser(wi_transform):
    def __init__(self, dm_max, nt_in, verbose=True):
        wi_transform.__init__(self, 'bz_fdmt_dedisperser')

        if not import_successful:
            # Rather than throw an exception, we let 'import bz_fdmt' throw an uncaught
            # exception, so that the caller can see what the problem is.
            import bz_fdmt
            raise RuntimeError("rf_pipelines.bz_fdmt_dedisperser internal error: 'import bz_fdmt' worked on the second try?!")

        if dm_max <= 0.0:
            raise RuntimeError("rf_pipelines.bz_fdmt_dedisperser constructor: expected dm_max > 0")

        if not is_power_of_two(nt_in):
            # Can be relaxed.
            raise RuntimeError("rf_pipelines.bz_fdmt_dedisperser: we currently require nt_in to be a power of two")
        
        self.dm_max = dm_max
        self.nt_in = nt_in
        self.verbose = verbose

        
    def _bind_transform(self, json_data):
        if not is_power_of_two(self.nfreq):
            raise RuntimeError("rf_pipelines.bz_fdmt_dedisperser: FDMT expects nfreq to be a power of two")

        for key in [ 'freq_lo_MHz', 'freq_hi_MHz', 'dt_sample' ]:
            if not json_data.has_key(key):
                raise RuntimeError("rf_pipelines.bz_fdmt_dedisperser: expected json_attrs to contain member '%s'" % key)
        
        self.freq_lo_MHz = float(json_data['freq_lo_MHz'])
        self.freq_hi_MHz = float(json_data['freq_hi_MHz'])
        self.dt_sample = float(json_data['dt_sample'])  # seconds

        # Delay at DM=1, across full band and in samples (in floating-point, not rounded to nearest integer sample)
        self.dm0 = 4.15e3 * (self.freq_lo_MHz**(-2) - self.freq_hi_MHz**(-2)) / self.dt_sample

        # FDMT uses one trial DM per time sample
        self.ndm = int(self.dm0 * self.dm_max) + 2

        if self.verbose:
            print 'bz_fdmt_dedisperser: dm_max=%s, ndm=%d' % (self.dm_max, self.ndm)

        if self.nt_in < self.ndm:
            raise RuntimeError("rf_pipelines.bz_fdmt_dedisperser: timestream length 'nt_in' is too small")


    def _allocate(self):
        self.in_buf = np.zeros((self.nfreq, self.nt_in), dtype=np.float32)
    

    def _start_pipeline(self, json_attrs):
        self.runflag = False
        self.max_trigger = 0.0
        self.max_trigger_dm = 0.0
        self.max_trigger_tfinal = 0.0


    def _process_chunk(self, intensity, weights, pos):
        assert intensity.shape == weights.shape == (self.nfreq, self.nt_chunk)

        nt_valid = self.nt_in - pos
        nt_valid = min(nt_valid, self.nt_chunk)
        nt_valid = max(nt_valid, 0)

        # Check weights
        assert(np.all(weights[:,:nt_valid] == 1.0))
        assert(np.all(weights[:,nt_valid:] == 0.0))

        # Note frequency channel ordering is reversed here!
        # rf_pipelines ordering is highest-to-lowest, but FDMT ordering is lowest-to-highest.
        self.in_buf[:,pos:(pos+nt_valid)] = intensity[::-1,:nt_valid]

        if self.runflag or (pos + self.nt_chunk < self.nt_in):
            return

        # If we get here, the dedisperser will be run in this iteration.
        # First, a little trick to get the normalization, by running FDMT on an all-ones array.

        nt0 = round_up_to_power_of_two(self.ndm + 50)
        a = np.ones((self.nfreq, nt0), dtype=np.float32)
        a = bz_fdmt.FDMT(a, self.freq_lo_MHz, self.freq_hi_MHz, self.ndm, np.float32, Verbose=False)
        var = a[:,-1]   # 1D array of length ndm

        # Sanity checks, just to sleep better at night.

        assert a.shape == (self.ndm, nt0)
        assert np.all(var > 0.0)
        assert np.all(np.abs(var - np.array(var+0.5, dtype=np.int)) < 1.0e-5)   # variances should all be integers

        for idm in xrange(self.ndm):
            assert np.all(np.abs(a[idm,idm:] - var[idm]) < 1.0e-5 * var[idm])

        # Now the real FDMT run.

        del a   # free memory before re-running FDMT
        a = bz_fdmt.FDMT(self.in_buf, self.freq_lo_MHz, self.freq_hi_MHz, self.ndm, np.float32, Verbose=False)

        # Now we just want to find the (DM,time) of the most significant pulse.

        for idm in xrange(self.ndm):
            it = int(np.argmax(a[idm,idm:])) + idm
            sigma = a[idm,it] / var[idm]**0.5

            if (idm == 0) or (sigma > self.max_trigger):
                self.max_trigger = sigma
                self.max_trigger_dm = idm / self.dm0
                self.max_trigger_tfinal = it * self.dt_sample

        self.runflag = True
        

    def _end_pipeline(self, json_outputs):
        if not self.runflag:
            # FIXME had to settle for a warning here instead of an exception, due to a bug in rf_pipelines which I'll fix later!
            # raise RuntimeError("rf_pipelines.bz_fdmt_dedisperser: end_pipeline() called prematurely")
            print >>sys.stderr, "rf_pipelines.bz_fdmt_dedisperser: WARNING: end_pipeline() called prematurely, no search will be run!"

        json_outputs['frb_global_max_trigger'] = self.max_trigger
        json_outputs['frb_global_max_trigger_dm'] = self.max_trigger_dm
        json_outputs['frb_global_max_trigger_tfinal'] = self.max_trigger_tfinal


    def _deallocate(self):
        del self.in_buf
