#!/usr/bin/env python
#
# I threw together this non-systematic test of some basics in the python interface.
# Some day I'll strengthen this to be a systematic unit test!


import sys
import numpy as np
import rf_pipelines

# Globals
g_nfreq = 64
g_nt_tot = 931

def intval(ifreq, it):
    return 1.3128*ifreq + 0.238*it

def wtval(ifreq, it):
    return 0.382*ifreq + 0.1893*it + 1.0

def almost_equal(x, y):
    epsilon = 1.0e-5 * max(1.0, abs(x)+abs(y))
    return abs(x-y) < epsilon


class test_stream(rf_pipelines.py_wi_stream):
    def __init__(self, nt_chunk):
        freq_lo_MHz = 400.
        freq_hi_MHz = 800.
        dt_sample = 1.0e-3

        rf_pipelines.py_wi_stream.__init__(self, g_nfreq, freq_lo_MHz, freq_hi_MHz, dt_sample, nt_chunk)
        self.nt_chunk = nt_chunk

    def stream_body(self, run_state):
        run_state.start_substream(0.0)

        tcurr = 0
        while tcurr < g_nt_tot:
            nt = min(g_nt_tot-tcurr, self.nt_chunk)
            intensity = np.zeros((g_nfreq, nt))
            weights = np.zeros((g_nfreq, nt))
            
            for ifreq in xrange(g_nfreq):
                for it in xrange(nt):
                    intensity[ifreq,it] = intval(ifreq, tcurr+it)
                    weights[ifreq,it] = wtval(ifreq, tcurr+it)
            
            run_state.write(intensity, weights)
            tcurr += nt

        run_state.end_substream()


class test_transform(rf_pipelines.py_wi_transform):
    def __init__(self, nt_chunk, nt_prepad, nt_postpad):
        self.nfreq = g_nfreq
        self.nt_chunk = nt_chunk
        self.nt_prepad = nt_prepad
        self.nt_postpad = nt_postpad


    def set_stream(self, s):
        self.dt_sample = s.dt_sample


    def start_substream(self, isubstream, t0):
        assert isubstream == 0
        self.tcurr = 0

        
    def process_chunk(self, t0, t1, intensity, weights, pp_intensity, pp_weights):
        assert almost_equal(t0, self.tcurr * self.dt_sample)
        assert intensity.shape == weights.shape == (g_nfreq, self.nt_chunk + self.nt_postpad)

        if self.nt_prepad > 0:
            assert pp_intensity.shape == pp_weights.shape == (g_nfreq, self.nt_prepad)
        else:
            assert pp_intensity == pp_weights == None

        for ifreq in xrange(g_nfreq):
            for it in xrange(self.nt_chunk + self.nt_postpad):
                if (self.tcurr + it) < g_nt_tot:
                    assert almost_equal(intensity[ifreq,it], intval(ifreq, self.tcurr + it))
                    assert almost_equal(weights[ifreq,it], wtval(ifreq, self.tcurr + it))
                else:
                    assert weights[ifreq,it] == 0.0

            for it in xrange(self.nt_prepad):
                if 0 <= (self.tcurr + it - self.nt_prepad) < g_nt_tot:
                    assert almost_equal(pp_intensity[ifreq,it], intval(ifreq, self.tcurr + it - self.nt_prepad))
                    assert almost_equal(pp_weights[ifreq,it], wtval(ifreq, self.tcurr + it - self.nt_prepad))
                else:
                    assert pp_weights[ifreq,it] == 0.0

        self.tcurr += self.nt_chunk


s = test_stream(55)
t1 = test_transform(87, 15, 13)
t2 = test_transform(53, 12, 16)
t3 = test_transform(28, 7, 12)
s.run([t1,t2,t3])
print 'Pass!'
