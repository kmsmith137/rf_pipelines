# This file is Cython boilerplate.

from libcpp.string cimport string

cdef extern from "rf_pipelines_cython.hpp":
    cdef cppclass _wi_stream:
        int get_nfreq() except +
        int get_nt_maxwrite() except +
        double get_freq_lo_MHz() except +
        double get_freq_hi_MHz() except +
        double get_dt_sample() except +

    _wi_stream *_make_psrfits_stream(string s)
