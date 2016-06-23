#include "rf_pipelines_internals.hpp"

#ifndef _RF_PIPELINES_CYTHON_HPP
#define _RF_PIPELINES_CYTHON_HPP

using namespace std;

struct _wi_stream {
    std::shared_ptr<rf_pipelines::wi_stream> p;

    _wi_stream(const shared_ptr<rf_pipelines::wi_stream> &p_) : p(p_)
    {
	rf_assert(p_);
    }

    int get_nfreq() { return p->nfreq; }
    int get_nt_maxwrite() { return p->nt_maxwrite; }
    double get_freq_lo_MHz() { return p->freq_lo_MHz; }
    double get_freq_hi_MHz() { return p->freq_hi_MHz; }
    double get_dt_sample() { return p->dt_sample; }    
};

_wi_stream *_make_psrfits_stream(string filename)
{
    return new _wi_stream(rf_pipelines::make_psrfits_stream(filename));
}


#endif  // _RF_PIPELINES_CYTHON
