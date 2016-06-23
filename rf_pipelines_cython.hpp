#include "rf_pipelines_internals.hpp"

#ifndef _RF_PIPELINES_CYTHON_HPP
#define _RF_PIPELINES_CYTHON_HPP

using namespace std;


// -------------------------------------------------------------------------------------------------
//
// Transform and stream objects


struct _wi_transform {
    shared_ptr<rf_pipelines::wi_transform> p;

    _wi_transform(const shared_ptr<rf_pipelines::wi_transform> &p_) : p(p_)
    {
	rf_assert(p_);
    }
};


struct _wi_stream {
    shared_ptr<rf_pipelines::wi_stream> p;
    vector<shared_ptr<rf_pipelines::wi_transform> > tbuf;

    _wi_stream(const shared_ptr<rf_pipelines::wi_stream> &p_) : p(p_)
    {
	rf_assert(p_);
    }

    int get_nfreq() { return p->nfreq; }
    int get_nt_maxwrite() { return p->nt_maxwrite; }
    double get_freq_lo_MHz() { return p->freq_lo_MHz; }
    double get_freq_hi_MHz() { return p->freq_hi_MHz; }
    double get_dt_sample() { return p->dt_sample; }    

    void clear_transforms() { tbuf.clear(); }

    void add_transform(_wi_transform *t)
    {
	rf_assert(t != nullptr);
	tbuf.push_back(t->p);
    }
    
    void run() { p->run(tbuf); }
};


// ------------------------------------------------------------------------------------------------
//
// Library


_wi_stream *_make_psrfits_stream(string filename)
{
    return new _wi_stream(rf_pipelines::make_psrfits_stream(filename));
}


_wi_transform *_make_simple_detrender(int nt_chunk)
{
    return new _wi_transform(rf_pipelines::make_simple_detrender(nt_chunk));
}


#endif  // _RF_PIPELINES_CYTHON
