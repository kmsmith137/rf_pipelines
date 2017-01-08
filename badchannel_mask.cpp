#include "rf_pipelines_internals.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif

// This file is a placeholder for a C++ implementation of the 'badchannel_mask' class.
//
// Right now, it just contains some general boilerplate needed to make a C++ badchannel_mask
// class available throughout the pipeline.  In particular, the C++ badchannel_mask can be
// obtained from python as 'rf_pipelines_c.badchannel_mask'.
//
// The actual implementation of the badchannel_mask class is still missing -- right now,
// all member functions throw exceptions.


struct badchannel_mask : public wi_transform {
    // Note: inherits { nfreq, nt_chunk, nt_prepad, nt_postpad } from base class wi_transform

    badchannel_mask(const string &maskpath, int nt_chunk)
    {
	throw runtime_error("oops, C++ badchannel_mask class was requested, but not implemented yet");
    }

    // As explaned in rf_pipelines.hpp, the following four virtual functions in the base class
    // must be overridden, in order to define the badchannel_mask subclass.

    virtual void set_stream(const wi_stream &stream) override
    {
	throw runtime_error("oops, C++ badchannel_mask class was requested, but not implemented yet");
    }

    virtual void start_substream(int isubstream, double t0) override
    {
	throw runtime_error("oops, C++ badchannel_mask class was requested, but not implemented yet");
    }

    virtual void process_chunk(double t0, double t1, float *intensity, float *weights, ssize_t stride, float *pp_intensity, float *pp_weights, ssize_t pp_stride) override
    {
	throw runtime_error("oops, C++ badchannel_mask class was requested, but not implemented yet");
    }

    virtual void end_substream() override
    {
	throw runtime_error("oops, C++ badchannel_mask class was requested, but not implemented yet");
    }
};


shared_ptr<wi_transform> make_badchannel_mask(const string &maskpath, int nt_chunk)
{
    return make_shared<badchannel_mask> (maskpath, nt_chunk);
}


}  // namespace rf_pipelines
