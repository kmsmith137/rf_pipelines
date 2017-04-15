#include "rf_pipelines_internals.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
};  // pacify emacs c-mode
#endif


// -------------------------------------------------------------------------------------------------
//
// online_mask_filler class
//
// Note: the online_mask_filler class declaration, and definitions of member functions
// are local to this file, but at the end we define the externally-visible factory
// function make_online_mask_filler(), which returns a pointer to a new online_mask_filler
// object.
//
// All member functions are currently placeholders!  All that is currently implemented is the
// bare minimum needed to run the pipeline without failing an assert (for example, we need to
// initialize online_mask_filler::nfreq).
//
// Recommended reading: the declaration of 'struct wi_transform' in rf_pipelines.hpp
// and comments contained therein.  This will explain (I hope!) what needs to be implemented,
// for example in online_mask_filler::process_chunk().


struct online_mask_filler : public wi_transform {
    // Specified at construction.
    // Note that nt_chunk is a member of the wi_transform base class.
    const int v1_chunk;
    const int v2_chunk;
    const float w_cutoff;

    // More fields can be added here!  ('v1_estimates', 'running_var', etc.)
    
    online_mask_filler(int v1_chunk, int v2_chunk, float w_cutoff, int nt_chunk);

    // Override pure virtual member functions in the wi_transform base class.
    // The definitions of these functions below will define the behavior of the online_mask_filler.
    virtual void set_stream(const wi_stream &stream) override;
    virtual void start_substream(int isubstream, double t0) override;
    virtual void process_chunk(double t0, double t1, float *intensity, float *weights, ssize_t stride, float *pp_intensity, float *pp_weights, ssize_t pp_stride) override;
    virtual void end_substream() override;
};


online_mask_filler::online_mask_filler(int v1_chunk_, int v2_chunk_, float w_cutoff_, int nt_chunk_) :
    v1_chunk(v1_chunk_),
    v2_chunk(v2_chunk_),
    w_cutoff(w_cutoff_)
{
    // Initialize members 'name', 'nt_chunk', which are inherited from wi_transform base class.

    stringstream ss;
    ss << "online_mask_filler(v1_chunk=" << v1_chunk << ",v2_chunk=" << v2_chunk
       << ",w_cutoff=" << w_cutoff << ",nt_chunk=" << nt_chunk_ << ",cpp=True)";

    this->name = ss.str();
    this->nt_chunk = nt_chunk_;
}


void online_mask_filler::set_stream(const wi_stream &stream)
{
    // Initialize wi_transform::nfreq from wi_stream::nfreq.
    this->nfreq = stream.nfreq;
}


void online_mask_filler::start_substream(int isubstream, double t0)
{
    // To be filled in!
}


void online_mask_filler::process_chunk(double t0, double t1, float *intensity, float *weights, ssize_t stride, float *pp_intensity, float *pp_weights, ssize_t pp_stride)
{
    // To be filled in!
}


void online_mask_filler::end_substream()
{
    // To be filled in!
}


// -------------------------------------------------------------------------------------------------
//
// Externally-visible factory function: returns pointer to newly constructed online_mask_filler_object


shared_ptr<wi_transform> make_online_mask_filler(int v1_chunk, int v2_chunk, float w_cutoff, int nt_chunk)
{
    return make_shared<online_mask_filler> (v1_chunk, v2_chunk, w_cutoff, nt_chunk);
}


}  // namespace rf_pipelines