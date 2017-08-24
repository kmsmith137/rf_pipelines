#include "rf_kernels/online_mask_filler.hpp"
#include "rf_pipelines_internals.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
};  // pacify emacs c-mode
#endif


// -------------------------------------------------------------------------------------------------
// Online mask filler class
//  
struct online_mask_filler : public wi_transform {
    const int v1_chunk;
    const float var_weight;
    const float w_clamp;
    const float w_cutoff;
    const bool modify_weights;
    const bool multiply_intensity_by_weights;
    const bool use_scalar_kernel;
    
    std::shared_ptr<rf_kernels::online_mask_filler> kernel;
    
    online_mask_filler(int v1_chunk, float var_weight, float w_clamp, float w_cutoff, int nt_chunk, 
		       bool overwrite_on_wt0, bool modify_weights, bool use_scalar_kernel);
 
    // Override pure virtual member functions in the wi_transform base class.
    // The definitions of these functions below will define the behavior of the online_mask_filler.
    virtual void set_stream(const wi_stream &stream) override;
    virtual void start_substream(int isubstream, double t0) override;
    virtual void process_chunk(double t0, double t1, float *intensity, float *weights, ssize_t stride, float *pp_intensity, float *pp_weights, ssize_t pp_stride) override;
    virtual void end_substream() override;
};


online_mask_filler::online_mask_filler(int v1_chunk_, float var_weight_, float w_clamp_, float w_cutoff_, int nt_chunk_, 
				       bool modify_weights_, bool multiply_intensity_by_weights_, bool use_scalar_kernel_) :
				       
    v1_chunk(v1_chunk_),
    var_weight(var_weight_),
    w_clamp(w_clamp_),
    w_cutoff(w_cutoff_),
    modify_weights(modify_weights_),
    multiply_intensity_by_weights(multiply_intensity_by_weights_),
    use_scalar_kernel(use_scalar_kernel_)
{
    // Initialize members 'name', 'nt_chunk', which are inherited from wi_transform base class.
    this->nt_chunk = nt_chunk_;

    stringstream ss;
    ss << "online_mask_filler(v1_chunk=" << v1_chunk
       << ",var_weight=" << var_weight
       << ",w_clamp=" << w_clamp
       << ",w_cutoff=" << w_cutoff
       << ",nt_chunk=" << nt_chunk
       << ",modify_weights=" << modify_weights
       << ",multiply_intensity_by_weights=" << multiply_intensity_by_weights
       << ",use_scalar_kernel=" << use_scalar_kernel
       << ")";

    this->name = ss.str();

    rf_assert (v1_chunk == 32);
    rf_assert (nt_chunk % 32 == 0);
    rf_assert (nt_chunk > 0);
    rf_assert (var_weight > 0);
    rf_assert (w_clamp > 0);
    rf_assert (w_cutoff >= 0);
    rf_assert (nt_chunk > 0);

    if (!modify_weights) {
	cout << "online_mask_filler.cpp warning: online_mask_filler is being run from rf_pipelines."
	     << " Weights will not be modified."
	     << " For old behaviour, set modify_weights to True." << endl;
    }

    if (multiply_intensity_by_weights) {
	cout << "online_mask_filler.cpp warning: online_mask_filler is being run from rf_pipelines."
	     << " Intensities will be multiplied by the running weights."
	     << " For old behaviour, set multiply_intensity_by_weights to False." << endl;
    }
}


void online_mask_filler::set_stream(const wi_stream &stream)
{
    // Initialize wi_transform::nfreq from wi_stream::nfreq.
    this->nfreq = stream.nfreq;

    // Now that nfreq is known, initialize kernel.
    this->kernel = make_unique<rf_kernels::online_mask_filler> (nfreq);

    kernel->v1_chunk = this->v1_chunk;
    kernel->var_weight = this->var_weight;
    kernel->w_clamp = this->w_clamp;
    kernel->w_cutoff = this->w_cutoff;
    kernel->modify_weights = this->modify_weights;
    kernel->multiply_intensity_by_weights = this->multiply_intensity_by_weights;
}


void online_mask_filler::start_substream(int isubstream, double t0)
{
    // Do nothing
}


void online_mask_filler::process_chunk(double t0, double t1, float *intensity, float *weights, ssize_t stride, float *pp_intensity, float *pp_weights, ssize_t pp_stride)
{
    rf_assert(kernel.get() != nullptr);

    if (use_scalar_kernel)
	kernel->scalar_mask_fill(nt_chunk, stride, intensity, weights);
    else
	kernel->mask_fill(nt_chunk, stride, intensity, weights);
}


void online_mask_filler::end_substream()
{
    // Do nothing
}


// -------------------------------------------------------------------------------------------------
//
// Externally-visible factory function: returns pointer to newly constructed online_mask_filler_object

shared_ptr<wi_transform> make_online_mask_filler(int v1_chunk, float var_weight, float w_clamp, float w_cutoff, int nt_chunk, bool modify_weights, bool multiply_intensity_by_weights)
{
    return make_shared<online_mask_filler> (v1_chunk, var_weight, w_clamp, w_cutoff, nt_chunk, modify_weights, multiply_intensity_by_weights, false);
}

shared_ptr<wi_transform> make_scalar_mask_filler(int v1_chunk, float var_weight, float w_clamp, float w_cutoff, int nt_chunk, bool modify_weights, bool multiply_intensity_by_weights)
{
    return make_shared<online_mask_filler> (v1_chunk, var_weight, w_clamp, w_cutoff, nt_chunk, modify_weights, multiply_intensity_by_weights, true);
}


}  // namespace rf_pipelines

