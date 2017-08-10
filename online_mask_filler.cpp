#include "rf_kernels.hpp"
#include "rf_pipelines_internals.hpp"
#include <random>  // for random_device


using namespace std;
static random_device rd; // accessible to both random number generators


namespace rf_pipelines {
#if 0
};  // pacify emacs c-mode
#endif


// -------------------------------------------------------------------------------------------------
// Online mask filler class
//  
struct online_mask_filler : public wi_transform {
    rf_kernels::online_mask_filler_params params{};
    vector<float> running_weights;
    vector<float> running_var;
    uint64_t rng_state[8];
    
    online_mask_filler(int v1_chunk, float var_weight, float var_clamp_add, float var_clamp_mult, float w_clamp, float w_cutoff, int nt_chunk, bool overwrite_on_wt0, bool modify_weights);
 
    // Override pure virtual member functions in the wi_transform base class.
    // The definitions of these functions below will define the behavior of the online_mask_filler.
    virtual void set_stream(const wi_stream &stream) override;
    virtual void start_substream(int isubstream, double t0) override;
    virtual void process_chunk(double t0, double t1, float *intensity, float *weights, ssize_t stride, float *pp_intensity, float *pp_weights, ssize_t pp_stride) override;
    virtual void end_substream() override;
};


online_mask_filler::online_mask_filler(int v1_chunk_, float var_weight_, float var_clamp_add_, float var_clamp_mult_, float w_clamp_, float w_cutoff_, int nt_chunk_, bool overwrite_on_wt0_, bool modify_weights_)
{
    // Initialize members 'name', 'nt_chunk', which are inherited from wi_transform base class.
    this->nt_chunk = nt_chunk_;

    stringstream ss;
    ss << "online_mask_filler(v1_chunk=" << v1_chunk_ 
       << ",var_weight=" << var_weight_
       << ",var_clamp_add=" << var_clamp_add_
       << ",var_clamp_mult=" << var_clamp_mult_
       << ",w_clamp=" << w_clamp_
       << ",w_cutoff=" << w_cutoff_ 
       << ",nt_chunk=" << nt_chunk_
       << ",overwrite_on_wt0=" << nt_chunk_
       << ",modify_weights=" << nt_chunk_ << ")";

    this->name = ss.str();

    rf_assert (v1_chunk_ == 32);
    rf_assert (nt_chunk_ % 32 == 0);
    rf_assert (nt_chunk_ > 0);
    rf_assert (var_weight_ > 0);
    rf_assert (var_clamp_add_ >= 0);
    rf_assert (var_clamp_mult_ >= 0);
    rf_assert (w_clamp_ > 0);
    rf_assert (w_cutoff_ > 0);
    rf_assert (nt_chunk_ > 0);
    
    params.v1_chunk = v1_chunk_;
    params.var_weight = var_weight_;
    params.var_clamp_add = var_clamp_add_;
    params.var_clamp_mult = var_clamp_mult_;
    params.w_clamp = w_clamp_;
    params.w_cutoff = w_cutoff_;
    params.overwrite_on_wt0 = overwrite_on_wt0_;
    params.modify_weights = modify_weights_;

    if (!params.modify_weights)
      cout << "online_mask_filler.cpp warning: online_mask_filler is being run from rf_pipelines. Weights will not be modified and intensities will be multiplied by the running variance and running weights."
	   << " For old behaviour, set modify_weights to True." << endl;

    random_device rd;
    for (int i = 0; i < 8; i++)
        rng_state[i] = rd();
}


void online_mask_filler::set_stream(const wi_stream &stream)
{
    // Initialize wi_transform::nfreq from wi_stream::nfreq.
    this->nfreq = stream.nfreq;
}


void online_mask_filler::start_substream(int isubstream, double t0)
{
    running_var.resize(nfreq);
    running_weights.resize(nfreq);
}


void online_mask_filler::process_chunk(double t0, double t1, float *intensity, float *weights, ssize_t stride, float *pp_intensity, float *pp_weights, ssize_t pp_stride)
{
    rf_kernels::online_mask_fill(params, nfreq, nt_chunk, stride, intensity, weights, &running_var[0], &running_weights[0], rng_state);
}


void online_mask_filler::end_substream()
{
    // Do nothing
}


// -------------------------------------------------------------------------------------------------
// Scalar online mask filler class
//  
struct scalar_mask_filler : public wi_transform {
    // Specified at construction.
    // Note that nt_chunk is a member of the wi_transform base class.
    rf_kernels::online_mask_filler_params params{};
    vector<float> running_weights;
    vector<float> running_var;
    rf_kernels::xorshift_plus rng;
    
    scalar_mask_filler(int v1_chunk, float var_weight, float var_clamp_add, float var_clamp_mult, float w_clamp, float w_cutoff, int nt_chunk, bool overwrite_on_wt0, bool modify_weights);
 
    // Override pure virtual member functions in the wi_transform base class.
    // The definitions of these functions below will define the behavior of the online_mask_filler.
    virtual void set_stream(const wi_stream &stream) override;
    virtual void start_substream(int isubstream, double t0) override;
    virtual void process_chunk(double t0, double t1, float *intensity, float *weights, ssize_t stride, float *pp_intensity, float *pp_weights, ssize_t pp_stride) override;
    virtual void end_substream() override;
};

scalar_mask_filler::scalar_mask_filler(int v1_chunk_, float var_weight_, float var_clamp_add_, float var_clamp_mult_, float w_clamp_, float w_cutoff_, int nt_chunk_, bool overwrite_on_wt0_, bool modify_weights_) 
{
    // Initialize members 'name', 'nt_chunk', which are inherited from wi_transform base class.
    this->nt_chunk = nt_chunk_;

    stringstream ss;
    ss << "scalar_mask_filler(v1_chunk=" << v1_chunk_ 
       << ",var_weight=" << var_weight_
       << ",var_clamp_add=" << var_clamp_add_
       << ",var_clamp_mult=" << var_clamp_mult_
       << ",w_clamp=" << w_clamp_
       << ",w_cutoff=" << w_cutoff_ 
       << ",nt_chunk=" << nt_chunk_
       << ",overwrite_on_wt0=" << nt_chunk_
       << ",modify_weights=" << nt_chunk_ << ")";

    this->name = ss.str();

    rf_assert (nt_chunk_ % v1_chunk_ == 0);
    rf_assert (nt_chunk_ > 0);
    rf_assert (var_weight_ > 0);
    rf_assert (var_clamp_add_ >= 0);
    rf_assert (var_clamp_mult_ >= 0);
    rf_assert (w_clamp_ > 0);
    rf_assert (w_cutoff_ > 0);
    rf_assert (nt_chunk_ > 0);
    
    params.v1_chunk = v1_chunk_;
    params.var_weight = var_weight_;
    params.var_clamp_add = var_clamp_add_;
    params.var_clamp_mult = var_clamp_mult_;
    params.w_clamp = w_clamp_;
    params.w_cutoff = w_cutoff_;
    params.overwrite_on_wt0 = overwrite_on_wt0_;
    params.modify_weights = modify_weights_;

    if (!params.modify_weights)
      cout << "scalar_mask_filler.cpp warning: online_mask_filler is being run from rf_pipelines. Weights will not be modified and intensities will be multiplied by the running variance and running weights."
	   << " For old behaviour, set modify_weights to True." << endl;
}


void scalar_mask_filler::set_stream(const wi_stream &stream)
{
    // Initialize wi_transform::nfreq from wi_stream::nfreq.
    this->nfreq = stream.nfreq;
}


void scalar_mask_filler::start_substream(int isubstream, double t0)
{
    // All are initialized to 0 this way
    running_var.resize(nfreq);
    running_weights.resize(nfreq);
}

void scalar_mask_filler::process_chunk(double t0, double t1, float *intensity, float *weights, ssize_t stride, float *pp_intensity, float *pp_weights, ssize_t pp_stride)
{
    rf_kernels::scalar_online_mask_fill(params, nfreq, nt_chunk, stride, intensity, weights, &running_var[0], &running_weights[0], rng);
}


void scalar_mask_filler::end_substream()
{
    // Do nothing
}




// -------------------------------------------------------------------------------------------------
//
// Externally-visible factory function: returns pointer to newly constructed online_mask_filler_object

shared_ptr<wi_transform> make_online_mask_filler(int v1_chunk, float var_weight, float var_clamp_add, float var_clamp_mult, float w_clamp, float w_cutoff, int nt_chunk, bool overwrite_on_wt0, bool modify_weights)
{
    return shared_ptr<online_mask_filler> (new online_mask_filler(v1_chunk, var_weight, var_clamp_add, var_clamp_mult, w_clamp, w_cutoff, nt_chunk, overwrite_on_wt0, modify_weights));
}

shared_ptr<wi_transform> make_scalar_mask_filler(int v1_chunk, float var_weight, float var_clamp_add, float var_clamp_mult, float w_clamp, float w_cutoff, int nt_chunk, bool overwrite_on_wt0, bool modify_weights)
{
    return make_shared<scalar_mask_filler> (v1_chunk, var_weight, var_clamp_add, var_clamp_mult, w_clamp, w_cutoff, nt_chunk, overwrite_on_wt0, modify_weights);
}


}  // namespace rf_pipelines

