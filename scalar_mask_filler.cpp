#include "rf_pipelines_internals.hpp"
#include <algorithm> // for max/min
#include <cmath> // for sqrt
#include <random> // for random_device, mt_19937, and random_distribution
#include <sys/time.h>

using namespace std;

namespace rf_pipelines {
#if 0
};  // pacify emacs c-mode
#endif


// -------------------------------------------------------------------------------------------------
//
// scalar_mask_filler class
//
// Note: the online_mask_filler class declaration, and definitions of member functions
// are local to this file, but at the end we define the externally-visible factory
// function make_online_mask_filler(), which returns a pointer to a new online_mask_filler
// object.
//
// Recommended reading: the declaration of 'struct wi_transform' in rf_pipelines.hpp
// and comments contained therein.  This will explain (I hope!) what needs to be implemented,
// for example in scalar_mask_filler::process_chunk().

struct xorshift_plus {
  // Initial seed given by arbitrary hardcoded constants.
  // This is fine for timing, but should come from a std::random_device in production.
  uint64_t s0 = 3289321;
  uint64_t s1 = 4328934;

  // Returns a random integer in the range (0, 2^64-1).
  inline uint64_t gen_u64()
  {
    uint64_t x = s0;
    uint64_t y = s1;
    s0 = y;
    x ^= (x << 23);
    s1 = x ^ y ^ (x >> 17) ^ (y >> 26);
    return s1 + y;
  }

  // Returns a random float in the range (0,1).
  inline float gen_float()
  {
    // The prefactor here is 2^(-64).
    return 5.421010862e-20f * float(gen_u64());
  }
};


struct scalar_mask_filler : public wi_transform {
    // Specified at construction.
    // Note that nt_chunk is a member of the wi_transform base class.
    const int v1_chunk;
    const float var_weight;
    const float var_clamp_add;
    const float var_clamp_mult;
    const float w_clamp;
    const float w_cutoff;
    vector<double> running_weights;
    vector<double> running_var;
  //    std::mt19937 mt_rand;
    xorshift_plus rand_x;
    
    scalar_mask_filler(int v1_chunk, float var_weight, float var_clamp_add, float var_clamp_mult, float w_clamp, float w_cutoff, int nt_chunk);
 
    // Override pure virtual member functions in the wi_transform base class.
    // The definitions of these functions below will define the behavior of the online_mask_filler.
    virtual void set_stream(const wi_stream &stream) override;
    virtual void start_substream(int isubstream, double t0) override;
    virtual void process_chunk(double t0, double t1, float *intensity, float *weights, ssize_t stride, float *pp_intensity, float *pp_weights, ssize_t pp_stride) override;
    virtual void end_substream() override;
};


scalar_mask_filler::scalar_mask_filler(int v1_chunk_, float var_weight_, float var_clamp_add_, float var_clamp_mult_, float w_clamp_, float w_cutoff_, int nt_chunk_) :
    v1_chunk(v1_chunk_),
    var_weight(var_weight_),
    var_clamp_add(var_clamp_add_),
    var_clamp_mult(var_clamp_mult_),
    w_clamp(w_clamp_),
    w_cutoff(w_cutoff_)
{
    // Initialize members 'name', 'nt_chunk', which are inherited from wi_transform base class.

    this->nt_chunk = nt_chunk_;

    stringstream ss;
    ss << "scalar_mask_filler(v1_chunk=" << v1_chunk 
       << ",var_weight=" << var_weight
       << ",var_clamp_add=" << var_clamp_add
       << ",var_clamp_mult=" << var_clamp_mult
       << ",w_clamp=" << w_clamp
       << ",w_cutoff=" << w_cutoff 
       << ",nt_chunk=" << nt_chunk << ")";

    this->name = ss.str();

    //std::random_device rd; 
    //std::mt19937 mt_rand(rd());

    rf_assert (nt_chunk % v1_chunk == 0);
    rf_assert (nt_chunk > 0);
    rf_assert (var_weight > 0);
    rf_assert (var_clamp_add >= 0);
    rf_assert (var_clamp_mult >= 0);
    rf_assert (w_clamp > 0);
    rf_assert (w_cutoff > 0);
    rf_assert (nt_chunk > 0);
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


inline bool get_v1(const float *intensity, const float *weights, double &v1, int v1_chunk)
{
    int zerocount = 0;
    double vsum = 0;
    double wsum = 0;

    for (int i=0; i < v1_chunk; ++i)
    {
        // I assume this is okay for checking whether the weight is 0?
        if (weights[i] < 1e-7)
	  ++zerocount;
        vsum += intensity[i] * intensity[i] * weights[i];
        wsum += weights[i];
    }

    // Check whether enough valid values were passed
    if (zerocount > v1_chunk * 0.75)
    {  
	v1 = 0;
        return false;
    }
    v1 = vsum / wsum;
    return true;
}


void scalar_mask_filler::process_chunk(double t0, double t1, float *intensity, float *weights, ssize_t stride, float *pp_intensity, float *pp_weights, ssize_t pp_stride)
{
    double v1;   // stores temporary v1 estimate before it is put into running_var    

    for (int ichunk=0; ichunk < nt_chunk-1; ichunk += v1_chunk)
    {
        for (int ifreq=0; ifreq < nfreq; ++ifreq)
        {
	    const float *iacc = &intensity[ifreq*stride + ichunk];
	    const float *wacc = &weights[ifreq*stride + ichunk];

	    // Get v1_chunk
	    if (get_v1(iacc, wacc, v1, v1_chunk))
	    {
	        // If the v1 was succesful, try to increase the weight, if possible
	        running_weights[ifreq] = min(2.0, running_weights[ifreq] + w_clamp);
	        // Then, restrict the change in variance estimate definted by the clamp parameters
	        v1 = min((1 - var_weight) * running_var[ifreq] + var_weight * v1, running_var[ifreq] + var_clamp_add + running_var[ifreq] * var_clamp_mult);
	        v1 = max(v1, running_var[ifreq] - var_clamp_add - running_var[ifreq] * var_clamp_mult);
	        // Finally, update the running variance
	        running_var[ifreq] = v1;
	    }
	    else
	    {
	        // For an unsuccessful v1, we decrease the weight if possible. We do not modify the running variance
	        running_weights[ifreq] = max(0.0, running_weights[ifreq] - w_clamp);
	    }
	      
	    // Do the mask filling for a particular frequency using our new variance estimate
	    //std::uniform_real_distribution<double> dist(-sqrt(3*running_var[ifreq]), sqrt(3*running_var[ifreq]));
	    //std::normal_distribution<double> dist(0, sqrt(running_var[ifreq]));
	    for (int i=0; i < v1_chunk; ++i)
	    {
	        if (running_weights[ifreq] != 0)
		{
		    if (weights[ifreq*stride+i+ichunk] < w_cutoff)
		      intensity[ifreq*stride+i+ichunk] = rand_x.gen_float() * 2 * (sqrt(3 * running_var[ifreq])) - sqrt(3 * running_var[ifreq]);
		}
	        weights[ifreq*stride+i+ichunk] = running_weights[ifreq];
	    }
	} // close the frequency loop
    } // close the ichunk loop
}


void scalar_mask_filler::end_substream()
{
    // Do nothing
}


// -------------------------------------------------------------------------------------------------
//
// Externally-visible factory function: returns pointer to newly constructed online_mask_filler_object


shared_ptr<wi_transform> make_scalar_mask_filler(int v1_chunk, float var_weight, float var_clamp_add, float var_clamp_mult, float w_clamp, float w_cutoff, int nt_chunk)
{
    return make_shared<scalar_mask_filler> (v1_chunk, var_weight, var_clamp_add, var_clamp_mult, w_clamp, w_cutoff, nt_chunk);
}


}  // namespace rf_pipelines
