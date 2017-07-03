#include "rf_pipelines_internals.hpp"
#include <sys/time.h>
#include <random>  // for random_device
#include "immintrin.h" // for the rest of the intrinsics

using namespace std;

namespace rf_pipelines {
#if 0
};  // pacify emacs c-mode
#endif


// -------------------------------------------------------------------------------------------------
// Class for random number generation
// Generates eight random 32-bit floats using a vectorized implementation of xorshift+
// between (-1, 1)
static random_device rd;
struct vec_xorshift_plus
{
  // Seed values
  __m256i s0; 
  __m256i s1;

  // Initialize seeds to random device, unless alternate seeds are specified
  vec_xorshift_plus(__m256i _s0 = _mm256_setr_epi64x(rd(), rd(), rd(), rd()), __m256i _s1 = _mm256_setr_epi64x(rd(), rd(), rd(), rd())) : s0(_s0), s1(_s1) {};

  // Generates 256 random bits (interpreted as 8 signed floats)
  // Returns an __m256 vector, so bits must be stored using _mm256_storeu_ps() intrinsic!
  inline __m256 gen_floats()
  {
    // x = s0
    __m256i x = s0;
    // y = s1
    __m256i y = s1;
    // s0 = y
    s0 = y;
    // x ^= (x << 23)
    x = _mm256_xor_si256(x, _mm256_slli_epi64(x, 23));
    // s1 = x ^ y ^ (x >> 17) ^ (y >> 26)
    s1 = _mm256_xor_si256(x, y);
    s1 = _mm256_xor_si256(s1, _mm256_srli_epi64(x, 17));
    s1 = _mm256_xor_si256(s1, _mm256_srli_epi64(y, 26));
      
    // Convert to 8 signed 32-bit floats in range (-1, 1), since we multiply by 
    // a prefactor of 2^(-31)
    return _mm256_mul_ps(_mm256_cvtepi32_ps(_mm256_add_epi64(y, s1)), _mm256_set1_ps(4.6566129e-10));
  }
};

// -------------------------------------------------------------------------------------------------
// Online mask filler class
//  
struct online_mask_filler : public wi_transform {
    const int v1_chunk;
    const float var_weight;
    const float var_clamp_add;
    const float var_clamp_mult;
    const float w_clamp;
    const float w_cutoff;
    vector<float> running_weights;
    vector<float> running_var;
    vec_xorshift_plus rn;
    
    online_mask_filler(int v1_chunk, float var_weight, float var_clamp_add, float var_clamp_mult, float w_clamp, float w_cutoff, int nt_chunk);
 
    // Override pure virtual member functions in the wi_transform base class.
    // The definitions of these functions below will define the behavior of the online_mask_filler.
    virtual void set_stream(const wi_stream &stream) override;
    virtual void start_substream(int isubstream, double t0) override;
    virtual void process_chunk(double t0, double t1, float *intensity, float *weights, ssize_t stride, float *pp_intensity, float *pp_weights, ssize_t pp_stride) override;
    virtual void end_substream() override;
};


online_mask_filler::online_mask_filler(int v1_chunk_, float var_weight_, float var_clamp_add_, float var_clamp_mult_, float w_clamp_, float w_cutoff_, int nt_chunk_) :
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
    ss << "online_mask_filler(v1_chunk=" << v1_chunk 
       << ",var_weight=" << var_weight
       << ",var_clamp_add=" << var_clamp_add
       << ",var_clamp_mult=" << var_clamp_mult
       << ",w_clamp=" << w_clamp
       << ",w_cutoff=" << w_cutoff 
       << ",nt_chunk=" << nt_chunk << ")";

    this->name = ss.str();

    rf_assert (nt_chunk % 32 == 0);
    rf_assert (nt_chunk > 0);
    rf_assert (var_weight > 0);
    rf_assert (var_clamp_add >= 0);
    rf_assert (var_clamp_mult >= 0);
    rf_assert (w_clamp > 0);
    rf_assert (w_cutoff > 0);
    rf_assert (nt_chunk > 0);
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


inline __m256 hadd(__m256 a)
{
    // Does a horizontal add of a __m256 register
    __m256 tmp0 = _mm256_add_ps(_mm256_permute2f128_ps(a, a, 0b00100001), a);
    __m256 tmp1 = _mm256_add_ps(_mm256_permute_ps(tmp0, 0b00111001), tmp0);
    return _mm256_add_ps(_mm256_permute_ps(tmp1, 0b01001110), tmp1);
}


inline __m256 check_weights(__m256 w0, __m256 w1, __m256 w2, __m256 w3)
{
    // Check the number of weights that are non-zero
    // Returns a constant register containing the number of _successful_ weights
    __m256 zero = _mm256_set1_ps(0.0);

    // Note here that _mm256_cmp_ps returns a _m256 in which groups of 8 bits are either all 1s if the cmp was true and all 0s if the cmp was false
    // If we reinterpret this as an _m256i, we get a register in which all values are either 0 (in the case of 8 zeros) or -1 (in the case of 8 ones
    // given that 11111111 is -1 in twos complement form). Thus, by doing an hadd at the end, we will get -1 * the number of successful intensities
    // Note that _mm256_castps_si256 just reinterprets bits whereas _mm256_cvtepi32_ps truncates a float into an int. 
    __m256i w0_mask = _mm256_castps_si256(_mm256_cmp_ps(w0, zero, _CMP_GT_OS));
    __m256i w1_mask = _mm256_castps_si256(_mm256_cmp_ps(w1, zero, _CMP_GT_OS));
    __m256i w2_mask = _mm256_castps_si256(_mm256_cmp_ps(w2, zero, _CMP_GT_OS));
    __m256i w3_mask = _mm256_castps_si256(_mm256_cmp_ps(w3, zero, _CMP_GT_OS));
    return hadd(_mm256_cvtepi32_ps(_mm256_add_epi32(_mm256_add_epi32(_mm256_add_epi32(w0_mask, w1_mask), w2_mask), w3_mask)));
}


inline __m256 var_est(__m256 w0, __m256 w1, __m256 w2, __m256 w3, __m256 i0, __m256 i1, __m256 i2, __m256 i3)
{
    // Does variance estimation for 32 intensity/weight values. Assumes mean=0.
    __m256 wi0 = _mm256_mul_ps(_mm256_mul_ps(i0, i0), w0);
    __m256 wi01 = _mm256_fmadd_ps(_mm256_mul_ps(i1, i1), w1, wi0);
    __m256 wi012 = _mm256_fmadd_ps(_mm256_mul_ps(i2, i2), w2, wi01);
    __m256 wi0123 = _mm256_fmadd_ps(_mm256_mul_ps(i3, i3), w3, wi012);
    __m256 vsum = hadd(wi0123);
    __m256 wsum = hadd(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(w0, w1), w2), w3));
    wsum = _mm256_max_ps(wsum, _mm256_set1_ps(1.0));
    return _mm256_div_ps(vsum, wsum);
}


inline __m256 update_var(__m256 tmp_var, __m256 prev_var, float var_weight, float var_clamp_add, float var_clamp_mult)
{
    // Does the update of the running variance (tmp_var) by checking the exponential update (normal_upd) doesn't exceed the bounds (high/low) specified 
    // by var_clamp_add and var_clamp_mult
    __m256 normal_upd = _mm256_fmadd_ps(_mm256_set1_ps(var_weight), tmp_var, _mm256_mul_ps(_mm256_set1_ps(1 - var_weight), prev_var)); 
    __m256 high = _mm256_add_ps(_mm256_add_ps(prev_var, _mm256_set1_ps(var_clamp_add)), _mm256_mul_ps(prev_var, _mm256_set1_ps(var_clamp_mult))); 
    __m256 low = _mm256_sub_ps(_mm256_sub_ps(prev_var, _mm256_set1_ps(var_clamp_add)), _mm256_mul_ps(prev_var, _mm256_set1_ps(var_clamp_mult))); 
    return _mm256_max_ps(_mm256_min_ps(normal_upd, high), low);
}


inline void print_arr(__m256 a)
{
    // Helper function to print __m256 register (float[8])
    float arr[8];
    _mm256_storeu_ps((float *) &arr, a);
    for (int i=0; i<8; ++i)
        cout << arr[i] << " ";
    cout << "\n";
}


inline void print_arri(__m256i a)
{
    // Helper function to print __m256i register (int[8])
    int arr[8];
    _mm256_storeu_si256((__m256i *) &arr, a);
    for (int i=0; i<8; ++i)
        cout << arr[i] << " ";
    cout << "\n";
}


void online_mask_filler::process_chunk(double t0, double t1, float *intensity, float *weights, ssize_t stride, float *pp_intensity, float *pp_weights, ssize_t pp_stride)
{
    __m256 tmp_var, prev_var, prev_w, w0, w1, w2, w3, i0, i1, i2, i3, res0, res1, res2, res3;
    __m256 c = _mm256_set1_ps(w_cutoff);
    __m256 root_three = _mm256_sqrt_ps(_mm256_set1_ps(3));
    __m256 two = _mm256_set1_ps(2);
    __m256 zero = _mm256_set1_ps(0.0);

    // Loop over frequencies first to avoid having to write the running_var and running_weights in each iteration of the ichunk loop
    for (int ifreq=0; ifreq<nfreq; ifreq++)
    {
        // Get the previous running_var and running_weights
        prev_var = _mm256_set1_ps(running_var[ifreq]);
	prev_w = _mm256_set1_ps(running_weights[ifreq]);

	for (int ichunk=0; ichunk<nt_chunk-1; ichunk += 32)
	{
	    // Load intensity and weight arrays
	    i0 = _mm256_loadu_ps(intensity + ifreq * stride + ichunk);
	    i1 = _mm256_loadu_ps(intensity + ifreq * stride + ichunk + 8);
	    i2 = _mm256_loadu_ps(intensity + ifreq * stride + ichunk + 16);
	    i3 = _mm256_loadu_ps(intensity + ifreq * stride + ichunk + 24);
	    w0 = _mm256_loadu_ps(weights + ifreq * stride + ichunk);
	    w1 = _mm256_loadu_ps(weights + ifreq * stride + ichunk + 8);
	    w2 = _mm256_loadu_ps(weights + ifreq * stride + ichunk + 16);
	    w3 = _mm256_loadu_ps(weights + ifreq * stride + ichunk + 24);

	    // First, we need to see how many of the weights are greater than zero. Note that npass contains 
	    // a number equal to -1 * the number of successful intensities (see check_weights comments)
	    __m256 npass = check_weights(w0, w1, w2, w3);

	    // If pass is less than -8, we treat this as a failed v1 case (since >75% of the data is masked, we can't use that to update 
	    // our running variance estimate). After doing the compare below, mask will contain we get a constant register that is either 
	    // all zeros if the variance estimate failed or all ones if the variance estimate was successful.
	    __m256 mask = _mm256_cmp_ps(npass, _mm256_set1_ps(-8), _CMP_LT_OS);
	    
	    // Here, we do the variance computation:
	    tmp_var = var_est(w0, w1, w2, w3, i0, i1, i2, i3);
	    
	    // Then, use update rules to update value we'll eventually set as our running variance (prevent it from changing too much)
	    tmp_var = update_var(tmp_var, prev_var, var_weight, var_clamp_add, var_clamp_mult);
	    prev_var = tmp_var;

	    // Finally, mask fill with the running variance -- if weights less than cutoff, fill
	    res0 = _mm256_blendv_ps(i0, _mm256_mul_ps(rn.gen_floats(), _mm256_mul_ps(root_three, _mm256_sqrt_ps(tmp_var))), _mm256_cmp_ps(w0, c, _CMP_LT_OS));
	    res1 = _mm256_blendv_ps(i1, _mm256_mul_ps(rn.gen_floats(), _mm256_mul_ps(root_three, _mm256_sqrt_ps(tmp_var))), _mm256_cmp_ps(w1, c, _CMP_LT_OS));
	    res2 = _mm256_blendv_ps(i2, _mm256_mul_ps(rn.gen_floats(), _mm256_mul_ps(root_three, _mm256_sqrt_ps(tmp_var))), _mm256_cmp_ps(w2, c, _CMP_LT_OS));
	    res3 = _mm256_blendv_ps(i3, _mm256_mul_ps(rn.gen_floats(), _mm256_mul_ps(root_three, _mm256_sqrt_ps(tmp_var))), _mm256_cmp_ps(w3, c, _CMP_LT_OS));
	    
	    // We also need to modify the weight values
	    __m256 w = _mm256_blendv_ps(_mm256_set1_ps(-w_clamp), _mm256_set1_ps(w_clamp), mask);    // either +w_clamp or -w_clamp
	    w = _mm256_min_ps(_mm256_max_ps(_mm256_add_ps(w, prev_w), zero), two);
	    prev_w = w;

	    // Store the new intensity values
	    _mm256_storeu_ps((float*) (intensity + ifreq * stride + ichunk), res0);
	    _mm256_storeu_ps((float*) (intensity + ifreq * stride + ichunk + 8), res1);
	    _mm256_storeu_ps((float*) (intensity + ifreq * stride + ichunk + 16), res2);
	    _mm256_storeu_ps((float*) (intensity + ifreq * stride + ichunk + 24), res3);
	    
	    // Store the new weight values
	    _mm256_storeu_ps((float*) (weights + ifreq * stride + ichunk), w);
	    _mm256_storeu_ps((float*) (weights + ifreq * stride + ichunk + 8), w);
	    _mm256_storeu_ps((float*) (weights + ifreq * stride + ichunk + 16), w);
	    _mm256_storeu_ps((float*) (weights + ifreq * stride + ichunk + 24), w);
	}
	// Since we've now completed all the variance estimation and filling for this frequency channel in this chunk, we must write our 
	// running variance and weight to the vector, which is a bit of a pain. Thanks for this hack, Kendrick!

	// First step: extract elements 0-3 into a 128-bit register.
	__m128 y = _mm256_extractf128_ps(prev_var, 0);
	__m128 z = _mm256_extractf128_ps(prev_w, 0);

	// The intrinsic _mm_extract_ps() extracts element 0 from the 128-bit register, but it has the wrong
	// return type (int32 instead of float32). The returned value is a "fake" int32 obtained by interpreting 
	// the bit pattern of the "real" float32 (in IEEE-754 representation) as an int32 (in twos-complement representation),
	// so it's not very useful. Nevertheless if we just write it to memory as an int32, and read back from 
	// the same memory location later as a float32, we'll get the right answer.
	int *i = reinterpret_cast<int *> (&running_var[ifreq]);   // hack: (int *) pointing to same memory location as q[0]
	int *j = reinterpret_cast<int *> (&running_weights[ifreq]);   
	*i = _mm_extract_ps(y, 0); // write a "fake" int32 to this memory location
	*j = _mm_extract_ps(z, 0);
    }
}


void online_mask_filler::end_substream()
{
    // Do nothing
}


// -------------------------------------------------------------------------------------------------
//
// Externally-visible factory function: returns pointer to newly constructed online_mask_filler_object


shared_ptr<wi_transform> make_online_mask_filler(int v1_chunk, float var_weight, float var_clamp_add, float var_clamp_mult, float w_clamp, float w_cutoff, int nt_chunk)
{
    return make_shared<online_mask_filler> (v1_chunk, var_weight, var_clamp_add, var_clamp_mult, w_clamp, w_cutoff, nt_chunk);
}


// Externally-visible function for unit testing
void run_online_mask_filler_unit_tests()
{
    cout << "run_online_mask_filler_unit_tests: placeholder for now!\n";
}


}  // namespace rf_pipelines
