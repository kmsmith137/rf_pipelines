#include "rf_pipelines_internals.hpp"
#include <algorithm> // for min/max for scalar mask filler
#include <sys/time.h>
#include <random>  // for random_device
#include "immintrin.h" // for the rest of the intrinsics

using namespace std;
static random_device rd; // accessible to both random number generators

namespace rf_pipelines {
#if 0
};  // pacify emacs c-mode
#endif


// -------------------------------------------------------------------------------------------------
// Class for random number generation
// Generates eight random 32-bit floats using a vectorized implementation of xorshift+
// between (-1, 1)
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


inline __m256 update_var(__m256 tmp_var, __m256 prev_var, float var_weight, float var_clamp_add, float var_clamp_mult, __m256 mask)
{
    // Does the update of the running variance (tmp_var) by checking the exponential update (normal_upd) doesn't exceed the bounds (high/low) specified 
    // by var_clamp_add and var_clamp_mult
    __m256 normal_upd = _mm256_fmadd_ps(_mm256_set1_ps(var_weight), tmp_var, _mm256_mul_ps(_mm256_set1_ps(1 - var_weight), prev_var)); 
    __m256 high = _mm256_add_ps(_mm256_add_ps(prev_var, _mm256_set1_ps(var_clamp_add)), _mm256_mul_ps(prev_var, _mm256_set1_ps(var_clamp_mult))); 
    __m256 low = _mm256_sub_ps(_mm256_sub_ps(prev_var, _mm256_set1_ps(var_clamp_add)), _mm256_mul_ps(prev_var, _mm256_set1_ps(var_clamp_mult)));     
    __m256 ideal_update = _mm256_max_ps(_mm256_min_ps(normal_upd, high), low); 
    return _mm256_blendv_ps(prev_var, ideal_update, mask);
}


inline void mask_filler(float *intensity, float *weights, vector<float> *running_var, vector<float> *running_weights, ssize_t stride, int nt_chunk, 
			float w_cutoff, float w_clamp, float var_clamp_add, float var_clamp_mult, float var_weight, int nfreq, vec_xorshift_plus &rn)
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
      prev_var = _mm256_set1_ps((*running_var)[ifreq]);
      prev_w = _mm256_set1_ps((*running_weights)[ifreq]);

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
	  __m256 mask = _mm256_cmp_ps(npass, _mm256_set1_ps(-8.1), _CMP_LT_OS);
	      
	  // Here, we do the variance computation:
	  tmp_var = var_est(w0, w1, w2, w3, i0, i1, i2, i3);
	  
	  // Then, use update rules to update value we'll eventually set as our running variance (prevent it from changing too much)
	  tmp_var = update_var(tmp_var, prev_var, var_weight, var_clamp_add, var_clamp_mult, mask);
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
      int *i = reinterpret_cast<int *> (&(*running_var)[ifreq]);   // hack: (int *) pointing to same memory location as q[0]
      int *j = reinterpret_cast<int *> (&(*running_weights)[ifreq]);   
      *i = _mm_extract_ps(y, 0); // write a "fake" int32 to this memory location
      *j = _mm_extract_ps(z, 0);
    }
}


void online_mask_filler::process_chunk(double t0, double t1, float *intensity, float *weights, ssize_t stride, float *pp_intensity, float *pp_weights, ssize_t pp_stride)
{
    mask_filler(intensity, weights, &running_var, &running_weights, stride, nt_chunk, w_cutoff, w_clamp, var_clamp_add, var_clamp_mult, var_weight, nfreq, rn);
}


void online_mask_filler::end_substream()
{
    // Do nothing
}



// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
//
// scalar_mask_filler class
// Everything below exists for the purpose of unit testing everything above! 

// -------------------------------------------------------------------------------------------------
// Class for random number generation
// Generates eight random 32-bit floats using a vectorized implementation of xorshift+
// between (-1, 1)
struct xorshift_plus
{
    vector<uint64_t> seeds;

    xorshift_plus(uint64_t _s0 = rd(), uint64_t _s1 = rd(), 
		  uint64_t _s2 = rd(), uint64_t _s3 = rd(), 
		  uint64_t _s4 = rd(), uint64_t _s5 = rd(), 
		  uint64_t _s6 = rd(), uint64_t _s7 = rd())
      : seeds{_s0, _s1, _s2, _s3, _s4, _s5, _s6, _s7} {};
  
    inline void gen_floats(float *rn)
    {
        // Generates 8 random floats and stores in rn
        for (int i=0; i<8; i+=2)
	{
	    uint64_t x = seeds[i];
	    uint64_t y = seeds[i+1];
	
	    seeds[i] = y;
	    x ^= (x << 23);
	    seeds[i+1] = x ^ y ^ (x >> 17) ^ (y >> 26);
	    
	    uint64_t tmp = seeds[i+1] + y;
	    uint32_t tmp0 = tmp; // low 32 bits
	    uint32_t tmp1 = tmp >> 32; // high 32
	    
	    rn[i] = float(int32_t(tmp0)) * 4.6566129e-10;
	    rn[i+1] = float(int32_t(tmp1)) * 4.6566129e-10;
	}
    }

    inline void gen_weights(float *weights, float pfailv1, float pallzero)
    {
        // This exists solely for the unit test of the online mask filler!
        // Essentially, this is gen_floats + 1 (to get the range right for weights)
        // with extra pfailv1 and pallzero parameters that dictate whether a group of 
        // 32 random numbers (it generates 32 random numbers at once) should result
        // in a failed v1 estimate (i.e. >=24 weights less than the cutoff) or whether
        // all weight values should be zero. Currently, the implementation is a little 
        // boneheaded and can definitely be imporved...
      
        // If the first random number generated is less than pallzero * 2, we make it all zero
        // If the first randon number generated is less than pallzero * 2 + pfailv1 * 2 but
        // greater than pallzero * 2, we make sure the v1 fails by distributing at least 
        // 24 zeros throughout the vector
        // Else, we just fill weights randomly!
        
        for (int i=0; i<32; i+=2)
    	{
    	    uint64_t x = seeds[i % 8];
    	    uint64_t y = seeds[(i+1) % 8];
	
    	    seeds[i % 8] = y;
    	    x ^= (x << 23);
    	    seeds[(i+1) % 8] = x ^ y ^ (x >> 17) ^ (y >> 26);
	    
    	    uint64_t tmp = seeds[(i+1) % 8] + y;
    	    uint32_t tmp0 = tmp; // low 32 bits
    	    uint32_t tmp1 = tmp >> 32; // high 32
	    
    	    if (i == 0)
    	    {
    	    	float rn = float(int32_t(tmp0)) * 4.6566129e-10 + 1;
    	    	if (rn < pallzero * 2)
    	    	{
    	    	    // Make all zero
    	    	    for (int j=0; j>32; j++)
    	    	        weights[j] = 0;
    	    	    return;
    	    	}
    	    	else if (rn < pallzero * 2 + pfailv1 * 2)
    	        {
    	    	    // Make the v1 fail -- this is not super good and can be improved! 
    	    	    for (int j=0; j<25; j++)
    	    	        weights[j] = 0;
    	    	    i = 24;
    	    	}
    	    }
	    
    	    weights[i] = float(int32_t(tmp0)) * 4.6566129e-10 + 1;
    	    weights[i+1] = float(int32_t(tmp1)) * 4.6566129e-10 + 1;
    	}
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
    vector<float> running_weights;
    vector<float> running_var;
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


inline bool get_v1(const float *intensity, const float *weights, float &v1)
{
    int zerocount = 0;
    float vsum = 0;
    float wsum = 0;

    for (int i=0; i < 32; ++i)
    {
        // I assume this is okay for checking whether the weight is 0?
        if (weights[i] < 1e-7)
	  ++zerocount;
        vsum += intensity[i] * intensity[i] * weights[i];
        wsum += weights[i];
    }

    wsum = max(wsum, 1.0f);

    // Check whether enough valid values were passed
    if (zerocount >= 23.9)
    {  
	v1 = 0;
        return false;
    }
    v1 = vsum / wsum;
    return true;
}


inline void scalar_mask_fill(int v1_chunk, float *intensity, float *weights, vector<float> &running_var, vector<float> &running_weights, int nt_chunk, int nfreq, ssize_t stride,
			     float w_cutoff, float w_clamp, float var_clamp_add, float var_clamp_mult, float var_weight, xorshift_plus &rand_x)
{
    float v1;   // stores temporary v1 estimate before it is put into running_var    
    float rn[8];

    rf_assert (v1_chunk == 32);

    for (int ifreq=0; ifreq < nfreq; ++ifreq)
    {
        for (int ichunk=0; ichunk < nt_chunk-1; ichunk += v1_chunk)
	{
	    const float *iacc = &intensity[ifreq*stride + ichunk];
	    const float *wacc = &weights[ifreq*stride + ichunk];

	    // Get v1_chunk
	    if (get_v1(iacc, wacc, v1))
	    {
	        // If the v1 was succesful, try to increase the weight, if possible
	        running_weights[ifreq] = min(2.0f, running_weights[ifreq] + w_clamp);
		// Then, restrict the change in variance estimate definted by the clamp parameters
		v1 = min((1 - var_weight) * running_var[ifreq] + var_weight * v1, running_var[ifreq] + var_clamp_add + running_var[ifreq] * var_clamp_mult);
		v1 = max(v1, running_var[ifreq] - var_clamp_add - running_var[ifreq] * var_clamp_mult);
		// Finally, update the running variance
		running_var[ifreq] = v1;
 	    }
	    else
	    {
	        // For an unsuccessful v1, we decrease the weight if possible. We do not modify the running variance
	        running_weights[ifreq] = max(0.0f, running_weights[ifreq] - w_clamp);
	    }

	    // Do the mask filling for a particular frequency using our new variance estimate
	    for (int i=0; i < v1_chunk; ++i)
	    {
	        if (i % 8 == 0)
		    rand_x.gen_floats(rn);
		if (running_weights[ifreq] != 0)
		{
	    	    if (weights[ifreq*stride+i+ichunk] < w_cutoff)
	    	      intensity[ifreq*stride+i+ichunk] = rn[i % 8] * sqrt(3 * running_var[ifreq]);
	    	}
	        weights[ifreq*stride+i+ichunk] = running_weights[ifreq];
	    }
	} // close the frequency loop
    } // close the ichunk loop
}

void scalar_mask_filler::process_chunk(double t0, double t1, float *intensity, float *weights, ssize_t stride, float *pp_intensity, float *pp_weights, ssize_t pp_stride)
{
  // The things are v1_chunk, *intensity, *weights, vector<float> &running_var, vector<float> &running_weights, int nt_chunk, int nfreq, ssize_t stride,
  // float w_cutoff, float w_clamp, float var_clamp_add, float var_clamp_mult, float var_weight, xorshift_plus &rn
  scalar_mask_fill(v1_chunk, intensity, weights, running_var, running_weights, nt_chunk, nfreq, stride, w_cutoff, w_clamp, var_clamp_add, var_clamp_mult, var_weight, rand_x);
}


void scalar_mask_filler::end_substream()
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

shared_ptr<wi_transform> make_scalar_mask_filler(int v1_chunk, float var_weight, float var_clamp_add, float var_clamp_mult, float w_clamp, float w_cutoff, int nt_chunk)
{
    return make_shared<scalar_mask_filler> (v1_chunk, var_weight, var_clamp_add, var_clamp_mult, w_clamp, w_cutoff, nt_chunk);
}




// -------------------------------------------------------------------------------------------------
//
// Unit test code! 
void print_vec(float *a)
{
    for (int i=0; i < 32; i++)
        cout << a[i] << " ";
    cout << "\n\n";
}

inline bool equality_checker(float a, float b, float epsilon)
{
    // I know this isn't great, but I think it should be sufficient
    return abs(a-b) < epsilon;
}

inline bool test_filler(int nfreq, int nt_chunk, float pfailv1, float pallzero, float w_cutoff=0.5, float w_clamp=3e-3, float var_clamp_add=3e-3, 
			float var_clamp_mult=3e-3, float var_weight=2e-3, int niter=1000)
{
    rf_assert (nfreq * nt_chunk % 8 == 0);
    rf_assert (nt_chunk % 8 == 0);
    rf_assert (pfailv1 < 1);
    rf_assert (pfailv1 >=0);
    rf_assert (pallzero < 1);
    rf_assert (pallzero >=0);
    
    for (int iter=0; iter<niter; iter++)
    {
    // First, we randomize the weights and intensity values
    // We need two copies to put through each processing function and compare
    xorshift_plus rn;
    
    float intensity[nfreq * nt_chunk];
    float weights[nfreq * nt_chunk];
    float intensity2[nfreq * nt_chunk];
    float weights2[nfreq * nt_chunk];

    // Generate intensities 8 at a time using vanilla prng
    for (int i=0; i < nfreq * nt_chunk; i+=8)
        rn.gen_floats(intensity + i);

    // Use custom function to generate weights 
    for (int i=0; i < nfreq * nt_chunk; i+=32)
      rn.gen_weights(weights + i, pfailv1, pallzero);

    // Copy
    for (int i=0; i < nfreq * nt_chunk; i++)
    {
        intensity2[i] = intensity[i];
	weights2[i] = weights[i];
    }

    // Now, we generate random values for the running variance and running weights
    // Using the vec_xorshift_plus functions was too much of a hassle due to vector issues
    // so I've opted to use the C++ rng suite
    mt19937 gen (rd());
    // Weights between 0 and 2
    uniform_real_distribution<float> w_dis(0.0, 2.0);
    // Variance between 0.0 and 0.02
    uniform_real_distribution<float> var_dis(0.0, 0.02);
    
    vector<float> running_var(nfreq);
    vector<float> running_var2(nfreq);
    vector<float> running_weights(nfreq);
    vector<float> running_weights2(nfreq);

    // Make two copies
    for (int i=0; i<nfreq; i++)
    {
        running_var[i] = var_dis(gen);
	running_var2[i] = running_var[i];
	running_weights[i] = w_dis(gen);
	running_weights2[i] = running_weights[i];
    }
    
    // As in the prng unit test, we need to ensure both random number generators are initialized with the same seed values!
    float rn1 = rd();
    float rn2 = rd();
    float rn3 = rd();
    float rn4 = rd();
    float rn5 = rd();
    float rn6 = rd();
    float rn7 = rd();
    float rn8 = rd();
    vec_xorshift_plus vec_rn(_mm256_setr_epi64x(rn1, rn3, rn5, rn7), _mm256_setr_epi64x(rn2, rn4, rn6, rn8));
    xorshift_plus sca_rn(rn1, rn2, rn3, rn4, rn5, rn6, rn7, rn8);
    
    // Process away! Note that the double instances of nt_chunk are for the "stride" parameter which is equal to nt_chunk for this test
    mask_filler(intensity, weights, &running_var, &running_weights, nt_chunk, nt_chunk, w_cutoff, w_clamp, var_clamp_add, var_clamp_mult, var_weight, nfreq, vec_rn);
    scalar_mask_fill(32, intensity2, weights2, running_var2, running_weights2, nt_chunk, nfreq, nt_chunk, w_cutoff, w_clamp, var_clamp_add, var_clamp_mult, var_weight, sca_rn);

    // I realize this next bit isn't the most effecient possible way of doing this comparison, but I think this order will be helpful
    // for debugging any future errors! So it's easy to see where things have gone wrong!
    for (int ifreq=0; ifreq<nfreq; ifreq++)
    {
        // Check running variance
        if (!equality_checker(running_var[ifreq], running_var2[ifreq], 10e-8))
	{
    	    cout << "Something's gone wrong! The running variances at frequency " << ifreq << " on iteration " << iter << " are unequal!" << endl;
	    cout << "Scalar output: " << running_var2[ifreq] << "\t\t Vectorized output: " << running_var[ifreq] << endl;
	    return false;
	}

	// Check running weights
	if (!equality_checker(running_weights[ifreq], running_weights2[ifreq], 10e-8))
	{
	    cout << "Something's gone wrong! The running weights at frequency " << ifreq << " on iteration " << iter << " are unequal!" << endl;
	    cout << "Scalar output: " << running_weights2[ifreq] << "\t\t Vectorized output: " << running_weights[ifreq] << endl;
	    return false;
	}

	for (int i=0; i<nt_chunk; i++)
	{
	    // Check intensity
	    if (!equality_checker(intensity[ifreq * nt_chunk + i], intensity2[ifreq * nt_chunk + i], 10e-5))
	    { 
		cout << "Something has gone wrong! The intensity array produced by the scalar mask filler does not match the intensity array produced by the vectorized mask filler!" << endl;
		cout << "Output terminated at time index " << i << " and frequency " << ifreq << " on iteration " << iter << endl;
		cout << "Scalar output: " << intensity2[ifreq * nt_chunk + i] << "\t\t Vectorized output: " << intensity[ifreq * nt_chunk + i] << endl;
		return false;
	    }
	    
	    // Check weights
	    if (!equality_checker(weights[ifreq * nt_chunk + i], weights2[ifreq * nt_chunk + i], 10e-5))
	    {
		cout << "Something has gone wrong! The weights array produced by the scalar mask filler does not match the weights array produced by the vectorized mask filler!" << endl;
		cout << "Output terminated at time index " << i << " and frequency " << ifreq << " on iteration " << iter << endl;
		cout << "Scalar output: " << weights2[ifreq * nt_chunk + i] << "\t\t Vectorized output: " << weights[ifreq * nt_chunk + i] << endl;
		return false;
	    }
	}
    }
    }
    cout << "***online_mask_filler unit test passed!" << endl;
    return true;
}


inline bool test_xorshift(int niter=10000)
{
    for (int iter=0; iter < niter; iter++)
    {
        // Make sure both prngs are initialized with the same random seeds
        float rn1 = rd();
	float rn2 = rd();
	float rn3 = rd();
	float rn4 = rd();
	float rn5 = rd();
	float rn6 = rd();
	float rn7 = rd();
	float rn8 = rd();

	vec_xorshift_plus a(_mm256_setr_epi64x(rn1, rn3, rn5, rn7), _mm256_setr_epi64x(rn2, rn4, rn6, rn8));
	float vrn_vec[8];
	
	xorshift_plus b(rn1, rn2, rn3, rn4, rn5, rn6, rn7, rn8);
	float srn_vec[8];
 
       __m256 vrn = a.gen_floats();
	_mm256_storeu_ps(&vrn_vec[0], vrn);
	b.gen_floats(srn_vec);

        for (int i=0; i<8; i++)
	{
	  if (!equality_checker(srn_vec[i], vrn_vec[i], 10e-7))
	    {
  	        cout << "S code outputs: ";
 	        print_vec(srn_vec);
		cout << "V code outputs: ";
		print_vec(vrn_vec);
		cout << "rng test failed on iteration " << iter << " and index " << i << ": scalar and vectorized prngs are out of sync!" << endl;
		return false;
	    }
	}
    }

    cout << "***vec_xorshift_plus unit test passed!" << endl;
    return true;
}


void run_online_mask_filler_unit_tests()
{
    // Externally-visible function for unit testing
    test_xorshift();
    test_filler(8, 32, 0.20, 0.20);
}
 

}  // namespace rf_pipelines

