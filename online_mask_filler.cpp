#include "rf_pipelines_internals.hpp"
#include <algorithm> // for min/max for scalar mask filler
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
// Class for random number generation
// Generates eight random 32-bit floats using a vectorized implementation of xorshift+
// between (-1, 1)
struct test_xorshift_plus
{
  uint64_t s0;
  uint64_t s1;

  test_xorshift_plus(uint64_t _s0 = rd(), uint64_t _s1 = rd()) : s0{_s0}, s1{_s1} {};
  
  inline void gen_floats(float &v1, float &v2)
  {
    uint64_t x = s0;
    uint64_t y = s1;

    s0 = y;
    x ^= (x << 23);
    s1 = x ^ y ^ (x >> 17) ^ (y >> 26);
    
    uint64_t tmp = s1 + y;
    uint32_t tmp0 = tmp; // low 32 bits
    uint32_t tmp1 = tmp >> 32; // high 32
    
    v1 = float(int32_t(tmp0)) * 4.6566129e-10;
    v2 = float(int32_t(tmp1)) * 4.6566129e-10;
  }
};

void print_vec(float *a)
{
    for (int i=0; i < 8; i++)
        cout << a[i] << " ";
    cout << "\n\n";
}

bool test_xorshift(uint64_t i=3289321, uint64_t j=4328934, int niter=100)
{
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

    test_xorshift_plus b(rn1, rn2);
    test_xorshift_plus c(rn3, rn4);
    test_xorshift_plus d(rn5, rn6);
    test_xorshift_plus e(rn7, rn8);
    float srn1, srn2, srn3, srn4, srn5, srn6, srn7, srn8;
    

    for (int iter=0; iter < niter; iter++)
    {
        __m256 vrn = a.gen_floats();
	_mm256_storeu_ps(&vrn_vec[0], vrn);
	b.gen_floats(srn1, srn2);
	c.gen_floats(srn3, srn4);
	d.gen_floats(srn5, srn6);
	e.gen_floats(srn7, srn8);
	float srn_vec[8] = {srn1, srn2, srn3, srn4, srn5, srn6, srn7, srn8};

        for (int i=0; i<8; i++)
	{
	    if (srn_vec[i] != vrn_vec[i])
	    {
  	        cout << "S code outputs: ";
 	        print_vec(srn_vec);
		cout << "V code outputs: ";
		print_vec(vrn_vec);
		cout << "rng test failed: scalar and vectorized prngs are out of sync!" << endl;
		return false;
	    }
	}
    }

    cout << "All rng tests passed." << endl;
    return true;
}


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

shared_ptr<wi_transform> make_online_mask_filler(int v1_chunk, float var_weight, float var_clamp_add, float var_clamp_mult, float w_clamp, float w_cutoff, int nt_chunk)
{
    return make_shared<online_mask_filler> (v1_chunk, var_weight, var_clamp_add, var_clamp_mult, w_clamp, w_cutoff, nt_chunk);
}

shared_ptr<wi_transform> make_scalar_mask_filler(int v1_chunk, float var_weight, float var_clamp_add, float var_clamp_mult, float w_clamp, float w_cutoff, int nt_chunk)
{
    return make_shared<scalar_mask_filler> (v1_chunk, var_weight, var_clamp_add, var_clamp_mult, w_clamp, w_cutoff, nt_chunk);
}



// Externally-visible function for unit testing
void run_online_mask_filler_unit_tests()
{
  test_xorshift();
}


}  // namespace rf_pipelines
