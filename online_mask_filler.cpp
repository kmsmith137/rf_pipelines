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
struct vec_xorshift_plus {
  random_device rd;
  __m256i s0 = _mm256_setr_epi64x(rd(), rd(), rd(), rd());
  __m256i s1 = _mm256_setr_epi64x(rd(), rd(), rd(), rd());
  
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
      
    // Let's convert to 8 signed 32-bit floats and return that
    // This will give us random numbers centered around 0!
    // We also need to multiply by a prefactor of 2^(-31) to get
    // the range to be (-1, 1)
    return _mm256_mul_ps(_mm256_cvtepi32_ps(_mm256_add_epi64(y, s1)), _mm256_set1_ps(4.6566129e-10));
  }
};

// -------------------------------------------------------------------------------------------------
// Online mask filler class
// Does good stuff
struct online_mask_filler : public wi_transform {
    const int v1_chunk;
    const float var_weight;
    const float var_clamp_add;
    const float var_clamp_mult;
    const float w_clamp;
    const float w_cutoff;
    vector<double> running_weights;
    vector<double> running_var;
    vec_xorshift_plus rand_x;
    
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

    rf_assert (nt_chunk % v1_chunk == 0);
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
    // All are initialized to 0 this way
    running_var.resize(nfreq);
    running_weights.resize(nfreq);
}


inline __m256 hadd(__m256 a)
{
    __m256 tmp0 = _mm256_add_ps(_mm256_permute2f128_ps(a, a, 0b00100001), a);
    __m256 tmp1 = _mm256_add_ps(_mm256_permute_ps(tmp0, 0b00111001), tmp0);
    return _mm256_add_ps(_mm256_permute_ps(tmp1, 0b01001110), tmp1);
}


void online_mask_filler::process_chunk(double t0, double t1, float *intensity, float *weights, ssize_t stride, float *pp_intensity, float *pp_weights, ssize_t pp_stride)
{
  __m256 vw = _mm256_set1_ps(var_weight);
  vec_xorshift_plus rn;
  __m256 vsum, wsum, tmp_var, tmp_w, tmp1, tmp2, w0, w1, w2, w3, i0, i1, i2, i3, res0, res1, res2, res3;
  __m256 c = _mm256_set1_ps(w_cutoff);
  __m256 root_three = _mm256_sqrt_ps(_mm256_set1_ps(3)); // handy for later
  __m256 two = _mm256_set1_ps(2);
  __m256 zero = _mm256_set1_ps(0);
    
  for (int ifreq=0; ifreq<nfreq; ifreq++)
  {
    cout << "Begin frequency " << ifreq << endl;
      tmp_var = _mm256_set1_ps(running_var[ifreq]);
      tmp_w = _mm256_set1_ps(running_weights[ifreq]);
      for (int ichunk=0; ichunk<nt_chunk-1; ichunk += 32)
      {
	cout << "Begin chunk " << ichunk << endl;
	  // Load intensity and weight arrays
	i0 = _mm256_load_ps((const float*) (intensity + ifreq * stride + ichunk));
	  cout << "Load literally one thing" << endl;
	  i1 = _mm256_load_ps(intensity + ifreq * stride + ichunk + 8);
	  i2 = _mm256_load_ps(intensity + ifreq * stride + ichunk + 16);
	  i3 = _mm256_load_ps(intensity + ifreq * stride + ichunk + 24);
	  w0 = _mm256_load_ps(weights + ifreq * stride + ichunk);
	  w1 = _mm256_load_ps(weights + ifreq * stride + ichunk + 8);
	  w2 = _mm256_load_ps(weights + ifreq * stride + ichunk + 16);
	  w3 = _mm256_load_ps(weights + ifreq * stride + ichunk + 24);

	  cout << "Load the things" << endl;
 
	  // First, we need to see how many of the weights are equal to zero
	  __m256i nfail = _mm256_add_epi32(_mm256_add_epi32(_mm256_add_epi32(_mm256_castps_si256(_mm256_cmp_ps(w0, zero, _CMP_GT_OS)), 
									     _mm256_castps_si256(_mm256_cmp_ps(w1, zero, _CMP_GT_OS))), 
							    _mm256_castps_si256(_mm256_cmp_ps(w2, zero, _CMP_GT_OS))), 
					   _mm256_castps_si256(_mm256_cmp_ps(w3, zero, _CMP_GT_OS)));

	  // If nfail is greater than -8, we treat this as a failed v1 case
	  // To do this, convert nfail to an __m256 and get a constant register with a horizontal sum
	  // Then, make a mask that is all 1 is the constant register is less than -8 (successful v1) or 0 if not (unsuccessful v1)
	  __m256 mask = _mm256_cmp_ps(_mm256_cvtepi32_ps(nfail), _mm256_set1_ps(-8), _CMP_LT_OS);

	  cout << "Calculate the mask" << endl;

	  // Here, we do the variance computation:
	  vsum = _mm256_fmadd_ps(_mm256_mul_ps(i3, i3), w3, _mm256_fmadd_ps(_mm256_mul_ps(i2, i2), w2, _mm256_fmadd_ps(_mm256_mul_ps(i1, i1), w1, _mm256_mul_ps(_mm256_mul_ps(i0, i0), w0))));
	  wsum = hadd(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(w0, w1), w2), w3));
	  wsum = _mm256_min_ps(wsum, _mm256_set1_ps(-0.0001));
	  tmp_var = hadd(vsum / wsum);

	  cout << "Get the variance estimate" << endl;
	  // Then, use update rules to update value we'll eventually set as our running variance (prevent it from changing too much)
	  tmp1 = _mm256_fmadd_ps(vw, tmp_var, _mm256_mul_ps(_mm256_set1_ps(1 - var_weight), tmp_var));
	  tmp2 = _mm256_add_ps(_mm256_add_ps(tmp_var, _mm256_set1_ps(var_clamp_add)), _mm256_mul_ps(tmp_var, _mm256_set1_ps(var_clamp_mult)));
	  tmp_var = _mm256_min_ps(tmp1, tmp2);
	  tmp_var = _mm256_max_ps(tmp_var, _mm256_set1_ps(running_var[ifreq] - var_clamp_add - running_var[ifreq] * var_clamp_mult));
	  
	  cout << "Do the variance update" << endl;

	  // Finally, mask fill with the running variance -- if weights less than cutoff AND the mask says our variance estimate passed, fill
	  res0 = _mm256_blendv_ps(i0, _mm256_mul_ps(rn.gen_floats(), _mm256_mul_ps(root_three, _mm256_sqrt_ps(tmp_var))), _mm256_and_ps(_mm256_cmp_ps(w0, c, _CMP_LT_OS), mask));
	  res1 = _mm256_blendv_ps(i1, _mm256_mul_ps(rn.gen_floats(), _mm256_mul_ps(root_three, _mm256_sqrt_ps(tmp_var))), _mm256_and_ps(_mm256_cmp_ps(w1, c, _CMP_LT_OS), mask));
	  res2 = _mm256_blendv_ps(i2, _mm256_mul_ps(rn.gen_floats(), _mm256_mul_ps(root_three, _mm256_sqrt_ps(tmp_var))), _mm256_and_ps(_mm256_cmp_ps(w2, c, _CMP_LT_OS), mask));
	  res3 = _mm256_blendv_ps(i3, _mm256_mul_ps(rn.gen_floats(), _mm256_mul_ps(root_three, _mm256_sqrt_ps(tmp_var))), _mm256_and_ps(_mm256_cmp_ps(w3, c, _CMP_LT_OS), mask));
	  
	  cout << "Blend the variance" << endl;

	  // We also need to modify the weight values
	  __m256 w = _mm256_blendv_ps(_mm256_set1_ps(w_clamp), _mm256_set1_ps(-w_clamp), mask);    // either +w_clamp or -w_clamp
	  tmp_w = _mm256_add_ps(w, tmp_w);
	  tmp_w = _mm256_max_ps(tmp_w, zero);
	  tmp_w = _mm256_min_ps(tmp_w, two);

	  cout << "Blend the weights" << endl;

	  // Store the new intensity values
	  _mm256_storeu_ps((float*) (intensity + ifreq * stride + ichunk), res0);
	  _mm256_storeu_ps((float*) (intensity + ifreq * stride + ichunk + 8), res1);
	  _mm256_storeu_ps((float*) (intensity + ifreq * stride + ichunk + 16), res2);
	  _mm256_storeu_ps((float*) (intensity + ifreq * stride + ichunk + 24), res3);
	  
	  cout << "Store the intensities" << endl;

	  // Store the new weight values
	  _mm256_storeu_ps((float*) (weights + ifreq * stride + ichunk), tmp_w);
	  _mm256_storeu_ps((float*) (weights + ifreq * stride + ichunk + 8), tmp_w);
	  _mm256_storeu_ps((float*) (weights + ifreq * stride + ichunk + 16), tmp_w);
	  _mm256_storeu_ps((float*) (weights + ifreq * stride + ichunk + 24), tmp_w);

	  cout << "Store the weights" << endl;
      }

      // Update the (scalar) running variance -- thanks Kendrick!
      // Convert the constant register 'x' to a scalar, and write to q[0].
      // First step: extract elements 0-3 into a 128-bit register.
      __m128 y = _mm256_extractf128_ps(tmp_var, 0);
      __m128 z = _mm256_extractf128_ps(tmp_w, 0);
      // The second step is really strange!  The intrinsic _mm_extract_ps()
      // extracts element 0 from the 128-bit register, but it has the wrong
      // return type (int32 instead of float32).  The returned value is a "fake"
      // int32 obtained by interpreting the bit pattern of the "real" float32
      // (in IEEE-754 representation) as an int32 (in twos-complement representation),
      // so it's not very useful.  Nevertheless if we just write it to memory as an int32,
      // and read back from the same memory location later as a float32, we'll get the
      // right answer.  This doesn't make much sense and there must be a better way to do it!
      int *i = reinterpret_cast<int *> (&running_var[ifreq]);   // hack: (int *) pointing to same memory location as q[0]
      *i = _mm_extract_ps(y, 0); // write a "fake" int32 to this memory location
      int *j = reinterpret_cast<int *> (&running_weights[ifreq]);   // hack: (int *) pointing to same memory location as q[0]
      *j = _mm_extract_ps(z, 0); // write a "fake" int32 to this memory location
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


}  // namespace rf_pipelines
