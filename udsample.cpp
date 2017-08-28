// Fast downsampling kernels.  Maybe I'll write fast upsampling kernels some day too!
//
// FIXME: currently we need to compile a new kernel for every (Df,Dt) pair, where
// Df,Dt are the frequency/time downsampling factors.  Eventually I'd like to 
// improve this by having special kernels to handle the large-Df and large-Dt cases.

#include <cassert>
#include "rf_kernels/downsample.hpp"
#include "rf_pipelines_internals.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif


// Externally visible downsampling routine declared in rf_pipelines.hpp
//
// Note that the normalization of the downsampled weights array differs (by a factor of Df*Dt) 
// from the python version of wi_downsample().  This is actually nontrivial to change, since
// the unit tests currently assume strict equivalence with the downsampling logic in the
// intensity_clipper, which we don't want to slow down by including an extra multiplication.

void wi_downsample(float *out_intensity, float *out_weights, int out_stride, const float *in_intensity, const float *in_weights, int in_nfreq, int in_nt, int in_stride, int Df, int Dt)
{
    if (_unlikely((in_nfreq <= 0) || (in_nt <= 0)))
	throw runtime_error("wi_downsample(): (in_nfreq,in_nt)=(" + to_string(in_nfreq) + "," + to_string(in_nt) + ") is invalid");

    if (_unlikely((Df <= 0) || (Dt <= 0)))
	throw runtime_error("wi_downsample(): (Df,Dt)=(" + to_string(Df) + "," + to_string(Dt) + ") is invalid");

    if (_unlikely(in_nfreq % Df))
	throw runtime_error("wi_downsample(): in_nfreq=" + to_string(in_nfreq) + " is not divisible by Df=" + to_string(Df));

    if (_unlikely(in_nt % Dt))
	throw runtime_error("wi_downsample(): in_nt=" + to_string(in_nfreq) + " is not divisible by Dt=" + to_string(Dt));

    int nfreq_out = in_nfreq / Df;
    int nt_out = in_nt / Dt;
    
    rf_kernels::wi_downsampler ds(Df, Dt);
    ds.downsample(nfreq_out, nt_out, out_intensity, out_weights, out_stride, in_intensity, in_weights, in_stride);
}


}  // namespace rf_pipelines
