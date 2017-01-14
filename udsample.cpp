// Fast downsampling kernels.  Maybe I'll write fast upsampling kernels some day too!
//
// FIXME: currently we need to compile a new kernel for every (Df,Dt) pair, where
// Df,Dt are the frequency/time downsampling factors.  Eventually I'd like to 
// improve this by having special kernels to handle the large-Df and large-Dt cases.

#include <array>
#include <cassert>

#include "rf_pipelines_internals.hpp"

#include "kernels/downsample.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif


// Note: kernels/downsample.hpp defines
//   _kernel_downsample_2d<T,S,Df,Dt> (out_intensity, out_weights, out_stride, in_intensity, in_nfreq, in_nt, in_stride)

using kernel_downsample_t = void (*) (float *, float *, int, const float *, const float *, int, int, int);


// -------------------------------------------------------------------------------------------------
//
// fill_2d_kernel_table<S,MaxDf,MaxDt>()


template<unsigned int S, unsigned int Df, unsigned int MaxDt, typename std::enable_if<(MaxDt==0),int>::type = 0>
inline void fill_1d_kernel_table(kernel_downsample_t *out) { }

template<unsigned int S, unsigned int Df, unsigned int MaxDt, typename std::enable_if<(MaxDt>0),int>::type = 0>
inline void fill_1d_kernel_table(kernel_downsample_t *out)
{
    constexpr int NDt = IntegerLog2<MaxDt>() + 1;

    fill_1d_kernel_table<S,Df,(MaxDt/2)> (out);
    out[NDt-1] = _kernel_downsample_2d<float,S,Df,MaxDt>;
}


template<unsigned int S, unsigned int MaxDf, unsigned int MaxDt, typename std::enable_if<(MaxDf==0),int>::type = 0>
inline void fill_2d_kernel_table(kernel_downsample_t *out) { }

template<unsigned int S, unsigned int MaxDf, unsigned int MaxDt, typename std::enable_if<(MaxDf>0),int>::type = 0>
inline void fill_2d_kernel_table(kernel_downsample_t *out) 
{ 
    constexpr int NDf = IntegerLog2<MaxDf>() + 1;
    constexpr int NDt = IntegerLog2<MaxDt>() + 1;

    fill_2d_kernel_table<S,(MaxDf/2),MaxDt> (out);
    fill_1d_kernel_table<S,MaxDf,MaxDt> (out + (NDf-1)*NDt);
}


// -------------------------------------------------------------------------------------------------


struct kernel_table {
    static constexpr int S = constants::single_precision_simd_length;
    static constexpr int MaxDf = constants::max_frequency_downsampling;
    static constexpr int MaxDt = constants::max_time_downsampling;
    static constexpr int NDf = IntegerLog2<MaxDf>() + 1;
    static constexpr int NDt = IntegerLog2<MaxDt>() + 1;

    // Weird: for some reason using std::max() here gives a clang linker (not compiler) error.
    static constexpr int MaxD = (MaxDf > MaxDt) ? MaxDf : MaxDt;

    std::vector<kernel_downsample_t> entries;

    integer_log2_lookup_table ilog2_lookup;


    kernel_table() : entries(NDf*NDt), ilog2_lookup(MaxD)
    {
	fill_2d_kernel_table<S,MaxDf,MaxDt> (&entries[0]);
    }

    // Before calling get_kernel(), caller must check:
    //    1 <= Df <= MaxDf 
    //    1 <= Dt <= MaxDt
    //    both Df,Dt are powers of two

    inline kernel_downsample_t get_kernel(int Df, int Dt)
    {
	int idf = ilog2_lookup(Df);
	int idt = ilog2_lookup(Dt);

	return entries[idf*NDt + idt];
    }
};


// -------------------------------------------------------------------------------------------------


static kernel_table global_kernel_table;


// Externally visible downsampling routine declared in rf_pipelines.hpp
void wi_downsample(float *out_intensity, float *out_weights, int out_stride, const float *in_intensity, const float *in_weights, int in_nfreq, int in_nt, int in_stride, int Df, int Dt)
{
    static constexpr int S = constants::single_precision_simd_length;
    static constexpr int MaxDf = constants::max_frequency_downsampling;
    static constexpr int MaxDt = constants::max_time_downsampling;

    // Part 1: Argument megacheck

    if (_unlikely((in_nfreq <= 0) || (in_nt <= 0)))
	throw runtime_error("wi_downsample(): (in_nfreq,in_nt)=(" + to_string(in_nfreq) + "," + to_string(in_nt) + ") is invalid");

    if (_unlikely((Df <= 0) || (Dt <= 0)))
	throw runtime_error("wi_downsample(): (Df,Dt)=(" + to_string(Df) + "," + to_string(Dt) + ") is invalid");

    if (_unlikely(!is_power_of_two(Df) || !is_power_of_two(Dt)))
	throw runtime_error("wi_downsample(): (Df,Dt)=(" + to_string(Df) + "," + to_string(Dt) + ") must be powers of two");

    if (_unlikely(in_nfreq % Df))
	throw runtime_error("wi_downsample(): in_nfreq=" + to_string(in_nfreq) + " is not divisible by Df=" + to_string(Df));

    if (_unlikely(in_nt % (Dt*S)))
	throw runtime_error("wi_downsample(): in_nt=" + to_string(in_nfreq) + " is not divisible by Dt*S, where Dt=" + to_string(Dt)
			    + " is the time downsampling factor and S=" + to_string(S) + " is the single-precision simd length on this machine");

    if (_unlikely(out_stride < in_nt/Dt))
	throw runtime_error("wi_downsample(): out_stride=" + to_string(out_stride) + " is < (in_nt/Dt) =" + to_string(in_nt/Dt));
    
    if (_unlikely(in_stride < in_nt))
	throw runtime_error("wi_downsample(): in_stride=" + to_string(in_stride) + " is < in_nt=" + to_string(in_nt));	

    if (_unlikely((Df > MaxDf) || (Dt > MaxDt))) {
	throw runtime_error("wi_downsample(): (Df,Dt)=(" + to_string(Df) + "," + to_string(Dt) + ")"
			    + " exceeds compile time limits; to fix this see 'constants' in rf_pipelines.hpp");
    }

    if (_unlikely(!out_intensity))
	throw runtime_error("wi_downsample(): 'out_intensity' argument is a NULL pointer");

    if (_unlikely(!out_weights))
	throw runtime_error("wi_downsample(): 'out_weights' argument is a NULL pointer");

    if (_unlikely(!in_intensity))
	throw runtime_error("wi_downsample(): 'in_intensity' argument is a NULL pointer");

    if (_unlikely(!in_weights))
	throw runtime_error("wi_downsample(): 'in_weights' argument is a NULL pointer");

    // Part 2: Get and apply downsampling kernel

    auto f = global_kernel_table.get_kernel(Df, Dt);

    f(out_intensity, out_weights, out_stride, in_intensity, in_weights, in_nfreq, in_nt, in_stride);
}


}  // namespace rf_pipelines
