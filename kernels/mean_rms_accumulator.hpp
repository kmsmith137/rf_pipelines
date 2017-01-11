#ifndef _RF_PIPELINES_KERNELS_MEAN_RMS_ACCUMULATOR_HPP
#define _RF_PIPELINES_KERNELS_MEAN_RMS_ACCUMULATOR_HPP

#include <simd_helpers/convert.hpp>

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif


template<typename T, unsigned int S> using simd_t = simd_helpers::simd_t<T,S>;
template<typename T, unsigned int S> using smask_t = simd_helpers::smask_t<T,S>;
template<typename T, unsigned int S, unsigned int N> using simd_ntuple = simd_helpers::simd_ntuple<T,S,N>;


// mean_rms_accumulator: helper class for computing weighted mean/rms (e.g. in clipper_transform)
//
// Each element of the simd_t<T,S> is processed independently.  If the desired behavior is to
// sum them together instead, then call mean_rms_accumulator::horizontal_sum() before calling
// mean_rms_accumulator::get_mean_rms() but after the calls to accumulate().
//
// An entry is invalid if the mean and rms cannot be computed, either because the
// sum of the weights is <= 0, or if the variance is too small compared to the mean.
// Invalid entries are indicated by rms=0 and arbitrary mean.
//
// The type simd_t<T,S> is the type which is used internally when accumulating samples.
// This need not be the same as the type simd_t<Td,Sd> of the data samples themselves.
// E.g. in the clipper_transform we use T=double and Td=float by default.
// The only requirement is that the converter simd_t<Td,Sd> -> simd_ntuple<T,S,N>,
// where N = Sd/S, is defined in simd_helpers/convert.hpp.

template<typename T, unsigned int S>
struct mean_rms_accumulator {
    simd_t<T,S> acc0;    // sum_i W_i
    simd_t<T,S> acc1;    // sum_i W_i I_i
    simd_t<T,S> acc2;    // sum_i W_i I_i^2

    mean_rms_accumulator()
    {
	acc0 = simd_t<T,S>::zero();
	acc1 = simd_t<T,S>::zero();
	acc2 = simd_t<T,S>::zero();
    }

    inline void accumulate(simd_t<T,S> ival, simd_t<T,S> wval)
    {
	simd_t<T,S> wi = wval * ival;
	acc0 += wval;
	acc1 += wi;
	acc2 += wi * ival;
    }

    template<unsigned int N, typename std::enable_if<(N>0),int>::type = 0>
    inline void accumulate(simd_ntuple<T,S,N> ival, simd_ntuple<T,S,N> wval)
    {
	accumulate(ival.v, wval.v);
	accumulate(ival.x, wval.x);
    }

    inline void accumulate(simd_ntuple<T,S,0> ival, simd_ntuple<T,S,0> wval)
    {
	return;
    }

    template<typename Td, unsigned int Sd>
    inline void accumulate(simd_t<Td,Sd> ival, simd_t<Td,Sd> wval)
    {
	static_assert(Sd % S == 0, "mean_rms_accumulator: \"data\" simd size Sd must be a multiple of \"accumulator\" simd size S");

	static constexpr unsigned int N = Sd / S;
	simd_ntuple<T,S,N> ival_c, wval_c;

	convert(ival_c, ival);
	convert(wval_c, wval);
	accumulate(ival_c, wval_c);
    }

    inline void horizontal_sum()
    {
	acc0 = acc0.horizontal_sum();
	acc1 = acc1.horizontal_sum();
	acc2 = acc2.horizontal_sum();
    }

    inline void get_mean_variance(simd_t<T,S> &mean, simd_t<T,S> &var) const
    {
	static constexpr T eps = 1.0e3 * simd_helpers::machine_epsilon<T> ();

	smask_t<T,S> valid = acc0.compare_gt(simd_t<T,S>::zero());

	simd_t<T,S> t0 = blendv(valid, acc0, simd_t<T,S>(1.0));
	mean = acc1 / t0;

	simd_t<T,S> mean2 = mean * mean;
	var = acc2/t0 - mean2;

	simd_t<T,S> thresh = simd_t<T,S>(eps) * mean2;
	valid = valid.bitwise_and(var.compare_gt(thresh));
	var = var.apply_mask(valid);
    }

    inline void get_mean_rms(simd_t<T,S> &mean, simd_t<T,S> &rms) const
    {
	simd_t<T,S> var;
	get_mean_variance(mean, var);
	rms = var.sqrt();
    }
};


}  // namespace rf_pipelines

#endif
