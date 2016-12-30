#ifndef _RF_PIPELINES_KERNELS_MEAN_RMS_ACCUMULATOR_HPP
#define _RF_PIPELINES_KERNELS_MEAN_RMS_ACCUMULATOR_HPP

#include <simd_helpers/simd_t.hpp>

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif


template<typename T, unsigned int S> using simd_t = simd_helpers::simd_t<T,S>;


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
// FIXME I think we should compute the mean/rms in double precision even when T=float.
// Suggested generalization: struct mean_rms_accumulator<float, 8, double>.


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

    inline void horizontal_sum()
    {
	acc0 = acc0.horizontal_sum();
	acc1 = acc1.horizontal_sum();
	acc2 = acc2.horizontal_sum();
    }

    inline void get_mean_rms(simd_t<T,S> &mean, simd_t<T,S> &rms) const
    {
	static constexpr T eps = 1.0e3 * simd_helpers::machine_epsilon<T> ();

	simd_t<int,S> valid = acc0.compare_gt(simd_t<T,S>::zero());

	simd_t<T,S> t0 = blendv(valid, acc0, simd_t<T,S>(1.0));
	mean = acc1 / t0;

	simd_t<T,S> mean2 = mean * mean;
	simd_t<T,S> var = acc2/t0 - mean2;

	simd_t<T,S> thresh = simd_t<T,S>(eps) * mean2;
	valid = valid.bitwise_and(var.compare_gt(thresh));
	var = var.bitwise_and(valid);
	rms = var.sqrt();
    }

    inline void get_mean_rms(simd_t<T,S> &mean, simd_t<T,S> &rms, simd_t<T,S> sigma) const
    {
	get_mean_rms(mean, rms);
	rms *= sigma;
    }
};


}  // namespace rf_pipelines

#endif
