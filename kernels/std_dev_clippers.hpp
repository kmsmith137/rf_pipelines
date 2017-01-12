#ifndef _RF_PIPELINES_KERNELS_STD_DEV_CLIPPERS_HPP
#define _RF_PIPELINES_KERNELS_STD_DEV_CLIPPERS_HPP

#include "downsample.hpp"
#include "mean_rms_accumulator.hpp"

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif


// -------------------------------------------------------------------------------------------------
//
// _kernel_std_dev_t()
//
// The output arrays 'out_sd' and 'out_valid' are 1D arrays of shape (nfreq/Df).
//
// Caller must check (nfreq % Df) == 0 and (nt % (Dt*S)) == 0.
// There is no requirement that (nfreq % (Df*S)) == 0.


template<typename T, unsigned int S, unsigned int Df, unsigned int Dt>
inline void _kernel_std_dev_t(T *out_sd, smask_t<T,1> *out_valid, const T *intensity, const T *weights, int nfreq, int nt, int stride)
{
    const simd_t<T,S> zero = simd_t<T,S>::zero();
    const simd_t<T,S> one = simd_t<T,S> (1.0);

    for (int ifreq = 0; ifreq < nfreq; ifreq += Df) {
	const T *irow = intensity + ifreq * stride;
	const T *wrow = weights + ifreq * stride;

	mean_rms_accumulator<T,S> acc;
	simd_t<T,S> mean, var;

	for (int it = 0; it < nt; it += Dt*S) {
	    simd_t<T,S> wival, wval;
	    _kernel_downsample<T,S,Df,Dt> (wival, wval, irow + it, wrow + it, stride);

	    simd_t<T,S> ival = wival / blendv(wval.compare_gt(zero), wval, one);
	    acc.accumulate(ival, wval);	    
	}

	acc.horizontal_sum();
	acc.get_mean_variance(mean, var);

	// scalar instructions should be fine here
	T sd = var.template extract<0> ();
	*out_sd++ = sd;
	*out_valid++ = (sd > 0.0) ? smask_t<T,1>(-1) : 0;
    }
}


}  // namespace rf_pipelines

#endif
