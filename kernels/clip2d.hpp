#ifndef _RF_PIPELINES_KERNELS_CLIP2D_HPP
#define _RF_PIPELINES_KERNELS_CLIP2D_HPP

#include <cassert>  // XXX remove
#include "mean_rms_accumulator.hpp"
#include "downsample.hpp"

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif

template<typename T, unsigned int S> using simd_t = simd_helpers::simd_t<T,S>;


// -------------------------------------------------------------------------------------------------
//
// The "bottom line" routine defined in this file is:
//
//    void _kernel_clip2d_wrms<T,S,Df,Dt> (simd_t<T,S> &mean, simd_t<T,S> &rms, const T *intensity, 
//                                         const T *weights, int nfreq, int nt, int stride);
//
// FIXME should add R argument


template<typename T, unsigned int S, unsigned int Df, unsigned int Dt>
void _kernel_clip2d_wrms(simd_t<T,S> &mean, simd_t<T,S> &rms, const T *intensity, const T *weights, int nfreq, int nt, int stride)
{
    // XXX assert -> throw
    assert(nfreq > 0);
    assert(nt > 0);
    assert(nfreq % Df == 0);
    assert(nt % (Dt*S) == 0);

    mean_rms_accumulator<T,S> acc;

    for (int ifreq = 0; ifreq < nfreq; ifreq += Df) {
	const T *irow = intensity + ifreq*stride;
	const T *wrow = weights + ifreq*stride;

	for (int it = 0; it < nt; it += Dt*S) {
	    simd_t<T,S> ival = _kernel_downsample<T,S,Df,Dt> (irow + it, stride);
	    simd_t<T,S> wval = _kernel_downsample<T,S,Df,Dt> (wrow + it, stride);
	    acc.accumulate(ival, wval);
	}
    }

    acc.horizontal_sum();
    acc.get_mean_rms(mean, rms);
}


}  // namespace rf_pipelines

#endif
