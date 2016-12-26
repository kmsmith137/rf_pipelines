#ifndef _RF_PIPELINES_KERNELS_CLIP2D_HPP
#define _RF_PIPELINES_KERNELS_CLIP2D_HPP

#include <cassert>  // XXX remove
#include "mean_rms_accumulator.hpp"
#include "downsample.hpp"
#include "mask.hpp"

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif

template<typename T, unsigned int S> using simd_t = simd_helpers::simd_t<T,S>;


// -------------------------------------------------------------------------------------------------
//
// _kernel_clip2d_wrms<T,S,Df,Dt>(): computes the weighted mean and rms of a 2D strided array,
// with downsampling factors (Df,Dt) in the (frequency,time) axes.
//
// The 'mean' and 'rms' outputs are simd vectors whose elements are all equal.
// If the weighted mean and rms cannot be computed (e.g. because all weights are zero), then
// rms=0 and mean is arbitrary.  (This behavior is inherited from 'struct ean_rms_accumulator'.)
//
// There are three variants of this routine, depending on whether the ds_intensity and ds_weights
// outputs are written.  (These are downsampled copies of the intensity and weights arrays.)
//
//    void _kernel_clip2d_wrms<T,S,Df,Dt> (simd_t<T,S> &mean, simd_t<T,S> &rms, const T *intensity, 
//                                         const T *weights, int nfreq, int nt, int stride);
//
//    void _kernel_clip2d_wrms<T,S,Df,Dt> (simd_t<T,S> &mean, simd_t<T,S> &rms, const T *intensity, 
//                                         const T *weights, int nfreq, int nt, int stride,
//                                         T *ds_intensity);
//
//    void _kernel_clip2d_wrms<T,S,Df,Dt> (simd_t<T,S> &mean, simd_t<T,S> &rms, const T *intensity, 
//                                         const T *weights, int nfreq, int nt, int stride,
//                                         T *ds_intensity, T *ds_weights);
//
// FIXME there is a lot of cut-and-pasted code between these three routines.  This could
// probably be fixed with some C++ template-ology.
//
// FIXME a more general version of these kernels is possible, with a template parameter R
// which controlls the number of rows read in each pass.


template<typename T, unsigned int S, unsigned int Df, unsigned int Dt>
inline void _kernel_clip2d_wrms(simd_t<T,S> &mean, simd_t<T,S> &rms, const T *intensity, const T *weights, int nfreq, int nt, int stride)
{
    // XXX assert -> throw
    assert(nfreq > 0);
    assert(nt > 0);
    assert(nfreq % Df == 0);
    assert(nt % (Dt*S) == 0);

    mean_rms_accumulator<T,S> acc;

    const simd_t<T,S> zero = simd_t<T,S>::zero();
    const simd_t<T,S> one = simd_t<T,S> (1.0);

    for (int ifreq = 0; ifreq < nfreq; ifreq += Df) {
	const T *irow = intensity + ifreq*stride;
	const T *wrow = weights + ifreq*stride;

	for (int it = 0; it < nt; it += Dt*S) {
	    simd_t<T,S> wival, wval;
	    _kernel_downsample<T,S,Df,Dt> (wival, wval, irow + it, wrow + it, stride);

	    simd_t<T,S> ival = wival / blendv(wval.compare_gt(zero), wval, one);
	    acc.accumulate(ival, wval);
	}
    }

    acc.horizontal_sum();
    acc.get_mean_rms(mean, rms);
}


template<typename T, unsigned int S, unsigned int Df, unsigned int Dt>
inline void _kernel_clip2d_wrms(simd_t<T,S> &mean, simd_t<T,S> &rms, const T *intensity, const T *weights, int nfreq, int nt, int stride, T *ds_intensity)
{
    // XXX assert -> throw
    assert(nfreq > 0);
    assert(nt > 0);
    assert(nfreq % Df == 0);
    assert(nt % (Dt*S) == 0);

    mean_rms_accumulator<T,S> acc;

    const simd_t<T,S> zero = simd_t<T,S>::zero();
    const simd_t<T,S> one = simd_t<T,S> (1.0);

    for (int ifreq = 0; ifreq < nfreq; ifreq += Df) {
	const T *irow = intensity + ifreq*stride;
	const T *wrow = weights + ifreq*stride;

	for (int it = 0; it < nt; it += Dt*S) {
	    simd_t<T,S> wival, wval;
	    _kernel_downsample<T,S,Df,Dt> (wival, wval, irow + it, wrow + it, stride);

	    simd_t<T,S> ival = wival / blendv(wval.compare_gt(zero), wval, one);
	    acc.accumulate(ival, wval);

	    ival.storeu(ds_intensity);
	    ds_intensity += S;
	}
    }

    acc.horizontal_sum();
    acc.get_mean_rms(mean, rms);
}

template<typename T, unsigned int S, unsigned int Df, unsigned int Dt>
inline void _kernel_clip2d_wrms(simd_t<T,S> &mean, simd_t<T,S> &rms, const T *intensity, const T *weights, int nfreq, int nt, int stride, T *ds_intensity, T *ds_weights)
{
    // XXX assert -> throw
    assert(nfreq > 0);
    assert(nt > 0);
    assert(nfreq % Df == 0);
    assert(nt % (Dt*S) == 0);

    mean_rms_accumulator<T,S> acc;

    const simd_t<T,S> zero = simd_t<T,S>::zero();
    const simd_t<T,S> one = simd_t<T,S> (1.0);

    for (int ifreq = 0; ifreq < nfreq; ifreq += Df) {
	const T *irow = intensity + ifreq*stride;
	const T *wrow = weights + ifreq*stride;

	for (int it = 0; it < nt; it += Dt*S) {
	    simd_t<T,S> wival, wval;
	    _kernel_downsample<T,S,Df,Dt> (wival, wval, irow + it, wrow + it, stride);

	    simd_t<T,S> ival = wival / blendv(wval.compare_gt(zero), wval, one);
	    acc.accumulate(ival, wval);

	    ival.storeu(ds_intensity);
	    ds_intensity += S;

	    wval.storeu(ds_weights);
	    ds_weights += S;
	}
    }

    acc.horizontal_sum();
    acc.get_mean_rms(mean, rms);
}


// -------------------------------------------------------------------------------------------------


template<typename T, unsigned int S, unsigned int Df, unsigned int Dt>
inline void _kernel_clip2d_mask(T *weights, const T *ds_intensity, simd_t<T,S> mean, simd_t<T,S> thresh, int nfreq, int nt, int stride, int ds_stride)
{
    const T *ds_irow = ds_intensity;

    // XXX assert -> throw
    assert(nfreq > 0);
    assert(nt > 0);
    assert(nfreq % Df == 0);
    assert(nt % (Dt*S) == 0);

    for (int ifreq = 0; ifreq < nfreq; ifreq += Df) {
	const T *ds_itmp = ds_irow;
	T *wrow = weights + ifreq * stride;

	for (int it = 0; it < nt; it += Dt*S) {
	    simd_t<T,S> ival = simd_t<T,S>::loadu(ds_itmp);
	    ds_itmp += S;

	    ival -= mean;
	    ival = ival.abs();

	    simd_t<int,S> valid = ival.compare_lt(thresh);
	    _kernel_mask<T,S,Df,Dt> (wrow + it, valid, stride);
	}

	ds_irow += ds_stride;
    }
}


}  // namespace rf_pipelines

#endif
