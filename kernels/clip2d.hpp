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
// _kernel_clip2d_wrms<T,S,Df,Dt,Iflag,Wflag>(simd_t<T,S> &mean, simd_t<T,S> &rms, const T *intensity, 
//                                            const T *weights, int nfreq, int nt, int stride, 
//                                            T *ds_intensity, T *ds_weights)
//
// Computes the weighted mean and rms of a 2D strided array,
// with downsampling factors (Df,Dt) in the (frequency,time) axes.
//
// The 'mean' and 'rms' outputs are simd vectors whose elements are all equal.
// If the weighted mean and rms cannot be computed (e.g. because all weights are zero), then
// rms=0 and mean is arbitrary.  (This behavior is inherited from 'struct mean_rms_accumulator'.)
//
// As the downsampled intensity and weights arrays are computed, they are written to
// 'ds_intensity' and 'ds_weights'.  These are unstrided arrays, i.e. the row stride
// is (nt/Dt).
//
// The Iflag, Wflag template arguments will omit writing the ds_intensity, ds_weights
// arrays if set to 'false'.  In this case, passing a NULL pointer is OK.
//
// FIXME a more general version of these kernels is possible, with a template parameter R
// which controlls the number of rows read in each pass.


// _clip2d_write_if(): helper for _kernel_clip2d_wrms()
template<typename T, unsigned int S, bool flag, typename std::enable_if<flag,int>::type = 0>
inline void _clip2d_write_if(T*& p, simd_t<T,S> x)
{
    x.storeu(p);
    p += S;
}

template<typename T, unsigned int S, bool flag, typename std::enable_if<(!flag),int>::type = 0>
inline void _clip2d_write_if(T*& p, simd_t<T,S> x)
{
    return;
}


template<typename T, unsigned int S, unsigned int Df, unsigned int Dt, bool Iflag, bool Wflag>
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

	    _clip2d_write_if<T,S,Iflag> (ds_intensity, ival);
	    _clip2d_write_if<T,S,Wflag> (ds_weights, wval);
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


// -------------------------------------------------------------------------------------------------


template<typename T, unsigned int S>
inline void _kernel_clip2d_iterate(simd_t<T,S> &out_mean, simd_t<T,S> &out_rms, const T *intensity, const T *weights,
				   simd_t<T,S> in_mean, simd_t<T,S> in_thresh, int nfreq, int nt, int stride)
{
    mean_rms_accumulator<T,S> acc;

    for (int ifreq = 0; ifreq < nfreq; ifreq++) {
	const T *irow = intensity + ifreq*stride;
	const T *wrow = weights + ifreq*stride;

	for (int it = 0; it < nt; it += S) {
	    simd_t<T,S> ival = simd_t<T,S>::loadu(irow + it);
	    simd_t<T,S> wval = simd_t<T,S>::loadu(wrow + it);

	    ival -= in_mean;
	    ival = ival.abs();

	    simd_t<int,S> valid = ival.compare_lt(in_thresh);
	    wval = wval.bitwise_and(valid);

	    acc.accumulate(ival, wval);
	}
    }

    acc.horizontal_sum();
    acc.get_mean_rms(out_mean, out_rms);
}


}  // namespace rf_pipelines

#endif
