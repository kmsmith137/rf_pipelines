// Kernels defined here:
//
// _kernel_noniterative_wrms_2d(): computes weighted mean/rms of a 2D array with optional downsampling 
//    (there is also an option to write the downsampled intensity/weights to auxiliary arrays)
//
// _kernel_clip2d_iterate(): computes weighted mean/rms of a 2D array, including only elements
//    in a certain range.  (there is no downsampling option here)
//
// _kernel_clip2d_mask(): sets weights to zero when intensity is outside a certain range.
//    Optionally, the intensity array can be downsampled relative to the weights.
//
// _kernel_iterative_wrms_2d():
//   This is the "bottom line" routine which is wrapped by weighted_mean_and_rms().
//
// _kernel_clip_2d(): 
//   This is the "bottom line" kernel which is wrapped by the intensity_clipper(AXIS_NONE).


#ifndef _RF_PIPELINES_KERNELS_INTENSITY_CLIPPERS_HPP
#define _RF_PIPELINES_KERNELS_INTENSITY_CLIPPERS_HPP

#include <cassert>  // XXX remove

#include "mean_rms_accumulator.hpp"
#include "downsample.hpp"
#include "mask.hpp"

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif

template<typename T, unsigned int S> using simd_t = simd_helpers::simd_t<T,S>;

// _extract_first<T,S,S2>(x): a general-purpose helper which extracts first "S2" elements from x.
template<typename T, unsigned int S, unsigned int S2, typename std::enable_if<(S==S2),int>::type = 0>
inline simd_t<T,S2> _extract_first(simd_t<T,S> x) { return x; }

template<typename T, unsigned int S, unsigned int S2, typename std::enable_if<(S==2*S2),int>::type = 0>
inline simd_t<T,S2> _extract_first(simd_t<T,S> x) { return x.template extract_half<0>(); }
    

// _write_and_advance_if(): helper which writes data and advances a pointer, 
// but only if a specified boolean predicate evaluates to true at compile-time.
template<typename T, unsigned int S, bool flag, typename std::enable_if<flag,int>::type = 0>
inline void _write_and_advance_if(T*& p, simd_t<T,S> x)
{
    x.storeu(p);
    p += S;
}

template<typename T, unsigned int S, bool flag, typename std::enable_if<(!flag),int>::type = 0>
inline void _write_and_advance_if(T*& p, simd_t<T,S> x) { return; }


// -------------------------------------------------------------------------------------------------
//
// _kernel_noniterative_wrms_2d<T, S, Df, Dt, Iflag, Wflag, Ti, Si>
//    (simd_t<T,S> &mean, simd_t<T,S> &rms, const T *intensity, const T *weights, 
//     int nfreq, int nt, int stride, T *ds_intensity, T *ds_weights)
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
// The template parameters (Ti,Si) are the type and simd stride of the floating-point
// type used internally when accumulating samples and computing the mean/rms.  For example,
// it often makes sense to take (T,S)=(float,8) and (Ti,Si)=(double,4).
//
// FIXME a more general version of this kernel is possible, with a template parameter R
// which controlls the number of rows read in each pass.


template<typename T, unsigned int S, unsigned int Df, unsigned int Dt, bool Iflag, bool Wflag, typename Ti, unsigned int Si>
inline void _kernel_noniterative_wrms_2d(simd_t<T,S> &mean, simd_t<T,S> &rms, const T *intensity, const T *weights, int nfreq, int nt, int stride, T *ds_intensity, T *ds_weights)
{
    // XXX assert -> throw
    assert(nfreq > 0);
    assert(nt > 0);
    assert(nfreq % Df == 0);
    assert(nt % (Dt*S) == 0);

    mean_rms_accumulator<Ti,Si> acc;

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

	    _write_and_advance_if<T,S,Iflag> (ds_intensity, ival);
	    _write_and_advance_if<T,S,Wflag> (ds_weights, wval);
	}
    }

    acc.horizontal_sum();

    simd_t<Ti,Si> mean_i, rms_i;
    acc.get_mean_rms(mean_i, rms_i);
    
    convert(mean, _extract_first<Ti,Si,S> (mean_i));
    convert(rms, _extract_first<Ti,Si,S> (rms_i));
}


// -------------------------------------------------------------------------------------------------
//
// _kernel_clip2d_mask(): Masks all intensity samples which differ from the mean by more than 
// 'thresh'.  The intensity array can be downsampled relative to the weights array.


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

	    smask_t<T,S> valid = ival.compare_lt(thresh);
	    _kernel_mask<T,S,Df,Dt> (wrow + it, valid, stride);
	}

	ds_irow += ds_stride;
    }
}


// -------------------------------------------------------------------------------------------------
//
// _kernel_clip2d_iterate(): Computes the weighted mean/rms of a 2D array, using only samples
// which differ by the mean by more than 'thresh'.  No downsampling in this one.


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

	    simd_t<T,S> ival_c = (ival - in_mean).abs();
	    smask_t<T,S> valid = ival_c.compare_lt(in_thresh);

	    wval = wval.apply_mask(valid);
	    acc.accumulate(ival, wval);
	}
    }

    acc.horizontal_sum();
    acc.get_mean_rms(out_mean, out_rms);
}


// -------------------------------------------------------------------------------------------------
//
// _kernel_iterative_wrms_2d():
//   This is the "bottom line" routine which is wrapped by weighted_mean_and_rms().
//
// _kernel_clip_2d(): 
//    This is the "bottom line" routine which is wrapped by intensity_clipper(AXIS_NONE).
//
// Note: caller must ensure that the IterFlag compile-time argument is 'true' iff (niter > 0).


template<typename T, unsigned int S, unsigned int Df, unsigned int Dt, bool IterFlag>
inline void _kernel_iterative_wrms_2d(simd_t<T,S> &mean, simd_t<T,S> &rms, const T *intensity, const T *weights, 
				      int nfreq, int nt, int stride, int niter, double sigma, T *ds_int, T *ds_wt)
{
    static constexpr bool DsiFlag = (Df > 1) || (Dt > 1);
    static constexpr bool DswFlag = IterFlag && ((Df > 1) || (Dt > 1));

    _kernel_noniterative_wrms_2d<T,S,Df,Dt,DsiFlag,DswFlag,T,S> (mean, rms, intensity, weights, nfreq, nt, stride, ds_int, ds_wt);

    const T *s_intensity = DsiFlag ? ds_int : intensity;
    const T *s_weights = DswFlag ? ds_wt : weights;
    int s_stride = DsiFlag ? (nt/Dt) : stride;   // must use DsiFlag here, not DswFlag
	
    for (int iter = 1; iter < niter; iter++) {
	simd_t<T,S> thresh = simd_t<T,S>(sigma) * rms;
	_kernel_clip2d_iterate<T,S> (mean, rms, s_intensity, s_weights, mean, thresh, nfreq/Df, nt/Dt, s_stride);
    }
}


template<typename T, unsigned int S, unsigned int Df, unsigned int Dt, bool IterFlag>
inline void _kernel_clip_2d(T *intensity, T *weights, int nfreq, int nt, int stride, int niter, double sigma, double iter_sigma, T *ds_int, T *ds_wt)
{
    static constexpr bool DsiFlag = (Df > 1) || (Dt > 1);

    simd_t<T,S> mean, rms;
    _kernel_iterative_wrms_2d<T,S,Df,Dt,IterFlag> (mean, rms, intensity, weights, nfreq, nt, stride, niter, iter_sigma, ds_int, ds_wt);

    const T *s_intensity = DsiFlag ? ds_int : intensity;
    int s_stride = DsiFlag ? (nt/Dt) : stride;   // must use DsiFlag here, not DswFlag

    // (s_intensity, weights, sigma) here
    simd_t<T,S> thresh = simd_t<T,S>(sigma) * rms;
    _kernel_clip2d_mask<T,S,Df,Dt> (weights, s_intensity, mean, thresh, nfreq, nt, stride, s_stride);
}


// -------------------------------------------------------------------------------------------------


template<typename T, unsigned int S, unsigned int Df, unsigned int Dt, bool Iflag, bool Wflag>
inline void _kernel_clip1d_f_wrms(simd_t<T,S> &mean, simd_t<T,S> &rms, const T *intensity, const T *weights, int nfreq, int stride, T *ds_intensity, T *ds_weights)
{
    mean_rms_accumulator<T,S> acc;

    const simd_t<T,S> zero = simd_t<T,S>::zero();
    const simd_t<T,S> one = simd_t<T,S> (1.0);

    for (int ifreq = 0; ifreq < nfreq; ifreq += Df) {
	simd_t<T,S> wival, wval;
	_kernel_downsample<T,S,Df,Dt> (wival, wval, intensity + ifreq*stride, weights + ifreq*stride, stride);

	simd_t<T,S> ival = wival / blendv(wval.compare_gt(zero), wval, one);
	acc.accumulate(ival, wval);

	_write_and_advance_if<T,S,Iflag> (ds_intensity, ival);
	_write_and_advance_if<T,S,Wflag> (ds_weights, wval);
    }

    acc.get_mean_rms(mean, rms);
}


// Iterates on a "strip" of shape (nfreq, S).
// There are no downsampling factors, so in the larger AXIS_FREQ kernel, it will be called with (nfreq/Df) instead of nfreq.
template<typename T, unsigned int S>
inline void _kernel_clip1d_f_iterate(simd_t<T,S> &out_mean, simd_t<T,S> &out_rms, const T *intensity, const T *weights,
				     simd_t<T,S> in_mean, simd_t<T,S> in_thresh, int nfreq, int stride)
{
    mean_rms_accumulator<T,S> acc;

    for (int ifreq = 0; ifreq < nfreq; ifreq++) {
	const T *irow = intensity + ifreq*stride;
	const T *wrow = weights + ifreq*stride;

	simd_t<T,S> ival = simd_t<T,S>::loadu(irow);
	simd_t<T,S> wval = simd_t<T,S>::loadu(wrow);

	simd_t<T,S> ival_c = (ival - in_mean).abs();
	smask_t<T,S> valid = ival_c.compare_lt(in_thresh);

	wval = wval.apply_mask(valid);
	acc.accumulate(ival, wval);
    }

    acc.get_mean_rms(out_mean, out_rms);
}


template<typename T, unsigned int S, unsigned int Df, unsigned int Dt>
inline void _kernel_clip1d_f_mask(T *weights, const T *ds_intensity, simd_t<T,S> mean, simd_t<T,S> thresh, int nfreq, int stride, int ds_stride)
{
    for (int ifreq = 0; ifreq < nfreq; ifreq += Df) {
	simd_t<T,S> ival = simd_t<T,S>::loadu(ds_intensity);
	ival -= mean;
	ival = ival.abs();

	smask_t<T,S> valid = ival.compare_lt(thresh);
	_kernel_mask<T,S,Df,Dt> (weights + ifreq*stride, valid, stride);

	ds_intensity += ds_stride;
    }
}


}  // namespace rf_pipelines

#endif
