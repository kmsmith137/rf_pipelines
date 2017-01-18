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


// -------------------------------------------------------------------------------------------------
//
// _kernel_noniterative_wrms_2d<T, S, Df, Dt, Iflag, Wflag>
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


template<typename T, unsigned int S, unsigned int Df, unsigned int Dt, bool Iflag, bool Wflag>
inline void _kernel_noniterative_wrms_2d(simd_t<T,S> &mean, simd_t<T,S> &rms, const T *intensity, const T *weights, int nfreq, int nt, int stride, T *ds_intensity, T *ds_weights)
{
    mean_rms_accumulator<T,S> acc;
    _kernel_mean_rms_accumulate_2d<T,S,Df,Dt,Iflag,Wflag> (acc, intensity, weights, nfreq, nt, stride, ds_intensity, ds_weights);

    acc.horizontal_sum();

    simd_t<T,S> mean_i, rms_i;
    acc.get_mean_rms(mean, rms);
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

    _kernel_noniterative_wrms_2d<T,S,Df,Dt,DsiFlag,DswFlag> (mean, rms, intensity, weights, nfreq, nt, stride, ds_int, ds_wt);

    const T *s_intensity = DsiFlag ? ds_int : intensity;
    const T *s_weights = DswFlag ? ds_wt : weights;
    int s_stride = DsiFlag ? (nt/Dt) : stride;   // must use DsiFlag here, not DswFlag
	
    for (int iter = 1; iter < niter; iter++) {
	simd_t<T,S> thresh = simd_t<T,S>(sigma) * rms;
	_kernel_clip2d_iterate<T,S> (mean, rms, s_intensity, s_weights, mean, thresh, nfreq/Df, nt/Dt, s_stride);
    }
}


template<unsigned int S, unsigned int Df, unsigned int Dt, bool IterFlag>
static void _kernel_nds_2d(int &nds_int, int &nds_wt, int nfreq, int nt)
{
    static constexpr bool DsiFlag = (Df > 1) || (Dt > 1);
    static constexpr bool DswFlag = IterFlag && ((Df > 1) || (Dt > 1));

    nds_int = DsiFlag ? ((nfreq/Df) * (nt/Dt)) : 0;
    nds_wt = DswFlag ? ((nfreq/Df) * (nt/Dt)) : 0;
}


template<typename T, unsigned int S, unsigned int Df, unsigned int Dt, bool IterFlag>
inline void _kernel_clip_2d(const T *intensity, T *weights, int nfreq, int nt, int stride, int niter, double sigma, double iter_sigma, T *ds_int, T *ds_wt)
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
//
// _kernel_clip_1d_t(): This is the "bottom line" routine which is wrapped by intensity_clipper(AXIS_TIME)
//
// Currently implemented by calling the 2d kernels many times with nt=1.
// FIXME there is a little extra overhead here, should improve by writing real 1D kernels.


template<unsigned int S, unsigned int Df, unsigned int Dt, bool IterFlag>
static void _kernel_nds_1d_t(int &nds_int, int &nds_wt, int nfreq, int nt)
{
    static constexpr bool DsiFlag = (Df > 1) || (Dt > 1);
    static constexpr bool DswFlag = IterFlag && ((Df > 1) || (Dt > 1));

    nds_int = DsiFlag ? (nt/Dt) : 0;
    nds_wt = DswFlag ? (nt/Dt) : 0;
}


template<unsigned int S, unsigned int Df, unsigned int Dt, bool IterFlag>
static void _kernel_clip_1d_t(const float *intensity, float *weights, int nfreq, int nt, int stride, int niter, double sigma, double iter_sigma, float *ds_int, float *ds_wt)
{
    static constexpr bool DsiFlag = (Df > 1) || (Dt > 1);
    static constexpr bool DswFlag = IterFlag && ((Df > 1) || (Dt > 1));

    for (int ifreq = 0; ifreq < nfreq; ifreq += Df) {
	const float *irow = intensity + ifreq * stride;
	float *wrow = weights + ifreq * stride;
	
	// We pass nfreq=Df to _kernel_clip2d_wrms, not the "true" nfreq
	simd_t<float,S> mean, rms;
	_kernel_noniterative_wrms_2d<float,S,Df,Dt,DsiFlag,DswFlag> (mean, rms, irow, wrow, Df, nt, stride, ds_int, ds_wt);
	
	const float *irow2 = DsiFlag ? ds_int : irow;
	const float *wrow2 = DswFlag ? ds_wt : wrow;
	
	for (int iter = 1; iter < niter; iter++) {
	    // Here we pass nfreq=1 and stride=0
	    simd_t<float,S> thresh = simd_t<float,S>(iter_sigma) * rms;
	    _kernel_clip2d_iterate<float,S> (mean, rms, irow2, wrow2, mean, thresh, 1, nt/Dt, 0);    // (irow2, wrow2, iter_sigma)
	}

	// Here we pass nfreq=Df.  Setting both strides to 'stride' is OK but this isn't completely obvious.
	simd_t<float,S> thresh = simd_t<float,S>(sigma) * rms;
	_kernel_clip2d_mask<float,S,Df,Dt> (wrow, irow2, mean, thresh, Df, nt, stride, stride);    // (wrow, irow2, sigma)
    }
}


// -------------------------------------------------------------------------------------------------


template<typename T, unsigned int S, unsigned int Df, unsigned int Dt, bool Iflag, bool Wflag>
inline void _kernel_clip1d_f_wrms(simd_t<T,S> &mean, simd_t<T,S> &rms, const T *intensity, const T *weights, int nfreq, int stride, T *ds_intensity, T *ds_weights)
{
    mean_rms_accumulator<T,S> acc;
    _kernel_mean_rms_accumulate_1d_f<T,S,Df,Dt,Iflag,Wflag> (acc, intensity, weights, nfreq, stride, ds_intensity, ds_weights);

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


// -------------------------------------------------------------------------------------------------
//
// _kernel_clip_1d_f(): This is the "bottom line" routine which is wrapped by intensity_clipper(AXIS_FREQ)
//
// Currently implemented by calling the 2d kernels many times with nfreq=1.
// FIXME there is a little extra overhead here, should improve by writing real 1D kernels.


template<unsigned int S, unsigned int Df, unsigned int Dt, bool IterFlag>
static void _kernel_nds_1d_f(int &nds_int, int &nds_wt, int nfreq, int nt)
{
    static constexpr bool DsiFlag = (Df > 1) || (Dt > 1);
    static constexpr bool DswFlag = IterFlag && ((Df > 1) || (Dt > 1));

    nds_int = DsiFlag ? ((nfreq/Df) * S) : 0;
    nds_wt = DswFlag ? ((nfreq/Df) * S) : 0;
}


template<unsigned int S, unsigned int Df, unsigned int Dt, bool IterFlag>
static void _kernel_clip_1d_f(const float *intensity, float *weights, int nfreq, int nt, int stride, int niter, double sigma, double iter_sigma, float *ds_int, float *ds_wt)
{
    static constexpr bool DsiFlag = (Df > 1) || (Dt > 1);
    static constexpr bool DswFlag = IterFlag && ((Df > 1) || (Dt > 1));

    for (int it = 0; it < nt; it += Dt*S) {
	const float *icol = intensity + it;
	float *wcol = weights + it;
	
	simd_t<float,S> mean, rms;	
	_kernel_clip1d_f_wrms<float,S,Df,Dt,DsiFlag,DswFlag> (mean, rms, icol, wcol, nfreq, stride, ds_int, ds_wt);
	
	const float *icol2 = DsiFlag ? ds_int : icol;
	const float *wcol2 = DswFlag ? ds_wt : wcol;
	int stride2 = DsiFlag ? S : stride;   // must use DsiFlag here, not DswFlag
	
	for (int iter = 1; iter < niter; iter++) {
	    simd_t<float,S> thresh = simd_t<float,S>(iter_sigma) * rms;
	    _kernel_clip1d_f_iterate<float,S> (mean, rms, icol2, wcol2, mean, thresh, nfreq/Df, stride2);    // (irow2, wrow2, iter_sigma)
	}

	simd_t<float,S> thresh = simd_t<float,S>(sigma) * rms;
	_kernel_clip1d_f_mask<float,S,Df,Dt> (wcol, icol2, mean, thresh, nfreq, stride, stride2);
    }
}


}  // namespace rf_pipelines

#endif
