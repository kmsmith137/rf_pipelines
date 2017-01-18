#ifndef _RF_PIPELINES_KERNELS_MEAN_RMS_ACCUMULATOR_HPP
#define _RF_PIPELINES_KERNELS_MEAN_RMS_ACCUMULATOR_HPP

#include <simd_helpers/convert.hpp>
#include "downsample.hpp"

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
//
// Warning: the case (Td,Sd) != (T,S) has not actually been tested!!

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

    // This version of accumulate() is better if the caller already has values 
    // of W, (W*I), and (W*I^2).  For example, if a "subaccumulator" is being accumulated.
    inline void accumulate(simd_t<T,S> wval, simd_t<T,S> wival, simd_t<T,S> wiival)
    {
	acc0 += wval;
	acc1 += wival;
	acc2 += wiival;
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

    // Warning: the case (Td,Sd) != (T,S) has not actually been tested!!
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

    inline void get_mean_variance(simd_t<T,S> &mean, simd_t<T,S> &var, smask_t<T,S> &valid) const
    {
	static constexpr T eps = 1.0e3 * simd_helpers::machine_epsilon<T> ();

	valid = acc0.compare_gt(simd_t<T,S>::zero());

	simd_t<T,S> t0 = blendv(valid, acc0, simd_t<T,S>(1.0));
	mean = acc1 / t0;

	simd_t<T,S> mean2 = mean * mean;
	var = acc2/t0 - mean2;

	simd_t<T,S> thresh = simd_t<T,S>(eps) * mean2;
	valid = valid.bitwise_and(var.compare_gt(thresh));
	var = var.apply_mask(valid);
    }

    inline void get_mean_variance(simd_t<T,S> &mean, simd_t<T,S> &var) const
    {
	smask_t<T,S> valid;
	get_mean_variance(mean, var, valid);
    }

    inline void get_mean_rms(simd_t<T,S> &mean, simd_t<T,S> &rms) const
    {
	simd_t<T,S> var;
	get_mean_variance(mean, var);
	rms = var.sqrt();
    }
};


// -------------------------------------------------------------------------------------------------
//
// _kernel_mean_rms_accumulate_2d<T, S, Df, Dt, Iflag, Wflag>
//     (mean_rms_accumulator<T,S> &acc, const T *intensity, const T *weights,
//      int nfreq, int nt, int stride, T *ds_intensity, T *ds_weights)
//
// Computes the weighted mean and rms of a 2D strided array,
// with downsampling factors (Df,Dt) in the (frequency,time) axes.
//
// As the downsampled intensity and weights arrays are computed, they are written to
// 'ds_intensity' and 'ds_weights'.  These are unstrided arrays, i.e. the row stride
// is (nt/Dt).
//
// The Iflag, Wflag template arguments will omit writing the ds_intensity, ds_weights
// arrays if set to 'false'.  In this case, passing a NULL pointer is OK.
//
// Caller must check that:
//    - nfreq is divisible by Df
//    - nt is divisible by (Dt*S)       NOTE (Dt*S) here, not (Dt) !!
//
// FIXME a more general version of this kernel is possible, with a template parameter R
// which controlls the number of rows read in each pass.



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


template<typename T, unsigned int S, unsigned int Df, unsigned int Dt, bool Iflag, bool Wflag, typename std::enable_if<((Df>1) || (Dt>1)),int>::type = 0>
inline void _kernel_mean_rms_accumulate_2d(mean_rms_accumulator<T,S> &acc, const T *intensity, const T *weights, int nfreq, int nt, int stride, T *ds_intensity, T *ds_weights)
{
    const simd_t<T,S> zero = simd_t<T,S>::zero();
    const simd_t<T,S> one = simd_t<T,S> (1.0);

    for (int ifreq = 0; ifreq < nfreq; ifreq += Df) {
	const T *irow = intensity + ifreq*stride;
	const T *wrow = weights + ifreq*stride;

	for (int it = 0; it < nt; it += Dt*S) {
	    simd_t<T,S> wival, wval;
	    _kernel_downsample<T,S,Df,Dt> (wival, wval, irow + it, wrow + it, stride);

	    simd_t<T,S> ival = wival / blendv(wval.compare_gt(zero), wval, one);
	    acc.accumulate(wval, wival, wival * ival);

	    _write_and_advance_if<T,S,Iflag> (ds_intensity, ival);
	    _write_and_advance_if<T,S,Wflag> (ds_weights, wval);
	}
    }
}


template<typename T, unsigned int S, unsigned int Df, unsigned int Dt, bool Iflag, bool Wflag, typename std::enable_if<((Df==1) && (Dt==1)),int>::type = 0>
inline void _kernel_mean_rms_accumulate_2d(mean_rms_accumulator<T,S> &acc, const T *intensity, const T *weights, int nfreq, int nt, int stride, T *ds_intensity, T *ds_weights)
{
    static_assert(!Iflag && !Wflag, "mean_rms_accumulate() with (Df,Dt)=(1,1) and Iflag/Wflag set to true: this is probably a bug");
    
    for (int ifreq = 0; ifreq < nfreq; ifreq++) {
	const T *irow = intensity + ifreq*stride;
	const T *wrow = weights + ifreq*stride;

	for (int it = 0; it < nt; it += S) {
	    simd_t<T,S> ival = simd_t<T,S>::loadu(irow + it);
	    simd_t<T,S> wval = simd_t<T,S>::loadu(wrow + it);

	    acc.accumulate(ival, wval);
	}
    }
}


// -------------------------------------------------------------------------------------------------
//
// _kernel_mean_rms_accumulate_1d_t<T, S, Df, Dt, Iflag, Wflag>
//     (mean_rms_accumulator<T,S> &acc, const T *intensity, const T *weights,
//      int nt, int stride, T *ds_intensity, T *ds_weights)
//
// For now, this is just the special case of _kernel_mean_rms_accumulate_2d() with nfreq=Df.
// It's just here as a placeholder, in case I decide to implement something different some day.


template<typename T, unsigned int S, unsigned int Df, unsigned int Dt, bool Iflag=false, bool Wflag=false>
inline void _kernel_mean_rms_accumulate_1d_t(mean_rms_accumulator<T,S> &acc, const T *intensity, const T *weights, int nt, int stride, T *ds_intensity = NULL, T *ds_weights = NULL)
{
    _kernel_mean_rms_accumulate_2d<T,S,Df,Dt,Iflag,Wflag> (acc, intensity, weights, Df, nt, stride, ds_intensity, ds_weights);
}


// -------------------------------------------------------------------------------------------------
//
// _kernel_mean_rms_acumulate_1d_f<T, S, Df, Dt, Iflag, Wflag>
//    (mean_rms_accumulator<T,S> &acc, const T *intensity, const T *weights,
//     int nfreq, int stride, T *ds_intensity, T *ds_weights)


template<typename T, unsigned int S, unsigned int Df, unsigned int Dt, bool Iflag=false, bool Wflag=false, typename std::enable_if<((Df>1) || (Dt>1)),int>::type = 0>
inline void _kernel_mean_rms_accumulate_1d_f(mean_rms_accumulator<T,S> &acc, const T *intensity, const T *weights, int nfreq, int stride, T *ds_intensity = NULL, T *ds_weights = NULL)
{
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
}


template<typename T, unsigned int S, unsigned int Df, unsigned int Dt, bool Iflag=false, bool Wflag=false, typename std::enable_if<((Df==1) && (Dt==1)),int>::type = 0>
inline void _kernel_mean_rms_accumulate_1d_f(mean_rms_accumulator<T,S> &acc, const T *intensity, const T *weights, int nfreq, int stride, T *ds_intensity = NULL, T *ds_weights = NULL)
{
    static_assert(!Iflag && !Wflag, "mean_rms_accumulate() with (Df,Dt)=(1,1) and Iflag/Wflag set to true: this is probably a bug");

    for (int ifreq = 0; ifreq < nfreq; ifreq++) {    
	simd_t<T,S> ival = simd_t<T,S>::loadu(intensity + ifreq*stride);
	simd_t<T,S> wval = simd_t<T,S>::loadu(weights + ifreq*stride);

	acc.accumulate(ival, wval);
    }
}


}  // namespace rf_pipelines

#endif
