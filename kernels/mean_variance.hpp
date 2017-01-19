#ifndef _RF_PIPELINES_KERNELS_MEAN_VARIANCE_HPP
#define _RF_PIPELINES_KERNELS_MEAN_VARIANCE_HPP

#include <simd_helpers/convert.hpp>
#include "downsample.hpp"

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif


#if 0
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
#endif


// -------------------------------------------------------------------------------------------------


template<typename T_, unsigned int S_, bool Iflag, bool Wflag>
struct _mean_variance_visitor {
    using T = T_;
    static constexpr unsigned int S = S_;

    const simd_t<T,S> zero;
    const simd_t<T,S> one;

    simd_t<T,S> acc0;
    simd_t<T,S> acc1;
    simd_t<T,S> acc2;

    T *ds_intensity;
    T *ds_weights;

    _mean_variance_visitor(T *ds_intensity_, T *ds_weights_) :
	zero(simd_t<T,S>::zero()), one(simd_t<T,S>(1.0))
    {
	acc0 = simd_t<T,S>::zero();
	acc1 = simd_t<T,S>::zero();
	acc2 = simd_t<T,S>::zero();

	ds_intensity = ds_intensity_;
	ds_weights = ds_weights_;
    }

    inline void accumulate_i(simd_t<T,S> ival, simd_t<T,S> wval) 
    { 
	simd_t<T,S> wival = wval * ival;

	acc0 += wval; 
	acc1 += wival;
	acc2 += wival * ival;

	_write_and_advance_if<T,S,Iflag> (ds_intensity, ival);
	_write_and_advance_if<T,S,Wflag> (ds_weights, wval);
    }

    inline void accumulate_wi(simd_t<T,S> wival, simd_t<T,S> wval)
    {
	simd_t<T,S> ival = wival / blendv(wval.compare_gt(zero), wval, one);

	acc0 += wval;
	acc1 += wival;
	acc2 += wival * ival;

	_write_and_advance_if<T,S,Iflag> (ds_intensity, ival);
	_write_and_advance_if<T,S,Wflag> (ds_weights, wval);
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

	valid = acc0.compare_gt(zero);

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
	simd_t<T,S> variance;
	get_mean_variance(mean, variance);
	rms = variance.sqrt();
    }
};


// -------------------------------------------------------------------------------------------------


template<unsigned int Df, unsigned int Dt, typename V, typename std::enable_if<((Df>1) || (Dt>1)),int>::type = 0>
inline void _kernel_visit_2d(V &v, const typename V::T *intensity, const typename V::T *weights, int nfreq, int nt, int stride)
{
    using T = typename V::T;
    constexpr int S = V::S;

    for (int ifreq = 0; ifreq < nfreq; ifreq += Df) {
	const T *irow = intensity + ifreq*stride;
	const T *wrow = weights + ifreq*stride;

	for (int it = 0; it < nt; it += Dt*S) {
	    simd_t<T,S> wival, wval;
	    _kernel_downsample<T,S,Df,Dt> (wival, wval, irow + it, wrow + it, stride);

	    v.accumulate_wi(wival, wval);
	}
    }

    v.horizontal_sum();
}


template<unsigned int Df, unsigned int Dt, typename V, typename std::enable_if<((Df==1) && (Dt==1)),int>::type = 0>
inline void _kernel_visit_2d(V &v, const typename V::T *intensity, const typename V::T *weights, int nfreq, int nt, int stride)
{
    using T = typename V::T;
    constexpr int S = V::S;

    for (int ifreq = 0; ifreq < nfreq; ifreq += Df) {
	const T *irow = intensity + ifreq*stride;
	const T *wrow = weights + ifreq*stride;

	for (int it = 0; it < nt; it += S) {
	    simd_t<T,S> ival = simd_t<T,S>::loadu(irow + it);
	    simd_t<T,S> wval = simd_t<T,S>::loadu(wrow + it);

	    v.accumulate_i(ival, wval);
	}
    }

    v.horizontal_sum();
}


template<unsigned int Df, unsigned int Dt, typename V, typename std::enable_if<((Df>1) || (Dt>1)),int>::type = 0>
inline void _kernel_visit_1d_f(V &v, const typename V::T *intensity, const typename V::T *weights, int nfreq, int stride)
{
    using T = typename V::T;
    constexpr int S = V::S;

    for (int ifreq = 0; ifreq < nfreq; ifreq += Df) {
	simd_t<T,S> wival, wval;
	_kernel_downsample<T,S,Df,Dt> (wival, wval, intensity + ifreq*stride, weights + ifreq*stride, stride);

	v.accumulate_wi(wival, wval);
    }
}


template<unsigned int Df, unsigned int Dt, typename V, typename std::enable_if<((Df==1) && (Dt==1)),int>::type = 0>
inline void _kernel_visit_1d_f(V &v, const typename V::T *intensity, const typename V::T *weights, int nfreq, int stride)
{
    using T = typename V::T;
    constexpr int S = V::S;

    for (int ifreq = 0; ifreq < nfreq; ifreq++) {
	simd_t<T,S> ival = simd_t<T,S>::loadu(intensity + ifreq*stride);
	simd_t<T,S> wval = simd_t<T,S>::loadu(weights + ifreq*stride);

	v.accumulate_i(ival, wval);
    }
}


}  // namespace rf_pipelines

#endif
