#ifndef _RF_PIPELINES_KERNELS_MEAN_VARIANCE_HPP
#define _RF_PIPELINES_KERNELS_MEAN_VARIANCE_HPP

#include <simd_helpers/convert.hpp>
#include "downsample.hpp"

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif

template<typename T, unsigned int S> using simd_t = simd_helpers::simd_t<T,S>;
template<typename T, unsigned int S> using smask_t = simd_helpers::smask_t<T,S>;



// -------------------------------------------------------------------------------------------------
//
// "Visit" kernels: the V arugment is an inline "visitor" class which will be defined next in this file.


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


// -------------------------------------------------------------------------------------------------
//
// The rest of this file defines "visitor" classes:
//   _mean_variance_visitor
//   _mean_variance_iterator
//   _mean_visitor
//   _variance_visitor


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

	// Branches will be optimized out at compile-time
	if (Iflag) {
	    ival.storeu(ds_intensity);
	    ds_intensity += S;
	}
	
	if (Wflag) {
	    wval.storeu(ds_weights);
	    ds_weights += S;
	}
    }

    inline void accumulate_wi(simd_t<T,S> wival, simd_t<T,S> wval)
    {
	simd_t<T,S> ival = wival / blendv(wval.compare_gt(zero), wval, one);

	acc0 += wval;
	acc1 += wival;
	acc2 += wival * ival;

	// Branches will be optimized out at compile-time
	if (Iflag) {
	    ival.storeu(ds_intensity);
	    ds_intensity += S;
	}
	
	if (Wflag) {
	    wval.storeu(ds_weights);
	    ds_weights += S;
	}
    }

    inline void horizontal_sum()
    {
	acc0 = acc0.horizontal_sum();
	acc1 = acc1.horizontal_sum();
	acc2 = acc2.horizontal_sum();
    }

    inline void get_mean_variance(simd_t<T,S> &mean, simd_t<T,S> &var) const
    {
	static constexpr T eps_3 = 1.0e3 * simd_helpers::machine_epsilon<T> ();

	smask_t<T,S> valid = acc0.compare_gt(zero);
	simd_t<T,S> t0 = blendv(valid, acc0, one);
	mean = acc1 / t0;

	simd_t<T,S> mean2 = mean * mean;
	var = acc2/t0 - mean2;

	simd_t<T,S> thresh = simd_t<T,S>(eps_3) * mean2;

	valid = valid.bitwise_and(var.compare_gt(thresh));
	var = var.apply_mask(valid);
    }

    inline void get_mean_rms(simd_t<T,S> &mean, simd_t<T,S> &rms) const
    {
	simd_t<T,S> variance;
	get_mean_variance(mean, variance);
	rms = variance.sqrt();
    }
};


// -------------------------------------------------------------------------------------------------


template<typename T, unsigned int S, bool TwoPass> 
struct _mean_variance_iterator;


// One-pass iterator
template<typename T_, unsigned int S_>
struct _mean_variance_iterator<T_,S_,false> {
    using T = T_;
    static constexpr unsigned int S = S_;

    simd_t<T,S> in_mean;
    simd_t<T,S> in_thresh;

    simd_t<T,S> acc0;
    simd_t<T,S> acc1;
    simd_t<T,S> acc2;

    _mean_variance_iterator(simd_t<T,S> in_mean_, simd_t<T,S> in_thresh_)
    {
	in_mean = in_mean_;
	in_thresh = in_thresh_;

	acc0 = simd_t<T,S>::zero();
	acc1 = simd_t<T,S>::zero();
	acc2 = simd_t<T,S>::zero();
    }

    inline void accumulate_i(simd_t<T,S> ival, simd_t<T,S> wval)
    {
	simd_t<T,S> ival_c = (ival - in_mean).abs();
	smask_t<T,S> valid = ival_c.compare_lt(in_thresh);
	
	wval = wval.apply_mask(valid);
	simd_t<T,S> wival = wval * ival;

	acc0 += wval;
	acc1 += wival;
	acc2 += wival * ival;
    }

    inline void horizontal_sum()
    {
	acc0 = acc0.horizontal_sum();
	acc1 = acc1.horizontal_sum();
	acc2 = acc2.horizontal_sum();
    }

    inline void get_mean_rms(simd_t<T,S> &out_mean, simd_t<T,S> &out_rms)
    {
	static constexpr T eps = 1.0e3 * simd_helpers::machine_epsilon<T> ();
	const simd_t<T,S> zero = simd_t<T,S>::zero();
	const simd_t<T,S> one = 1.0;

	smask_t<T,S> valid = acc0.compare_gt(zero);

	simd_t<T,S> t0 = blendv(valid, acc0, one);
	out_mean = acc1/t0;

	simd_t<T,S> out_mean2 = out_mean * out_mean;
	simd_t<T,S> var = acc2/t0 - out_mean2;

	simd_t<T,S> thresh = simd_t<T,S>(eps) * out_mean2;
	valid = valid.bitwise_and(var.compare_gt(thresh));
	var = var.apply_mask(valid);

	out_rms = var.sqrt();
    }
};


// -------------------------------------------------------------------------------------------------


template<typename T_, unsigned int S_, bool Iflag, bool Wflag>
struct _mean_visitor {
    using T = T_;
    static constexpr unsigned int S = S_;

    const simd_t<T,S> zero;
    const simd_t<T,S> one;

    simd_t<T,S> acc0;
    simd_t<T,S> acc1;

    T *ds_intensity;
    T *ds_weights;

    _mean_visitor(T *ds_intensity_, T *ds_weights_) :
	zero(simd_t<T,S>::zero()), one(simd_t<T,S>(1.0))
    {
	acc0 = simd_t<T,S>::zero();
	acc1 = simd_t<T,S>::zero();

	ds_intensity = ds_intensity_;
	ds_weights = ds_weights_;
    }

    inline void accumulate_i(simd_t<T,S> ival, simd_t<T,S> wval) 
    { 
	acc0 += wval; 
	acc1 += wval * ival;

	if (Iflag) {
	    ival.storeu(ds_intensity);
	    ds_intensity += S;
	}
	
	if (Wflag) {
	    wval.storeu(ds_weights);
	    ds_weights += S;
	}
    }

    inline void accumulate_wi(simd_t<T,S> wival, simd_t<T,S> wval)
    {
	acc0 += wval;
	acc1 += wival;

	if (Iflag) {
	    simd_t<T,S> ival = wival / blendv(wval.compare_gt(zero), wval, one);
	    ival.storeu(ds_intensity);
	    ds_intensity += S;
	}
	
	if (Wflag) {
	    wval.storeu(ds_weights);
	    ds_weights += S;
	}
    }

    inline void horizontal_sum()
    {
	acc0 = acc0.horizontal_sum();
	acc1 = acc1.horizontal_sum();
    }

    inline simd_t<T,S> get_mean() const
    {
	smask_t<T,S> valid = acc0.compare_gt(zero);
	return acc1 / blendv(valid, acc0, one);
    }
};


// -------------------------------------------------------------------------------------------------


template<typename T_, unsigned int S_>
struct _variance_visitor {
    using T = T_;
    static constexpr unsigned int S = S_;

    const simd_t<T,S> zero;
    const simd_t<T,S> one;

    simd_t<T,S> acc0;
    simd_t<T,S> acc2;
    simd_t<T,S> in_mean;

    _variance_visitor(simd_t<T,S> in_mean_) :
	zero(simd_t<T,S>::zero()), one(simd_t<T,S>(1.0))
    {
	acc0 = simd_t<T,S>::zero();
	acc2 = simd_t<T,S>::zero();
	in_mean = in_mean_;
    }

    inline void accumulate_i(simd_t<T,S> ival, simd_t<T,S> wval) 
    { 
	ival -= in_mean;
	acc0 += wval; 
	acc2 += wval * ival * ival;
    }

    inline void horizontal_sum()
    {
	acc0 = acc0.horizontal_sum();
	acc2 = acc2.horizontal_sum();
    }

    inline simd_t<T,S> get_variance() const
    {
	static constexpr T eps_2 = 1.0e2 * simd_helpers::machine_epsilon<T> ();

	smask_t<T,S> valid = acc0.compare_gt(zero);
	simd_t<T,S> t0 = blendv(valid, acc0, one);
	simd_t<T,S> var = acc2/t0;

	simd_t<T,S> thresh = simd_t<T,S>(eps_2) * in_mean;
	thresh = thresh * thresh;

	valid = valid.bitwise_and(var.compare_gt(thresh));
	return var.apply_mask(valid);
    }
};


// -------------------------------------------------------------------------------------------------


// Two-pass iterator (really a misnomer!)
template<typename T_, unsigned int S_>
struct _mean_variance_iterator<T_,S_,true> {
    using T = T_;
    static constexpr unsigned int S = S_;

    const simd_t<T,S> zero;
    const simd_t<T,S> one;

    simd_t<T,S> in_mean;
    simd_t<T,S> in_thresh;

    simd_t<T,S> acc0;
    simd_t<T,S> acc1;
    simd_t<T,S> acc2;

    _mean_variance_iterator(simd_t<T,S> in_mean_, simd_t<T,S> in_thresh_) :
	zero(simd_t<T,S>::zero()), one(simd_t<T,S>(1.0))
    {
	in_mean = in_mean_;
	in_thresh = in_thresh_;

	acc0 = simd_t<T,S>::zero();
	acc1 = simd_t<T,S>::zero();
	acc2 = simd_t<T,S>::zero();
    }

    inline void accumulate_i(simd_t<T,S> ival, simd_t<T,S> wval)
    {
	ival -= in_mean;
	smask_t<T,S> valid = ival.abs().compare_lt(in_thresh);
	
	wval = wval.apply_mask(valid);
	simd_t<T,S> wival = wval * ival;

	acc0 += wval;
	acc1 += wival;
	acc2 += wival * ival;
    }

    inline void horizontal_sum()
    {
	acc0 = acc0.horizontal_sum();
	acc1 = acc1.horizontal_sum();
	acc2 = acc2.horizontal_sum();
    }

    inline void get_mean_rms(simd_t<T,S> &out_mean, simd_t<T,S> &out_rms)
    {
	static constexpr T eps_2 = 1.0e2 * simd_helpers::machine_epsilon<T> ();
	static constexpr T eps_3 = 1.0e3 * simd_helpers::machine_epsilon<T> ();

	smask_t<T,S> valid = acc0.compare_gt(zero);

	simd_t<T,S> t0 = blendv(valid, acc0, one);
	simd_t<T,S> dmean = acc1/t0;

	simd_t<T,S> dmean2 = dmean * dmean;
	simd_t<T,S> var = acc2/t0 - dmean2;

	simd_t<T,S> thresh1 = simd_t<T,S>(eps_2) * in_mean;
	simd_t<T,S> thresh2 = simd_t<T,S>(eps_3) * dmean2;
	thresh2 = thresh2.max(thresh1 * thresh1);

	valid = valid.bitwise_and(var.compare_gt(thresh2));
	var = var.apply_mask(valid);

	out_mean = in_mean + dmean;
	out_rms = var.sqrt();
    }
};


// -------------------------------------------------------------------------------------------------


// _kernel_mean_variance_2d() case 1: one-pass version
template<typename T, unsigned int S, unsigned int Df, unsigned int Dt, bool Iflag, bool Wflag, bool TwoPass, typename std::enable_if<(!TwoPass),int>::type = 0>
inline void _kernel_mean_variance_2d(simd_t<T,S> &mean, simd_t<T,S> &var, const T *intensity, const T *weights, int nfreq, int nt, int stride, T *ds_intensity, T *ds_weights)
{
    _mean_variance_visitor<T,S,Iflag,Wflag> v(ds_intensity, ds_weights);
    _kernel_visit_2d<Df,Dt> (v, intensity, weights, nfreq, nt, stride);
    v.get_mean_variance(mean, var);
}

// _kernel_mean_variance_2d() case 2: two-pass version, downsampled
template<typename T, unsigned int S, unsigned int Df, unsigned int Dt, bool Iflag, bool Wflag, bool TwoPass, typename std::enable_if<(TwoPass && ((Df>1) || (Dt>1))),int>::type = 0>
inline void _kernel_mean_variance_2d(simd_t<T,S> &mean, simd_t<T,S> &var, const T *intensity, const T *weights, int nfreq, int nt, int stride, T *ds_intensity, T *ds_weights)
{
    _mean_visitor<T,S,true,true> v(ds_intensity, ds_weights);
    _kernel_visit_2d<Df,Dt> (v, intensity, weights, nfreq, nt, stride);
    mean = v.get_mean();

    _variance_visitor<T,S> vv(mean);
    _kernel_visit_2d<1,1> (vv, ds_intensity, ds_weights, nfreq/Df, nt/Dt, nt/Dt);
    var = vv.get_variance();
}

// _kernel_mean_variance_2d() case 3: two-pass version, non-downsampled
template<typename T, unsigned int S, unsigned int Df, unsigned int Dt, bool Iflag, bool Wflag, bool TwoPass, typename std::enable_if<(TwoPass && (Df==1) && (Dt==1)),int>::type = 0>
inline void _kernel_mean_variance_2d(simd_t<T,S> &mean, simd_t<T,S> &var, const T *intensity, const T *weights, int nfreq, int nt, int stride, T *ds_intensity, T *ds_weights)
{
    _mean_visitor<T,S,Iflag,Wflag> v(ds_intensity, ds_weights);
    _kernel_visit_2d<1,1> (v, intensity, weights, nfreq, nt, stride);
    mean = v.get_mean();

    _variance_visitor<T,S> vv(mean);
    _kernel_visit_2d<1,1> (vv, intensity, weights, nfreq, nt, stride);
    var = vv.get_variance();
}


// _kernel_mean_variance_1d_f() case 1: one-pass version
template<typename T, unsigned int S, unsigned int Df, unsigned int Dt, bool Iflag, bool Wflag, bool TwoPass, typename std::enable_if<(!TwoPass),int>::type = 0>
inline void _kernel_mean_variance_1d_f(simd_t<T,S> &mean, simd_t<T,S> &var, const T *intensity, const T *weights, int nfreq, int stride, T *ds_intensity, T *ds_weights)
{
    _mean_variance_visitor<T,S,Iflag,Wflag> v(ds_intensity, ds_weights);
    _kernel_visit_1d_f<Df,Dt> (v, intensity, weights, nfreq, stride);
    v.get_mean_variance(mean, var);
}

// _kernel_mean_variance_1d_f() case 2: two-pass version, downsampled
template<typename T, unsigned int S, unsigned int Df, unsigned int Dt, bool Iflag, bool Wflag, bool TwoPass, typename std::enable_if<(TwoPass && ((Df>1) || (Dt>1))),int>::type = 0>
inline void _kernel_mean_variance_1d_f(simd_t<T,S> &mean, simd_t<T,S> &var, const T *intensity, const T *weights, int nfreq, int stride, T *ds_intensity, T *ds_weights)
{
    _mean_visitor<T,S,true,true> v(ds_intensity, ds_weights);
    _kernel_visit_1d_f<Df,Dt> (v, intensity, weights, nfreq, stride);
    mean = v.get_mean();

    _variance_visitor<T,S> vv(mean);
    _kernel_visit_1d_f<1,1> (vv, ds_intensity, ds_weights, nfreq/Df, S);
    var = vv.get_variance();
}

// _kernel_mean_variance_1d_f() case 3: two-pass version, non-downsampled
template<typename T, unsigned int S, unsigned int Df, unsigned int Dt, bool Iflag, bool Wflag, bool TwoPass, typename std::enable_if<(TwoPass && (Df==1) && (Dt==1)),int>::type = 0>
inline void _kernel_mean_variance_1d_f(simd_t<T,S> &mean, simd_t<T,S> &var, const T *intensity, const T *weights, int nfreq, int stride, T *ds_intensity, T *ds_weights)
{
    _mean_visitor<T,S,Iflag,Wflag> v(ds_intensity, ds_weights);
    _kernel_visit_1d_f<1,1> (v, intensity, weights, nfreq, stride);
    mean = v.get_mean();

    _variance_visitor<T,S> vv(mean);
    _kernel_visit_1d_f<1,1> (vv, intensity, weights, nfreq, stride);    
    var = vv.get_variance();
}


// Placeholder for future expansion
template<typename T, unsigned int S, unsigned int Df, unsigned int Dt, bool Iflag, bool Wflag, bool TwoPass>
inline void _kernel_mean_variance_1d_t(simd_t<T,S> &mean, simd_t<T,S> &var, const T *intensity, const T *weights, int nt, int stride, T *ds_intensity, T *ds_weights)
{
    _kernel_mean_variance_2d<T,S,Df,Dt,Iflag,Wflag,TwoPass> (mean, var, intensity, weights, Df, nt, stride, ds_intensity, ds_weights);
}


}  // namespace rf_pipelines

#endif
