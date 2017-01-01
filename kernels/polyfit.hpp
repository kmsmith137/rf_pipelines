#ifndef _RF_PIPELINES_KERNELS_POLYFIT_HPP
#define _RF_PIPELINES_KERNELS_POLYFIT_HPP

#include <cstring>
#include <simd_helpers/simd_t.hpp>
#include <simd_helpers/simd_ntuple.hpp>
#include <simd_helpers/simd_trimatrix.hpp>

// Branch predictor hint
#ifndef _unlikely
#define _unlikely(cond)  (__builtin_expect(cond,0))
#endif


namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif

template<typename T, unsigned int S> using simd_t = simd_helpers::simd_t<T,S>;
template<typename T, unsigned int S, unsigned int D> using simd_ntuple = simd_helpers::simd_ntuple<T,S,D>;
template<typename T, unsigned int S, unsigned int N> using simd_trimatrix = simd_helpers::simd_trimatrix<T,S,N>;


// -------------------------------------------------------------------------------------------------


template<typename T, unsigned int S>
inline void _kernel_legpoly_eval(simd_ntuple<T,S,1> &pl, simd_t<T,S> z)
{
    pl.x = 1.0;
}

template<typename T, unsigned int S>
inline void _kernel_legpoly_eval(simd_ntuple<T,S,2> &pl, simd_t<T,S> z)
{
    pl.v.x = 1.0;
    pl.x = z;
}

template<typename T, unsigned int S, unsigned int N, typename std::enable_if<(N > 2),int>::type = 0>
inline void _kernel_legpoly_eval(simd_ntuple<T,S,N> &pl, simd_t<T,S> z)
{
     // P_N(z) = a z P_{N-1}(z) + b P_{N-2}(z)
    constexpr T a = double(2*N-3) / double(N-1);
    constexpr T b = -double(N-2) / double(N-1);

    _kernel_legpoly_eval(pl.v, z);
    pl.x = a * z * pl.v.x + b * pl.v.v.x;
}


// -------------------------------------------------------------------------------------------------
//
// _kernel_detrend_accum_mv(outm, outv, pvec, ival, wval)
//
// Accumulates contribution to (M,v) in first pass of detrender.


template<typename T, unsigned int S>
inline void _kernel_detrend_accum_mv(simd_trimatrix<T,S,0> &outm, simd_ntuple<T,S,0> &outv, const simd_ntuple<T,S,0> &pl, simd_t<T,S> ival, simd_t<T,S> wval) { }

template<typename T, unsigned int S, unsigned int N, typename std::enable_if<(N>0),int>::type = 0>
inline void _kernel_detrend_accum_mv(simd_trimatrix<T,S,N> &outm, simd_ntuple<T,S,N> &outv, const simd_ntuple<T,S,N> &pvec, simd_t<T,S> ival, simd_t<T,S> wval)
{
    _kernel_detrend_accum_mv(outm.m, outv.v, pvec.v, ival, wval);

    simd_t<T,S> wp = wval * pvec.x;

    outm.v += wp * pvec;
    outv.x += wp * ival;
}


// -------------------------------------------------------------------------------------------------
//
// _kernel_detrend_t<T,S,N> (nfreq, nt, intensity, weights, stride, epsilon)
//
// Detrend along time (=fastest varying) axis of 2D strided array.
// Note: the degree of the polynomial fit is (N-1), not N!


template<typename T, unsigned int S, unsigned int N>
inline void _kernel_detrend_t_pass1(simd_trimatrix<T,S,N> &outm, simd_ntuple<T,S,N> &outv, int nt, const T *ivec, const T *wvec)
{
    outm.setzero();
    outv.setzero();

    simd_t<T,S> z0 = simd_t<T,S>::range();
    z0 -= simd_t<T,S>(0.5 * (nt-1));
    z0 *= simd_t<T,S>(2.0 / T(nt));

    T dz = 2.0 / T(nt);

    for (int i = 0; i < nt; i += S) {
	simd_t<T,S> z = z0 + dz * simd_t<T,S>(i);

	simd_t<T,S> ival = simd_t<T,S>::loadu(ivec+i);
	simd_t<T,S> wval = simd_t<T,S>::loadu(wvec+i);

	simd_ntuple<T,S,N> pvec;
	_kernel_legpoly_eval(pvec, z);

	_kernel_detrend_accum_mv(outm, outv, pvec, ival, wval);
    }

    outm.horizontal_sum_in_place();
    outv.horizontal_sum_in_place();
}


template<typename T, unsigned int S, unsigned int N>
inline void _kernel_detrend_t_pass2(float *ivec, int nt, const simd_ntuple<T,S,N> &coeffs)
{
    simd_t<T,S> z0 = simd_t<T,S>::range();
    z0 -= simd_t<T,S>(0.5 * (nt-1));
    z0 *= simd_t<T,S>(2.0 / T(nt));

    T dz = 2.0 / T(nt);

    for (int i = 0; i < nt; i += S) {
	simd_t<T,S> z = z0 + dz * simd_t<T,S>(i);

	simd_ntuple<T,S,N> pvec;
	_kernel_legpoly_eval(pvec, z);

	simd_t<T,S> ival = simd_t<T,S>::loadu(ivec + i);

	ival = pvec._vertical_dotn(coeffs, ival);
	ival.storeu(ivec + i);
    }
}


template<typename T, unsigned int S, unsigned int N>
inline void _kernel_detrend_t(int nfreq, int nt, T *intensity, T *weights, int stride, double epsilon=1.0e-2)
{
    // Caller should have asserted this already, but rechecking here should have negligible overhead
    if (_unlikely((nt % S) != 0))
	throw std::runtime_error("rf_pipelines internal error: nt is not divisible by S in _kernel_detrend_t()");

    for (int ifreq = 0; ifreq < nfreq; ifreq++) {
	T *ivec = intensity + ifreq * stride;
	T *wvec = weights + ifreq * stride;

	simd_trimatrix<T,S,N> xmat;
	simd_ntuple<T,S,N> xvec;

	_kernel_detrend_t_pass1(xmat, xvec, nt, ivec, wvec);

	simd_t<int,S> flags = xmat.cholesky_in_place_checked(epsilon);

	if (flags.is_all_zeros()) {
	    // Case 1: Cholesky factorization was badly conditioned
	    memset(weights + ifreq*stride, 0, nt * sizeof(T));
	    continue;
	}

	// Case 2: Cholesky factorization is numerically stable, polynomial fitting can be performed
	xmat.solve_lower_in_place(xvec);
	xmat.solve_upper_in_place(xvec);

	_kernel_detrend_t_pass2(ivec, nt, xvec);
    }
}


// -------------------------------------------------------------------------------------------------
//
// _kernel_detrend_f<T,S,N> (nfreq, nt, intensity, weights, stride, epsilon)
//
// Detrend along frequency (=slowest varying) axis of 2D strided array.
// Note: the degree of the polynomial fit is (N-1), not N!


template<typename T, unsigned int S, unsigned int N>
inline void _kernel_detrend_f_pass1(simd_trimatrix<T,S,N> &outm, simd_ntuple<T,S,N> &outv, int nfreq, const T *ivec, const T *wvec, int stride)
{
    outm.setzero();
    outv.setzero();

    T z0 = -(nfreq-1) / T(nfreq);
    T dz = 2.0 / T(nfreq);

    for (int i = 0; i < nfreq; i++) {
	T z = z0 + i*dz;

	simd_t<T,S> ival = simd_t<T,S>::loadu(ivec + i*stride);
	simd_t<T,S> wval = simd_t<T,S>::loadu(wvec + i*stride);

	simd_ntuple<T,S,N> pvec;
	_kernel_legpoly_eval(pvec, simd_t<T,S>(z));
	_kernel_detrend_accum_mv(outm, outv, pvec, ival, wval);
    }
}

template<typename T, unsigned int S, unsigned int N>
inline void _kernel_detrend_f_pass2(float *ivec, int nfreq, const simd_ntuple<T,S,N> &coeffs, int stride)
{
    T z0 = -(nfreq-1) / T(nfreq);
    T dz = 2.0 / T(nfreq);

    for (int i = 0; i < nfreq; i++) {
	T z = z0 + i*dz;
    
	simd_ntuple<T,S,N> pvec;
	_kernel_legpoly_eval(pvec, simd_t<T,S>(z));

	simd_t<T,S> ival = simd_t<T,S>::loadu(ivec + i*stride);

	ival = pvec._vertical_dotn(coeffs, ival);
	ival.storeu(ivec + i*stride);
    }
}


// Zeros a complete block of S columns in the 'weights' array.
template<typename T, unsigned int S>
inline void _kernel_colzero_full(float *weights, int nfreq, int stride)
{
    simd_t<T,S> z = simd_t<T,S>::zero();

    for (int i = 0; i < nfreq; i++)
	z.storeu(weights + i*stride);
}


// Zeros a partial block of S columns in the 'weights' array.
// Each word in the 'mask' array should be either 0 or -1=0xff..
template<typename T, unsigned int S>
inline void _kernel_colzero_partial(float *weights, int nfreq, int stride, simd_t<int,S> mask)
{
    for (int i = 0; i < nfreq; i++) {
	simd_t<T,S> w = simd_t<T,S>::loadu(weights + i*stride);
	w = w.bitwise_and(mask);
	w.storeu(weights + i*stride);
    }
}


template<typename T, unsigned int S, unsigned int N>
inline void _kernel_detrend_f(int nfreq, int nt, T *intensity, T *weights, int stride, double epsilon=1.0e-2)
{
    // Caller should have asserted this already, but rechecking here should have negligible overhead
    if (_unlikely((nt % S) != 0))
	throw std::runtime_error("rf_pipelines internal error: nt is not divisible by S in _kernel_detrend_f()");

    for (int it = 0; it < nt; it += S) {
	T *ivec = intensity + it;
	T *wvec = weights + it;

	simd_trimatrix<T,S,N> xmat;
	simd_ntuple<T,S,N> xvec;

	_kernel_detrend_f_pass1(xmat, xvec, nfreq, ivec, wvec, stride);

	simd_t<int,S> flags = xmat.cholesky_in_place_checked(epsilon);

	if (flags.is_all_zeros()) {
	    // If we get here, then all columns of the weights array should be zeroed,
	    // but the intensity array can be left unmodified.
	    _kernel_colzero_full<T,S> (wvec, nfreq, stride);
	    continue;
	}

	// In columns where the fit is poorly conditioned, we leave the intensity array unmodified.
	xvec = xvec.bitwise_and(flags);

	xmat.solve_lower_in_place(xvec);
	xmat.solve_upper_in_place(xvec);
	
	_kernel_detrend_f_pass2(ivec, nfreq, xvec, stride);

	if (flags.is_all_ones())
	    continue;

	// If we get here, then partial zeroing of the weights array is needed.
	_kernel_colzero_partial<T,S> (wvec, nfreq, stride, flags);
    }
}


}  // namespace rf_pipelines

#endif
