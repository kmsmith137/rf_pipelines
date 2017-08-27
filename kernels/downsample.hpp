#ifndef _RF_PIPELINES_KERNELS_DOWNSAMPLE_HPP
#define _RF_PIPELINES_KERNELS_DOWNSAMPLE_HPP

#include <simd_helpers/simd_float32.hpp>
#include <simd_helpers/simd_ntuple.hpp>
#include <simd_helpers/udsample.hpp>


namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif

template<typename T, int S> using simd_t = simd_helpers::simd_t<T,S>;
template<typename T, int S, int D> using simd_ntuple = simd_helpers::simd_ntuple<T,S,D>;


// -------------------------------------------------------------------------------------------------
//
// _kernel_downsample1<T,S,R,N> (simd_t<T,S> &ds_wi, simd_t<T,S> &ds_w, const T *intensity, const T *weights, int stride)
//
// Reads a strided array of shape (R,N*S), and sums the result over outer index r
// and middle index N, returning a simd_t<T,S>.


// The "1a" variant accumulates its result
template<typename T, int S, int R, int N, typename std::enable_if<(R==0 || N==0),int>::type = 0>
inline void _kernel_downsample1a(simd_t<T,S> &ds_wi, simd_t<T,S> &ds_w, const T *intensity, const T *weights, int stride)
{
    return;
}


template<typename T, int S, int R, int N, typename std::enable_if<(R > 0 && N > 0),int>::type = 0>
inline void _kernel_downsample1a(simd_t<T,S> &ds_wi, simd_t<T,S> &ds_w, const T *intensity, const T *weights, int stride)
{
    simd_t<T,S> ival = simd_helpers::simd_load<T,S> (intensity);
    simd_t<T,S> wval = simd_helpers::simd_load<T,S> (weights);

    ds_wi += wval * ival;
    ds_w += wval;

    _kernel_downsample1a<T,S,R-1,1> (ds_wi, ds_w, intensity+stride, weights+stride, stride);
    _kernel_downsample1a<T,S,R,N-1> (ds_wi, ds_w, intensity+S, weights+S, stride);
}


template<typename T, int S, int R, int N, typename std::enable_if<(R > 0 && N > 0),int>::type = 0>
inline void _kernel_downsample1(simd_t<T,S> &ds_wi, simd_t<T,S> &ds_w, const T *intensity, const T *weights, int stride)
{
    simd_t<T,S> ival = simd_helpers::simd_load<T,S> (intensity);
    simd_t<T,S> wval = simd_helpers::simd_load<T,S> (weights);
    
    ds_wi = wval * ival;
    ds_w = wval;
    
    _kernel_downsample1a<T,S,R-1,1> (ds_wi, ds_w, intensity+stride, weights+stride, stride);
    _kernel_downsample1a<T,S,R,N-1> (ds_wi, ds_w, intensity+S, weights+S, stride);
}


// -------------------------------------------------------------------------------------------------
//
// _kernel_downsample2<T,S,R,D,N> (simd_ntuple<T,S,D> &ds_wi, simd_ntuple<T,S,D> &ds_w, const T *intensity, const T *weights, int stride)
//
// Reads a strided array of shape (R,D*N*S), and sums the result over outer index r
// and middle index n, returning a simd_ntuple<T,S,D>.


template<typename T, int S, int R, int D, int N, typename std::enable_if<(D==0),int>::type = 0>
inline void _kernel_downsample2(simd_ntuple<T,S,D> &ds_wi, simd_ntuple<T,S,D> &ds_w, const T *intensity, const T *weights, int stride)
{
    return;
}


template<typename T, int S, int R, int D, int N, typename std::enable_if<(D>0),int>::type = 0>
inline void _kernel_downsample2(simd_ntuple<T,S,D> &ds_wi, simd_ntuple<T,S,D> &ds_w, const T *intensity, const T *weights, int stride)
{
    _kernel_downsample2<T,S,R,D-1,N> (ds_wi.v, ds_w.v, intensity, weights, stride);
    _kernel_downsample1<T,S,R,N> (ds_wi.x, ds_w.x, intensity + (D-1)*N*S, weights + (D-1)*N*S, stride);
}


// -------------------------------------------------------------------------------------------------
//
// _kernel_downsample<T,S,R,D> (simd_t<T,S> &ds_wi, simd_t<T,S> &ds_w, const T *intensity, const T *weights, int stride)
//
// Reads a strided array of shape (R,D*S), and sums the result over outer index r and inner index d, 
// returning a simd_t<T,S>.

 
template<typename T, int S, int R, int D, typename std::enable_if<(D==1),int>::type = 0>
inline void _kernel_downsample(simd_t<T,S> &ds_wi, simd_t<T,S> &ds_w, const T *intensity, const T *weights, int stride)
{
    _kernel_downsample1<T,S,R,1> (ds_wi, ds_w, intensity, weights, stride);
}


template<typename T, int S, int R, int D, typename std::enable_if<(D>1 && D<=S),int>::type = 0>
inline void _kernel_downsample(simd_t<T,S> &ds_wi, simd_t<T,S> &ds_w, const T *intensity, const T *weights, int stride)
{
    simd_ntuple<T,S,D> dsn_wi, dsn_w;
    _kernel_downsample2<T,S,R,D,1> (dsn_wi, dsn_w, intensity, weights, stride);

    ds_wi = simd_helpers::downsample(dsn_wi);   // defined in simd_helpers/udsample.hpp
    ds_w = simd_helpers::downsample(dsn_w);
}


template<typename T, int S, int R, int D, typename std::enable_if<(D>S),int>::type = 0>
inline void _kernel_downsample(simd_t<T,S> &ds_wi, simd_t<T,S> &ds_w, const T *intensity, const T *weights, int stride)
{
    simd_ntuple<T,S,S> dsn_wi, dsn_w;
    _kernel_downsample2<T,S,R,S,D/S> (dsn_wi, dsn_w, intensity, weights, stride);

    ds_wi = simd_helpers::downsample(dsn_wi);
    ds_w = simd_helpers::downsample(dsn_w);
}


// -------------------------------------------------------------------------------------------------
//
// _kernel_downsample_2d<T,S,Df,Dt> (out_intensity, out_weights, out_stride, in_intensity, in_nfreq, in_nt, in_stride)
//
// Caller must check that nfreq is divisible by Df, and nt is divisible by (Dt*S).
//
// This is the kernel which gets called in the externally visible function wi_downsample().


template<typename T, int S, int Df, int Dt>
inline void _kernel_downsample_2d(T *out_intensity, T *out_weights, int out_stride, const T *in_intensity, const T *in_weights, int in_nfreq, int in_nt, int in_stride)
{
    const simd_t<T,S> zero = simd_t<T,S>::zero();
    const simd_t<T,S> one = simd_t<T,S> (1.0);

    int out_nfreq = in_nfreq / Df;
    int out_nt = in_nt / Dt;

    for (int ifreq = 0; ifreq < out_nfreq; ifreq++) {
	T *out_irow = out_intensity + ifreq * out_stride;
	T *out_wrow = out_weights + ifreq * out_stride;

	const T *in_irow = in_intensity + (ifreq*Df) * in_stride;
	const T *in_wrow = in_weights + (ifreq*Df) * in_stride;

	for (int it = 0; it < out_nt; it += S) {
	    simd_t<T,S> ds_wival, ds_wval;
	    _kernel_downsample<T,S,Df,Dt> (ds_wival, ds_wval, in_irow + it*Dt, in_wrow + it*Dt, in_stride);

	    simd_t<T,S> ds_ival = ds_wival / blendv(ds_wval.compare_gt(zero), ds_wval, one);
	    ds_ival.storeu(out_irow + it);
	    ds_wval.storeu(out_wrow + it);
	}
    }
}


}  // namespace rf_pipelines

#endif
