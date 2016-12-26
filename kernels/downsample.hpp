#ifndef _RF_PIPELINES_KERNELS_DOWNSAMPLE_HPP
#define _RF_PIPELINES_KERNELS_DOWNSAMPLE_HPP

#include <simd_helpers/simd_t.hpp>
#include <simd_helpers/simd_ntuple.hpp>
#include <simd_helpers/udsample.hpp>


namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif

template<typename T, unsigned int S> using simd_t = simd_helpers::simd_t<T,S>;
template<typename T, unsigned int S, unsigned int D> using simd_ntuple = simd_helpers::simd_ntuple<T,S,D>;


// -------------------------------------------------------------------------------------------------
//
// _kernel_downsample1<T,S,R,N> (simd_t<T,S> &ds_wi, simd_t<T,S> &ds_w, const T *intensity, const T *weights, int stride)
//
// Reads a strided array of shape (R,N*S), and sums the result over outer index r
// and middle index N, returning a simd_t<T,S>.


// The "1a" variant accumulates its result
template<typename T, unsigned int S, unsigned int R, unsigned int N, typename std::enable_if<(R==0 || N==0),int>::type = 0>
inline void _kernel_downsample1a(simd_t<T,S> &ds_wi, simd_t<T,S> &ds_w, const T *intensity, const T *weights, int stride)
{
    return;
}


template<typename T, unsigned int S, unsigned int R, unsigned int N, typename std::enable_if<(R > 0 && N > 0),int>::type = 0>
inline void _kernel_downsample1a(simd_t<T,S> &ds_wi, simd_t<T,S> &ds_w, const T *intensity, const T *weights, int stride)
{
    simd_t<T,S> ival = simd_t<T,S>::loadu(intensity);
    simd_t<T,S> wval = simd_t<T,S>::loadu(weights);

    ds_wi += wval * ival;
    ds_w += wval;

    _kernel_downsample1a<T,S,R-1,1> (ds_wi, ds_w, intensity+stride, weights+stride, stride);
    _kernel_downsample1a<T,S,R,N-1> (ds_wi, ds_w, intensity+S, weights+S, stride);
}


template<typename T, unsigned int S, unsigned int R, unsigned int N, typename std::enable_if<(R > 0 && N > 0),int>::type = 0>
inline void _kernel_downsample1(simd_t<T,S> &ds_wi, simd_t<T,S> &ds_w, const T *intensity, const T *weights, int stride)
{
    simd_t<T,S> ival = simd_t<T,S>::loadu(intensity);
    simd_t<T,S> wval = simd_t<T,S>::loadu(weights);
    
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


template<typename T, unsigned int S, unsigned int R, unsigned int D, unsigned int N, typename std::enable_if<(D==0),int>::type = 0>
inline void _kernel_downsample2(simd_ntuple<T,S,D> &ds_wi, simd_ntuple<T,S,D> &ds_w, const T *intensity, const T *weights, int stride)
{
    return;
}


template<typename T, unsigned int S, unsigned int R, unsigned int D, unsigned int N, typename std::enable_if<(D>0),int>::type = 0>
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

 
template<typename T, unsigned int S, unsigned int R, unsigned int D, typename std::enable_if<(D==1),int>::type = 0>
inline void _kernel_downsample(simd_t<T,S> &ds_wi, simd_t<T,S> &ds_w, const T *intensity, const T *weights, int stride)
{
    _kernel_downsample1<T,S,R,1> (ds_wi, ds_w, intensity, weights, stride);
}


template<typename T, unsigned int S, unsigned int R, unsigned int D, typename std::enable_if<(D>1 && D<=S),int>::type = 0>
inline void _kernel_downsample(simd_t<T,S> &ds_wi, simd_t<T,S> &ds_w, const T *intensity, const T *weights, int stride)
{
    simd_ntuple<T,S,D> dsn_wi, dsn_w;
    _kernel_downsample2<T,S,R,D,1> (dsn_wi, dsn_w, intensity, weights, stride);

    ds_wi = downsample(dsn_wi);   // defined in simd_helpers/udsample.hpp
    ds_w = downsample(dsn_w);
}


template<typename T, unsigned int S, unsigned int R, unsigned int D, typename std::enable_if<(D>S),int>::type = 0>
inline void _kernel_downsample(simd_t<T,S> &ds_wi, simd_t<T,S> &ds_w, const T *intensity, const T *weights, int stride)
{
    simd_ntuple<T,S,S> dsn_wi, dsn_w;
    _kernel_downsample2<T,S,R,S,D/S> (dsn_wi, dsn_w, intensity, weights, stride);

    ds_wi = downsample(dsn_wi);
    ds_w = downsample(dsn_w);
}


}  // namespace rf_pipelines

#endif
