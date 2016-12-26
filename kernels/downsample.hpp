#ifndef _RF_PIPELINES_KERNELS_DOWNSAMPLE_HPP
#define _RF_PIPELINES_KERNELS_DOWNSAMPLE_HPP

#include <simd_helpers/simd_t.hpp>
#include <simd_helpers/simd_ntuple.hpp>
#include <simd_helpers/udsample.hpp>


namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif


// -------------------------------------------------------------------------------------------------
//
// simd_t<T,S> _kernel_downsample1<T,S,R,N> (const T *p, int stride)
//
// Reads a strided array of shape (R,N*S), and sums the result over outer index r
// and middle index N, returning a simd_t<T,S>.


template<typename T, unsigned int S, unsigned int R, unsigned int N, typename enable_if<(R==0 || N==0),int>::type = 0>
inline simd_t<T,S> _kernel_downsample1(const T *p, int stride, simd_t<T,S> x)
{
    return x;
}

template<typename T, unsigned int S, unsigned int R, unsigned int N, typename enable_if<(R > 0 && N > 0),int>::type = 0>
inline simd_t<T,S> _kernel_downsample1(const T *p, int stride, simd_t<T,S> x)
{
    x = _kernel_downsample1<T,S,R-1,1> (p+stride, stride, x + simd_t<T,S>::loadu(p));
    x = _kernel_downsample1<T,S,R,N-1> (p+S, stride, x);
    return x;
}

template<typename T, unsigned int S, unsigned int R, unsigned int N>
inline simd_t<T,S> _kernel_downsample1(const T *p, int stride)
{
    simd_t<T,S> x = simd_t<T,S>::loadu(p);
    x = _kernel_downsample1<T,S,R-1,1> (p+stride, stride, x);
    x = _kernel_downsample1<T,S,R,N-1> (p+S, stride, x);
    return x;
}


// -------------------------------------------------------------------------------------------------
//
// simd_ntuple<T,S> _kernel_downsample2<T,S,R,D,N> (const T *p, int stride)
//
// Reads a strided array of shape (R,D*N*S), and sums the result over outer index r
// and middle index n, returning a simd_ntuple<T,S,D>.


template<typename T, unsigned int S, unsigned int R, unsigned int D, unsigned int N, typename enable_if<(D==0),int>::type = 0>
inline simd_ntuple<T,S,D> _kernel_downsample2(const T *p, int stride)
{
    return simd_ntuple<T,S,0> ();
}


template<typename T, unsigned int S, unsigned int R, unsigned int D, unsigned int N, typename enable_if<(D>0),int>::type = 0>
inline simd_ntuple<T,S,D> _kernel_downsample2(const T *p, int stride)
{
    simd_ntuple<T,S,D> ret;
    ret.v = _kernel_downsample2<T,S,R,D-1,N> (p, stride);
    ret.x = _kernel_downsample1<T,S,R,N> (p+(D-1)*N*S, stride);
    return ret;
}



// -------------------------------------------------------------------------------------------------
//
// simd_t<T,S> _kernel_downsample<T,S,R,D> (const T *p, int stride)
//
// Reads a strided array of shape (R,D*S), and sums the result over outer index r and inner index d, 
// returning a simd_t<T,S,D>.


template<typename T, unsigned int S, unsigned int R, unsigned int D, typename enable_if<(D==1),int>::type = 0>
inline simd_t<T,S> _kernel_downsample(const T *p, int stride)
{
    return _kernel_downsample1<T,S,R,1> (p, stride);
}

template<typename T, unsigned int S, unsigned int R, unsigned int D, typename enable_if<(D>1 && D<=S),int>::type = 0>
inline simd_t<T,S> _kernel_downsample(const T *p, int stride)
{
    simd_ntuple<T,S,D> t = _kernel_downsample2<T,S,R,D,1> (p, stride);
    return downsample(t);  // defined in simd_helpers.hpp
}

template<typename T, unsigned int S, unsigned int R, unsigned int D, typename enable_if<(D>S,int>::type = 0>
inline simd_t<T,S> _kernel_downsample(const T *p, int stride)
{
    simd_ntuple<T,S,S> t = _kernel_downsample2<T,S,R,S,D/S> (p, stride);
    return downsample(t);  // defined in simd_helpers.hpp
}


}  // namespace rf_pipelines

#endif
