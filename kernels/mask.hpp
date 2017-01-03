#ifndef _RF_PIPELINES_KERNELS_MASK_HPP
#define _RF_PIPELINES_KERNELS_MASK_HPP

#include <simd_helpers/simd_float32.hpp>
#include <simd_helpers/simd_ntuple.hpp>
#include <simd_helpers/udsample.hpp>


namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif

template<typename T, unsigned int S> using simd_t = simd_helpers::simd_t<T,S>;
template<typename T, unsigned int S, unsigned int D> using simd_ntuple = simd_helpers::simd_ntuple<T,S,D>;
template<typename T, unsigned int S, unsigned int D> using smask_ntuple = simd_helpers::smask_ntuple<T,S,D>;


// -------------------------------------------------------------------------------------------------
//
// _kernel_mask1<T,S,R,N> (T *weights, smask_t<T,S> mask, int stride)
//
// Masks a strided array of shape (R,N*S), by repeating a single simd_t<T,S>.


template<typename T, unsigned int S, unsigned int R, unsigned int N, typename std::enable_if<(R==0 || N==0),int>::type = 0>
inline void _kernel_mask1(T *weights, smask_t<T,S> mask, int stride)
{
    return;
}


template<typename T, unsigned int S, unsigned int R, unsigned int N, typename std::enable_if<(R>0 && N>0),int>::type = 0>
inline void _kernel_mask1(T *weights, smask_t<T,S> mask, int stride)
{
    simd_t<T,S> wval = simd_t<T,S>::loadu(weights);
    wval = wval.apply_mask(mask);
    wval.storeu(weights);

    _kernel_mask1<T,S,R-1,1> (weights+stride, mask, stride);
    _kernel_mask1<T,S,R,N-1> (weights+S, mask, stride);
}


// -------------------------------------------------------------------------------------------------
//
// _kernel_mask2<T,S,R,D,N> (T *weights, smask_ntuple<T,S,D> &mask, int stride)
//
// Masks a strided array of shape (R,D*N*S), by repeating the mask over the (r,n) axes.


template<typename T, unsigned int S, unsigned int R, unsigned int D, unsigned int N, typename std::enable_if<(D==0),int>::type = 0>
inline void _kernel_mask2(T *weights, const smask_ntuple<T,S,D> &mask, int stride)
{
    return;
}


template<typename T, unsigned int S, unsigned int R, unsigned int D, unsigned int N, typename std::enable_if<(D>0),int>::type = 0>
inline void _kernel_mask2(T *weights, const smask_ntuple<T,S,D> &mask, int stride)
{
    _kernel_mask2<T,S,R,D-1,N> (weights, mask.v, stride);
    _kernel_mask1<T,S,R,N> (weights + (D-1)*N*S, mask.x, stride);
}


// -------------------------------------------------------------------------------------------------
//
// _kernel_mask<T,S,R,D> (T *weights, smask_t<T,S> mask, int stride)
//
// Upsamples mask by a factor (R,D) in (freq,time) directions, and applies it to the shape-(R,D*S)
// strided array based at 'weights'.


template<typename T, unsigned int S, unsigned int R, unsigned int D, typename std::enable_if<(D==1),int>::type = 0>
inline void _kernel_mask(T *weights, smask_t<T,S> mask, int stride)
{
    _kernel_mask1<T,S,R,1> (weights, mask, stride);
}


template<typename T, unsigned int S, unsigned int R, unsigned int D, typename std::enable_if<(D>1 && D<=S),int>::type = 0>
inline void _kernel_mask(T *weights, smask_t<T,S> mask, int stride)
{
    smask_ntuple<T,S,D> masku;
    simd_helpers::upsample(masku, mask);

    _kernel_mask2<T,S,R,D,1> (weights, masku, stride);
}


template<typename T, unsigned int S, unsigned int R, unsigned int D, typename std::enable_if<(D>S),int>::type = 0>
inline void _kernel_mask(T *weights, smask_t<T,S> mask, int stride)
{
    smask_ntuple<T,S,S> masku;
    simd_helpers::upsample(masku, mask);

    _kernel_mask2<T,S,R,S,D/S> (weights, masku, stride);
}


}  // namespace rf_pipelines

#endif
