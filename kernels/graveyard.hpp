// This file is dead code!
//
// It contains some experiments with alternate kernels, which 
// didn't end up being faster than the originals.
//
// Warning: nothing here is debugged!


// -------------------------------------------------------------------------------------------------
//
// _kernel_detrend_xt_pass1(): similar to _kernel_t_pass1(), but reads precomputed P_l's
//  from an array, rather than computing them on-the-fly.


template<typename T, int S, int N>
inline void _kernel_detrend_xt_pass1(simd_trimatrix<T,S,N> &outm, simd_ntuple<T,S,N> &outv, int nt, const T *ivec, const T *wvec, const T *tmp)
{
    outm.setzero();
    outv.setzero();

    for (int i = 0; i < nt; i += S) {
	simd_t<T,S> ival = simd_t<T,S>::loadu(ivec+i);
	simd_t<T,S> wval = simd_t<T,S>::loadu(wvec+i);

	simd_ntuple<T,S,N> pvec;
	pvec.loadu(tmp + i*N);

	_kernel_detrend_accum_mv(outm, outv, pvec, ival, wval);
    }

    outm.horizontal_sum_in_place();
    outv.horizontal_sum_in_place();
}


// -------------------------------------------------------------------------------------------------
//
// _kernel_merge(): This is a helper for _kernel_unpack_legpoly()
//
// Inputs: simd_ntuple<J> src1, simd_ntuple<K> src2
//
// Output: simd_ntuple<I> dst, defined by
//   dst[:J] = src1[:]
//   dst[J:] = src2[:(I-J)]
//
// The condition J <= I <= (J+K) must be satisfied



template<typename T, int S, int I, int J, int K, typename std::enable_if<(J==I),int>::type = 0>
inline void _kernel_merge(simd_ntuple<T,S,I> &dst, const simd_ntuple<T,S,J> &src1, const simd_ntuple<T,S,K> &src2)
{
    dst = src1;
}

template<typename T, int S, int I, int J, int K, typename std::enable_if<(J<I && I==J+K),int>::type = 0>
inline void _kernel_merge(simd_ntuple<T,S,I> &dst, const simd_ntuple<T,S,J> &src1, const simd_ntuple<T,S,K> &src2)
{
    _kernel_merge(dst.v, src1, src2.v);
    dst.x = src2.x;
}

template<typename T, int S, int I, int J, int K, typename std::enable_if<(J<I && I<J+K),int>::type = 0>
inline void _kernel_merge(simd_ntuple<T,S,I> &dst, const simd_ntuple<T,S,J> &src1, const simd_ntuple<T,S,K> &src2)
{
    _kernel_merge(dst, src1, src2.v);
}


// -------------------------------------------------------------------------------------------------
//
// _kernel_unpack_legpoly(): This is a helper for _kernel_xf_pass1().
//
// It reads an array of (N-1) floats, and fills a simd_ntuple<N> by setting the first
// vector in the tuple to 1, and the remaining vectors to appropriate constant values.


template<typename T, int S, int N, typename std::enable_if<(N==1),int>::type = 0>
inline void _kernel_unpack_legpoly(simd_ntuple<T,S,N> &dst, const T *p)
{
    dst.x = 1.0;
}

template<typename T, int S, int N, typename std::enable_if<(N>1),int>::type = 0>
inline void _kernel_unpack_legpoly(simd_ntuple<T,S,N> &dst, const T *p)
{
    constexpr int M = N - (((N-2) % S) + 1);

    simd_ntuple<T,S,M> tmp;
    _kernel_unpack_legpoly(tmp, p);

    simd_ntuple<T,S,S> tmp2;
    simd_helpers::upsample(tmp2, simd_t<T,S>::loadu(p+M-1));

    _kernel_merge(dst, tmp, tmp2);
}



// -------------------------------------------------------------------------------------------------
//
// _kernel_detrend_xf_pass1(): similar to _kernel_f_pass1(), but reads precomputed P_l's
//  from an array, rather than computing them on-the-fly.


template<typename T, int S, int N>
inline void _kernel_detrend_xf_pass1(simd_trimatrix<T,S,N> &outm, simd_ntuple<T,S,N> &outv, int nfreq, const T *ivec, const T *wvec, int stride, const T *tmp)
{
    constexpr int NP = N - 1 + (S+1 - (N%S)) % S;

    outm.setzero();
    outv.setzero();

    for (int i = 0; i < nfreq; i++) {
	simd_t<T,S> ival = simd_t<T,S>::loadu(ivec + i*stride);
	simd_t<T,S> wval = simd_t<T,S>::loadu(wvec + i*stride);

	simd_ntuple<T,S,N> pvec;
        _kernel_unpack_legpoly(pvec, tmp + i*NP);

	_kernel_detrend_accum_mv(outm, outv, pvec, ival, wval);
    }
}
