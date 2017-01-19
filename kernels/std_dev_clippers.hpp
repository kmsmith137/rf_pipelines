#ifndef _RF_PIPELINES_KERNELS_STD_DEV_CLIPPERS_HPP
#define _RF_PIPELINES_KERNELS_STD_DEV_CLIPPERS_HPP

#include "downsample.hpp"
#include "mean_variance.hpp"

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif


// Defined in std_dev_clippers.cpp
extern void clip_1d(int n, float *tmp_sd, smask_t<float,1> *tmp_valid, double sigma);


// -------------------------------------------------------------------------------------------------
//
// _kernel_std_dev_t()
//
// The output arrays 'out_sd' and 'out_valid' are 1D arrays of shape (nfreq/Df).
//
// Caller must check (nfreq % Df) == 0 and (nt % (Dt*S)) == 0.
// There is no requirement that (nfreq % (Df*S)) == 0.


template<typename T, unsigned int S, unsigned int Df, unsigned int Dt>
inline void _kernel_std_dev_t(T *out_sd, smask_t<T,1> *out_valid, const T *intensity, const T *weights, int nfreq, int nt, int stride)
{
    for (int ifreq = 0; ifreq < nfreq; ifreq += Df) {
	const T *irow = intensity + ifreq * stride;
	const T *wrow = weights + ifreq * stride;

	simd_t<T,S> mean, var;
	_kernel_mean_variance_1d_t<T,S,Df,Dt,false,false,false> (mean, var, irow, wrow, nt, stride, NULL, NULL);

	// scalar instructions should be fine here
	T sd = var.template extract<0> ();
	*out_sd++ = sd;
	*out_valid++ = (sd > 0.0) ? smask_t<T,1>(-1) : 0;
    }
}


template<unsigned int S, unsigned int Df, unsigned int Dt>
inline void _kernel_std_dev_clip_time_axis(const float *intensity, float *weights, int nfreq, int nt, int stride, double sigma, float *tmp_sd, smask_t<float,1> *tmp_valid)
{
    _kernel_std_dev_t<float,S,Df,Dt> (tmp_sd, tmp_valid, intensity, weights, nfreq, nt, stride);

    clip_1d(nfreq/Df, tmp_sd, tmp_valid, sigma);

    for (int i = 0; i < nfreq/Df; i++) {
	if (tmp_valid[i])
	    continue;

	for (int ifreq = i*Df; ifreq < (i+1)*Df; ifreq++)
	    memset(weights + ifreq*stride, 0, nt * sizeof(float));
    }
}


// -------------------------------------------------------------------------------------------------
//
// _kernel_std_dev_f()
//
// The output arrays 'out_sd' and 'out_valid' are 1D arrays of shape (nt/Dt).
//
// Caller must check (nfreq % Df) == 0 and (nt % (Dt*S)) == 0.
// There is no requirement that (nfreq % (Df*S)) == 0.


template<typename T, unsigned int S, unsigned int Df, unsigned int Dt>
inline void _kernel_std_dev_f(T *out_sd, smask_t<T,1> *out_valid, const T *intensity, const T *weights, int nfreq, int nt, int stride)
{
    for (int it = 0; it < nt; it += Dt*S) {
	const T *icol = intensity + it;
	const T *wcol = weights + it;

	simd_t<T,S> mean, var;
	_kernel_mean_variance_1d_f<T,S,Df,Dt,false,false,false> (mean, var, icol, wcol, nfreq, stride, NULL, NULL);
	
	smask_t<T,S> valid = var.compare_gt(simd_t<T,S>::zero());

	var.storeu(out_sd);
	valid.storeu(out_valid);

	out_sd += S;
	out_valid += S;
    }
}


template<unsigned int S, unsigned int Df, unsigned int Dt>
inline void _kernel_std_dev_clip_freq_axis(const float *intensity, float *weights, int nfreq, int nt, int stride, double sigma, float *tmp_sd, smask_t<float,1> *tmp_valid)
{
    _kernel_std_dev_f<float,S,Df,Dt> (tmp_sd, tmp_valid, intensity, weights, nfreq, nt, stride);

    clip_1d(nt/Dt, tmp_sd, tmp_valid, sigma);

    _kernel_mask_columns<float,S,Dt> (weights, tmp_valid, nfreq, nt, stride);
}


}  // namespace rf_pipelines

#endif
