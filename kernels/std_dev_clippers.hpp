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


// Stores temporary buffers which are needed in the std_dev_clipper kernels.
// Note that the allocation logic is not in this file -- it's in allocate_buffers(), defined in std_dev_clippers.cpp
template<typename T>
struct std_dev_clipper_buffers {
    T *sd = NULL;
    smask_t<T,1> *sd_valid = NULL;

    T *ds_intensity = NULL;
    T *ds_weights = NULL;

    ~std_dev_clipper_buffers()
    {
	free(sd);
	free(sd_valid);
	free(ds_intensity);
	free(ds_weights);

	sd = ds_intensity = ds_weights = NULL;
	sd_valid = NULL;
    }
};


// -------------------------------------------------------------------------------------------------
//
// _kernel_std_dev_t()
//
// This is the "bottom line" routine called by std_dev_clipper(AXIS_TIME).
//
// Caller must check (nfreq % Df) == 0 and (nt % (Dt*S)) == 0.
// There is no requirement that (nfreq % (Df*S)) == 0.


template<typename T, unsigned int S, unsigned int Df, unsigned int Dt, bool TwoPass>
inline void _kernel_std_dev_t(const std_dev_clipper_buffers<T> &buf, const T *intensity, const T *weights, int nfreq, int nt, int stride)
{
    T *out_sd = buf.sd;
    smask_t<T,1> *out_valid = buf.sd_valid;

    for (int ifreq = 0; ifreq < nfreq; ifreq += Df) {
	const T *irow = intensity + ifreq * stride;
	const T *wrow = weights + ifreq * stride;

	simd_t<T,S> mean, var;
	_kernel_mean_variance_1d_t<T,S,Df,Dt,false,false,TwoPass> (mean, var, irow, wrow, nt, stride, buf.ds_intensity, buf.ds_weights);

	// scalar instructions should be fine here
	T sd = var.template extract<0> ();
	*out_sd++ = sd;
	*out_valid++ = (sd > 0.0) ? smask_t<T,1>(-1) : 0;
    }
}


template<typename T, unsigned int S, unsigned int Df, unsigned int Dt, bool TwoPass>
inline void _kernel_std_dev_clip_time_axis(const std_dev_clipper_buffers<T> &buf, const T *intensity, T *weights, int nfreq, int nt, int stride, double sigma)
{
    _kernel_std_dev_t<T,S,Df,Dt,TwoPass> (buf, intensity, weights, nfreq, nt, stride);

    clip_1d(nfreq/Df, buf.sd, buf.sd_valid, sigma);

    for (int i = 0; i < nfreq/Df; i++) {
	if (buf.sd_valid[i])
	    continue;

	for (int ifreq = i*Df; ifreq < (i+1)*Df; ifreq++)
	    memset(weights + ifreq*stride, 0, nt * sizeof(T));
    }
}


// -------------------------------------------------------------------------------------------------
//
// _kernel_std_dev_f()
//
// This is the "bottom line" routine called by std_dev_clipper(AXIS_FREQ).
//
// Caller must check (nfreq % Df) == 0 and (nt % (Dt*S)) == 0.
// There is no requirement that (nfreq % (Df*S)) == 0.


template<typename T, unsigned int S, unsigned int Df, unsigned int Dt, bool TwoPass>
inline void _kernel_std_dev_f(const std_dev_clipper_buffers<T> &buf, const T *intensity, const T *weights, int nfreq, int nt, int stride)
{
    const simd_t<T,S> zero = simd_t<T,S>::zero();

    T *out_sd = buf.sd;
    smask_t<T,1> *out_valid = buf.sd_valid;

    for (int it = 0; it < nt; it += Dt*S) {
	const T *icol = intensity + it;
	const T *wcol = weights + it;

	simd_t<T,S> mean, var;
	_kernel_mean_variance_1d_f<T,S,Df,Dt,false,false,TwoPass> (mean, var, icol, wcol, nfreq, stride, buf.ds_intensity, buf.ds_weights);
	
	smask_t<T,S> valid = var.compare_gt(zero);

	var.storeu(out_sd);
	valid.storeu(out_valid);

	out_sd += S;
	out_valid += S;
    }
}


template<typename T, unsigned int S, unsigned int Df, unsigned int Dt, bool TwoPass>
inline void _kernel_std_dev_clip_freq_axis(const std_dev_clipper_buffers<T> &buf, const T *intensity, T *weights, int nfreq, int nt, int stride, double sigma)
{
    _kernel_std_dev_f<T,S,Df,Dt,TwoPass> (buf, intensity, weights, nfreq, nt, stride);

    clip_1d(nt/Dt, buf.sd, buf.sd_valid, sigma);

    _kernel_mask_columns<T,S,Dt> (weights, buf.sd_valid, nfreq, nt, stride);
}


}  // namespace rf_pipelines

#endif
