#include <rf_kernels/intensity_clipper.hpp>

#include "rf_pipelines_internals.hpp"
#include "kernels/downsample.hpp"
#include "kernels/mean_variance.hpp"

// _kernel_noniterative_wrms, maybe other things?
#include "kernels/intensity_clippers.hpp"


using namespace std;

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif


struct intensity_clipper_transform : public wi_transform 
{
    // (Frequency, time) downsampling factors and axis.
    const int Df;
    const int Dt;
    const rf_kernels::axis_type axis;
    
    // Clipping thresholds.
    const int niter;
    const double sigma;
    const double iter_sigma;
    const bool two_pass;

    // Created in set_stream()
    unique_ptr<rf_kernels::intensity_clipper> kernel;

    
    intensity_clipper_transform(int Df_, int Dt_, rf_kernels::axis_type axis_, int nt_chunk_, double sigma_, int niter_, double iter_sigma_, bool two_pass_)
	: Df(Df_),
	  Dt(Dt_),
	  axis(axis_),
	  niter(niter_),
	  sigma(sigma_),
	  iter_sigma(iter_sigma_ ? iter_sigma_ : sigma_),
	  two_pass(two_pass_)
    {
	stringstream ss;
        ss << "intensity_clipper_cpp(nt_chunk=" << nt_chunk_ << ", axis=" << axis 
	   << ", sigma=" << sigma << ", niter=" << niter << ", iter_sigma=" << iter_sigma 
	   << ", Df=" << Df << ", Dt=" << Dt << ", two_pass=" << two_pass << ")";

	this->name = ss.str();
	this->nt_chunk = nt_chunk_;
	this->nt_prepad = 0;
	this->nt_postpad = 0;

	// Can't construct the kernel yet, since 'nfreq' is not known until set_stream()
	// However, for argument checking purposes, we construct a dummy kernel with Df=nfreq.
	// FIXME eventaully there will be a constructor argument 'allocate=false' that will make sense here.
	rf_kernels::intensity_clipper dummy(Df, nt_chunk, axis, sigma, Df, Dt, niter, iter_sigma, two_pass);
    }

    virtual ~intensity_clipper_transform() { }

    virtual void set_stream(const wi_stream &stream) override
    {
	if (stream.nfreq % Df)
	    throw runtime_error("rf_pipelines::intensity_clipper: stream nfreq (=" + to_string(stream.nfreq) 
				+ ") is not divisible by frequency downsampling factor Df=" + to_string(Df));

	this->nfreq = stream.nfreq;
	this->kernel = make_unique<rf_kernels::intensity_clipper> (nfreq, nt_chunk, axis, sigma, Df, Dt, niter, iter_sigma, two_pass);
    }

    virtual void process_chunk(double t0, double t1, float *intensity, float *weights, ssize_t stride, float *pp_intensity, float *pp_weights, ssize_t pp_stride) override
    {
	this->kernel->clip(intensity, weights, stride);
    }

    virtual void start_substream(int isubstream, double t0) override { }
    virtual void end_substream() override { }
};


// -------------------------------------------------------------------------------------------------


static void check_params(const char *name, int Df, int Dt, rf_kernels::axis_type axis, int nfreq, int nt, int stride, double sigma, int niter, double iter_sigma)
{
    static constexpr int S = constants::single_precision_simd_length;
    static constexpr int MaxDf = constants::max_frequency_downsampling;
    static constexpr int MaxDt = constants::max_time_downsampling;

    if (_unlikely((Df <= 0) || !is_power_of_two(Df)))
	throw runtime_error(string(name) + ": Df=" + to_string(Df) + " must be a power of two");

    if (_unlikely((Dt <= 0) || !is_power_of_two(Dt)))
	throw runtime_error(string(name) + ": Dt=" + to_string(Dt) + " must be a power of two");

    if (_unlikely((axis != rf_kernels::AXIS_FREQ) && (axis != rf_kernels::AXIS_TIME) && (axis != rf_kernels::AXIS_NONE)))
	throw runtime_error(string(name) + ": axis=" + stringify(axis) + " is not defined for this transform");

    if (_unlikely(nfreq <= 0))
	throw runtime_error(string(name) + ": nfreq=" + to_string(nfreq) + ", positive value was expected");

    if (_unlikely(nt <= 0))
	throw runtime_error(string(name) + ": nt=" + to_string(nt) + ", positive value was expected");

    if (_unlikely(abs(stride) < nt))
	throw runtime_error(string(name) + ": stride=" + to_string(stride) + " must be >= nt");

    if (_unlikely(sigma < 1.0))
	throw runtime_error(string(name) + ": sigma=" + to_string(sigma) + " must be >= 1.0");

    if (_unlikely(niter < 1))
	throw runtime_error(string(name) + ": niter=" + to_string(niter) + " must be >= 1");

    if (_unlikely((nfreq % Df) != 0))
	throw runtime_error(string(name) + ": nfreq=" + to_string(nfreq)
			    + " must be a multiple of the downsampling factor Df=" + to_string(Df));
    
    if (_unlikely((nt % (Dt*S)) != 0))
	throw runtime_error(string(name) + ": nt=" + to_string(nt)
			    + " must be a multiple of the downsampling factor Dt=" + to_string(Dt)
			    + " multiplied by constants::single_precision_simd_length=" + to_string(S));

    if (_unlikely((Df > MaxDf) || (Dt > MaxDt)))
	throw runtime_error(string(name) + ": (Df,Dt)=(" + to_string(Df) + "," + to_string(Dt) + ")"
			    + " exceeds compile time limits; to fix this see 'constants' in rf_pipelines.hpp");
}


// Externally visible
shared_ptr<wi_transform> make_intensity_clipper(int nt_chunk, rf_kernels::axis_type axis, double sigma, int niter, double iter_sigma, int Df, int Dt, bool two_pass)
{
    return make_shared<intensity_clipper_transform> (Df, Dt, axis, nt_chunk, sigma, niter, iter_sigma, two_pass);
}


// Externally visible
void apply_intensity_clipper(const float *intensity, float *weights, int nfreq, int nt, int stride, rf_kernels::axis_type axis, double sigma, int niter, double iter_sigma, int Df, int Dt, bool two_pass)
{
    rf_kernels::intensity_clipper ic(nfreq, nt, axis, sigma, Df, Dt, niter, iter_sigma, two_pass);
    ic.clip(intensity, weights, stride);
}


template<typename T, int S>
inline void _weighted_mean_and_rms(simd_t<T,S> &mean, simd_t<T,S> &rms, const float *intensity, const float *weights, int nfreq, int nt, int stride, int niter, double sigma, bool two_pass)
{
    check_params("rf_pipelines: weighted_mean_and_rms()", 1, 1, rf_kernels::AXIS_NONE, nfreq, nt, stride, sigma, niter, sigma);

    if (two_pass)
	_kernel_noniterative_wrms_2d<T,S,1,1,false,false,true> (mean, rms, intensity, weights, nfreq, nt, stride, NULL, NULL);
    else
	_kernel_noniterative_wrms_2d<T,S,1,1,false,false,false> (mean, rms, intensity, weights, nfreq, nt, stride, NULL, NULL);

    _kernel_wrms_iterate_2d<T,S> (mean, rms, intensity, weights, nfreq, nt, stride, niter, sigma);
}

// Externally visible
void weighted_mean_and_rms(float &mean, float &rms, const float *intensity, const float *weights, int nfreq, int nt, int stride, int niter, double sigma, bool two_pass)
{
    static constexpr int S = constants::single_precision_simd_length;

    if (_unlikely(!intensity))
	throw runtime_error("rf_pipelines: weighted_mean_and_rms(): NULL intensity pointer");
    if (_unlikely(!weights))
	throw runtime_error("rf_pipelines: weighted_mean_and_rms(): NULL weights pointer");

    simd_t<float,S> mean_x, rms_x;
    _weighted_mean_and_rms(mean_x, rms_x, intensity, weights, nfreq, nt, stride, niter, sigma, two_pass);

    rms = rms_x.template extract<0> ();
    mean = (rms > 0.0) ? (mean_x.template extract<0> ()) : 0.0;
}



}  // namespace rf_pipelines
