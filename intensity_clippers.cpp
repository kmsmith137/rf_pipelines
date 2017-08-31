#include <rf_kernels/intensity_clipper.hpp>
#include <rf_kernels/mean_rms.hpp>

#include "rf_pipelines_internals.hpp"


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


// Externally visible
void weighted_mean_and_rms(float &mean, float &rms, const float *intensity, const float *weights, int nfreq, int nt, int stride, int niter, double sigma, bool two_pass)
{
    rf_kernels::weighted_mean_rms w(nfreq, nt, niter, sigma, two_pass);
    w.compute_wrms(mean, rms, intensity, weights, stride);
}



}  // namespace rf_pipelines
