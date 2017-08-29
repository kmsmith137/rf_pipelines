// FIXME: currently we need to compile a new kernel for every (Df,Dt) pair, where
// Df,Dt are the frequency/time downsampling factors.  Eventually I'd like to 
// improve this by having special kernels to handle the large-Df and large-Dt cases.

#include <cassert>
#include <rf_kernels/std_dev_clipper.hpp>
#include "rf_pipelines_internals.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif


struct std_dev_clipper_transform : public wi_transform 
{
    // (Frequency, time) downsampling factors and axis.
    const int Df;
    const int Dt;
    const rf_kernels::axis_type axis;
    const bool two_pass;
    
    // Clipping threshold.
    const double sigma;

    unique_ptr<rf_kernels::std_dev_clipper> kernel;

    std_dev_clipper_transform(int Df_, int Dt_, rf_kernels::axis_type axis_, int nt_chunk_, double sigma_, bool two_pass_)
	: Df(Df_), Dt(Dt_), axis(axis_), two_pass(two_pass_), sigma(sigma_)
    {
	stringstream ss;
        ss << "std_dev_clipper_transform_cpp(nt_chunk=" << nt_chunk_ << ", axis=" << axis
           << ", sigma=" << sigma << ", Df=" << Df << ", Dt=" << Dt << ", two_pass=" << two_pass << ")";

	// FIXME more argument checking would be good here
	
        this->name = ss.str();
	this->nt_chunk = nt_chunk_;
	this->nt_prepad = 0;
	this->nt_postpad = 0;
    }

    virtual void set_stream(const wi_stream &stream) override
    {
	if (stream.nfreq % Df)
	    throw runtime_error("rf_pipelines std_dev_clipper: stream nfreq (=" + to_string(stream.nfreq) 
				+ ") is not divisible by frequency downsampling factor Df=" + to_string(Df));

	this->nfreq = stream.nfreq;
	this->kernel = make_unique<rf_kernels::std_dev_clipper> (nfreq, nt_chunk, axis, Df, Dt, sigma, two_pass);
    }

    virtual void process_chunk(double t0, double t1, float *intensity, float *weights, ssize_t stride, float *pp_intensity, float *pp_weights, ssize_t pp_stride) override
    {
	this->kernel->clip(intensity, weights, stride);
    }

    virtual void start_substream(int isubstream, double t0) override { }
    virtual void end_substream() override { }
};


// Externally callable
shared_ptr<wi_transform> make_std_dev_clipper(int nt_chunk, rf_kernels::axis_type axis, double sigma, int Df, int Dt, bool two_pass)
{
    return make_shared<std_dev_clipper_transform> (Df, Dt, axis, nt_chunk, sigma, two_pass);
}


// Externally callable
void apply_std_dev_clipper(const float *intensity, float *weights, int nfreq, int nt, int stride, rf_kernels::axis_type axis, double sigma, int Df, int Dt, bool two_pass)
{
    rf_kernels::std_dev_clipper sd(nfreq, nt, axis, Df, Dt, sigma, two_pass);
    sd.clip(intensity, weights, stride);
}


}  // namespace rf_pipelines
