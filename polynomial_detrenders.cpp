#include "rf_kernels/polynomial_detrender.hpp"
#include "rf_pipelines_internals.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif


// Just being pedantic, since rf_pipelines defines enum { AXIS_FREQ, AXIS_TIME, ... },
// whereas rf_kernels uses axis=0,1 for frequency,time respectively.
static_assert(AXIS_FREQ==0, "This implementation assumes AXIS_FREQ==0");
static_assert(AXIS_TIME==1, "This implementation assumes AXIS_TIME==1");


struct polynomial_detrender : public wi_transform
{
    rf_kernels::polynomial_detrender kernel;
    const double epsilon;

    polynomial_detrender(int axis, int nt_chunk_, int polydeg, double epsilon_) :
	kernel(axis, polydeg),
	epsilon(epsilon_)
    {
	stringstream ss;
        ss << "polynomial_detrender_cpp(nt_chunk=" << nt_chunk_ << ", axis=" << axis << ", polydeg=" << polydeg << ", epsilon=" << epsilon_ << ")";

	this->name = ss.str();
	this->nt_chunk = nt_chunk_;
	this->nt_prepad = 0;
	this->nt_postpad = 0;
    }
    
    virtual void set_stream(const wi_stream &stream) override
    {
	this->nfreq = stream.nfreq;
    }

    virtual void process_chunk(double t0, double t1, float *intensity, float *weights, ssize_t stride, float *pp_intensity, float *pp_weights, ssize_t pp_stride) override
    {
	this->kernel.detrend(nfreq, nt_chunk, intensity, weights, stride, epsilon);
    }

    virtual void start_substream(int isubstream, double t0) override { }
    virtual void end_substream() override { }


    virtual Json::Value serialize_to_json() const override
    {
	Json::Value ret;

	ret["transform_name"] = "polynomial_detrender";
	ret["nt_chunk"] = int(this->nt_chunk);
	ret["axis"] = axis_type_to_string(this->kernel.axis);
	ret["polydeg"] = this->kernel.polydeg;
	ret["epsilon"] = this->epsilon;

	return ret;
    }
};


// Externally callable factory function
shared_ptr<wi_transform> make_polynomial_detrender(int nt_chunk, axis_type axis, int polydeg, double epsilon)
{
    return make_shared<polynomial_detrender> (axis, nt_chunk, polydeg, epsilon);
}


void apply_polynomial_detrender(float *intensity, float *weights, int nfreq, int nt, int stride, axis_type axis, int polydeg, double epsilon)
{
    rf_kernels::polynomial_detrender kernel(axis, polydeg);
    kernel.detrend(nfreq, nt, intensity, weights, stride, epsilon);
}


}  // namespace rf_pipelines
