#include "rf_kernels/spline_detrender.hpp"
#include "rf_pipelines_internals.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif


struct spline_detrender : public wi_transform
{
    const int nbins;
    const double epsilon;
    std::unique_ptr<rf_kernels::spline_detrender> kernel;

    spline_detrender(int nt_chunk_, int nbins_, double epsilon_) :
	nbins(nbins_),
	epsilon(epsilon_)
    {
	stringstream ss;
        ss << "spline_detrender(nt_chunk=" << nt_chunk_ << ", axis=AXIS_FREQ, nbins=" << nbins << ", epsilon=" << epsilon << ")";

	this->name = ss.str();
	this->nt_chunk = nt_chunk_;
	this->nt_prepad = 0;
	this->nt_postpad = 0;
    }
    
    virtual void set_stream(const wi_stream &stream) override
    {
	this->nfreq = stream.nfreq;
	this->kernel = make_unique<rf_kernels::spline_detrender> (nfreq, nbins, epsilon);
    }

    virtual void process_chunk(double t0, double t1, float *intensity, float *weights, ssize_t stride, float *pp_intensity, float *pp_weights, ssize_t pp_stride) override
    {
	rf_assert(kernel.get() != nullptr);
	kernel->detrend(nt_chunk, stride, intensity, weights);
    }

    virtual void start_substream(int isubstream, double t0) override { }
    virtual void end_substream() override { }


    virtual Json::Value serialize_to_json() const override
    {
	Json::Value ret;

	ret["transform_name"] = "spline_detrender";
	ret["nt_chunk"] = int(this->nt_chunk);
	ret["axis"] = rf_kernels::axis_type_to_string(rf_kernels::AXIS_FREQ);
	ret["nbins"] = this->nbins;
	ret["epsilon"] = this->epsilon;

	return ret;
    }
};


// Externally callable factory function
shared_ptr<wi_transform> make_spline_detrender(int nt_chunk, rf_kernels::axis_type axis, int nbins, double epsilon)
{
    if (axis != rf_kernels::AXIS_FREQ)
	throw runtime_error("rf_pipelines::spline_detrender: only AXIS_FREQ is currently implemented");

    return make_shared<spline_detrender> (nt_chunk, nbins, epsilon);
}


}  // namespace rf_pipelines
