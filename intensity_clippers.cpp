#include <rf_kernels/mean_rms.hpp>
#include <rf_kernels/intensity_clipper.hpp>

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
	: wi_transform("intensity_clipper"),
	  Df(Df_),
	  Dt(Dt_),
	  axis(axis_),
	  niter(niter_),
	  sigma(sigma_),
	  iter_sigma(iter_sigma_ ? iter_sigma_ : sigma_),
	  two_pass(two_pass_)
    {
	stringstream ss;
        ss << "intensity_clipper(nt_chunk=" << nt_chunk_ << ", axis=" << axis 
	   << ", sigma=" << sigma << ", niter=" << niter << ", iter_sigma=" << iter_sigma 
	   << ", Df=" << Df << ", Dt=" << Dt << ", two_pass=" << two_pass << ")";

	this->name = ss.str();
	this->nt_chunk = nt_chunk_;
	this->kernel_chunk_size = 8 * Dt;
	this->nds = 0;   // allows intensity_clipper to run in a wi_sub_pipeline.
	
	if ((nt_chunk == 0) && (axis != rf_kernels::AXIS_FREQ))
	    throw runtime_error("rf_pipelines::intensity_clipper: nt_chunk must be specified (unless axis=AXIS_FREQ)");

	// Can't construct the kernel yet, since 'nfreq' and 'nds' are not known until bind().
	// However, for argument checking purposes, we construct a dummy kernel with (nfreq,nt_chunk)=(Df,8*Dt).
	// FIXME eventually there will be a constructor argument 'allocate=false' that will make sense here.
	
	rf_kernels::intensity_clipper dummy(Df, 8*Dt, axis, sigma, Df, Dt, niter, iter_sigma, two_pass);
    }

    virtual ~intensity_clipper_transform() { }

    // Called after (nfreq, nds) are initialized.
    virtual void _bind_transform(Json::Value &json_attrs) override
    {
	if (nfreq % Df)
	    throw runtime_error("rf_pipelines::intensity_clipper: nfreq (=" + to_string(nfreq) 
				+ ") is not divisible by frequency downsampling factor Df=" + to_string(Df));

	// Note xdiv(nt_chunk, nds) here.
	this->kernel = make_unique<rf_kernels::intensity_clipper> (nfreq, xdiv(nt_chunk,nds), axis, sigma, Df, Dt, niter, iter_sigma, two_pass);
    }

    virtual void _process_chunk(float *intensity, ssize_t istride, float *weights, ssize_t wstride, ssize_t pos) override
    {
	if (pos >= 10000)
	    throw runtime_error("DOH");
	this->kernel->clip(intensity, istride, weights, wstride);
    }

    virtual Json::Value jsonize() const override
    {
	Json::Value ret;

	ret["class_name"] = "intensity_clipper";
	ret["Df"] = Df;
	ret["Dt"] = Dt;
	ret["axis"] = rf_kernels::axis_type_to_string(axis);
	ret["nt_chunk"] = int(this->get_prebind_nt_chunk());
	ret["sigma"] = sigma;
	ret["niter"] = niter;
	ret["iter_sigma"] = iter_sigma;
	ret["two_pass"] = two_pass;

	return ret;
    }

    virtual void _unbind_transform() override
    {
	this->kernel.reset();
    }

    static shared_ptr<intensity_clipper_transform> from_json(const Json::Value &j)
    {
	int Df = int_from_json(j, "Df");
	int Dt = int_from_json(j, "Dt");
	int niter = int_from_json(j, "niter");
	int nt_chunk = int_from_json(j, "nt_chunk");
	bool two_pass = bool_from_json(j, "two_pass");
	double sigma = double_from_json(j, "sigma");
	double iter_sigma = double_from_json(j, "iter_sigma");
	rf_kernels::axis_type axis = axis_type_from_json(j, "axis");

	return make_shared<intensity_clipper_transform> (Df, Dt, axis, nt_chunk, sigma, niter, iter_sigma, two_pass);
    }
};


namespace {
    struct _init {
	_init() {
	    pipeline_object::register_json_deserializer("intensity_clipper", intensity_clipper_transform::from_json);
	}
    } init;
}


// -------------------------------------------------------------------------------------------------


// Externally visible
shared_ptr<wi_transform> make_intensity_clipper(int nt_chunk, rf_kernels::axis_type axis, double sigma, int niter, double iter_sigma, int Df, int Dt, bool two_pass)
{
    return make_shared<intensity_clipper_transform> (Df, Dt, axis, nt_chunk, sigma, niter, iter_sigma, two_pass);
}



}  // namespace rf_pipelines
