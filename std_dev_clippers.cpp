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

    std_dev_clipper_transform(int Df_, int Dt_, rf_kernels::axis_type axis_, int nt_chunk_, double sigma_, bool two_pass_) :
	wi_transform("std_dev_clipper", nt_chunk_),
	Df(Df_),
	Dt(Dt_),
	axis(axis_),
	two_pass(two_pass_),
	sigma(sigma_)
    {	
	stringstream ss;
        ss << "std_dev_clipper(nt_chunk=" << nt_chunk_ << ", axis=" << axis << ", sigma=" << sigma
           << ", Df=" << Df << ", Dt=" << Dt << ", two_pass=" << two_pass << ")";
	
        this->name = ss.str();
	
	if (nt_chunk == 0)
	    throw runtime_error("rf_pipelines::std_dev_clipper: nt_chunk must be specified");
	
	// Can't construct the kernel yet, since 'nfreq' is not known until set_stream()
	// However, for argument checking purposes, we construct a dummy kernel with nfreq=max(Df,8).
	// FIXME eventaully there will be a constructor argument 'allocate=false' that will make sense here.

	int nfreq_dummy = max(Df,8);
	rf_kernels::std_dev_clipper dummy(nfreq_dummy, nt_chunk, axis, sigma, Df, Dt, two_pass);
    }

    virtual ~std_dev_clipper_transform() { }

    // Called after (nfreq, nds) are initialized.
    virtual void _bind_transform(Json::Value &json_data) override
    {
	if (nfreq % Df)
	    throw runtime_error("rf_pipelines std_dev_clipper: nfreq (=" + to_string(nfreq) 
				+ ") is not divisible by frequency downsampling factor Df=" + to_string(Df));
	
	this->kernel = make_unique<rf_kernels::std_dev_clipper> (nfreq, xdiv(nt_chunk,nds), axis, sigma, Df, Dt, two_pass);
    }

    virtual void _process_chunk(float *intensity, ssize_t istride, float *weights, ssize_t wstride, ssize_t pos) override
    {
	this->kernel->clip(intensity, istride, weights, wstride);
    }

    virtual Json::Value jsonize() const override
    {
	Json::Value ret;

	ret["class_name"] = "std_dev_clipper";
	ret["Df"] = Df;
	ret["Dt"] = Dt;
	ret["sigma"] = sigma;
	ret["two_pass"] = two_pass;
	ret["nt_chunk"] = int(this->get_orig_nt_chunk());
	ret["axis"] = rf_kernels::axis_type_to_string(axis);
	
	return ret;
    }

    static shared_ptr<std_dev_clipper_transform> from_json(const Json::Value &j)
    {
	int Df = int_from_json(j, "Df");
	int Dt = int_from_json(j, "Dt");
	bool two_pass = bool_from_json(j, "two_pass");
	double sigma = double_from_json(j, "sigma");
	ssize_t nt_chunk = int_from_json(j, "nt_chunk");
	rf_kernels::axis_type axis = axis_type_from_json(j, "axis");

	return make_shared<std_dev_clipper_transform> (Df, Dt, axis, nt_chunk, sigma, two_pass);	
    }
};


namespace {
    struct _init {
	_init() {
	    pipeline_object::register_json_constructor("std_dev_clipper", std_dev_clipper_transform::from_json);
	}
    } init;
}


// Externally callable
shared_ptr<wi_transform> make_std_dev_clipper(int nt_chunk, rf_kernels::axis_type axis, double sigma, int Df, int Dt, bool two_pass)
{
    return make_shared<std_dev_clipper_transform> (Df, Dt, axis, nt_chunk, sigma, two_pass);
}


}  // namespace rf_pipelines
