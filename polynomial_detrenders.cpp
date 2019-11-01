#include "rf_pipelines_internals.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif


polynomial_detrender::polynomial_detrender(rf_kernels::axis_type axis, int nt_chunk_, int polydeg, double epsilon_) :
    wi_transform("polynomial_detrender"),
    kernel(axis, polydeg),
    epsilon(epsilon_)
{
    stringstream ss;
    ss << "polynomial_detrender(nt_chunk=" << nt_chunk_ << ", axis=" << axis << ", polydeg=" << polydeg << ", epsilon=" << epsilon_ << ")";

    this->name = ss.str();
    this->nt_chunk = nt_chunk_;
    this->kernel_chunk_size = 8;
    this->nds = 0;  // allows polynomial_detrender to run in a wi_sub_pipeline.
	
    if ((nt_chunk == 0) && (axis != rf_kernels::AXIS_FREQ))
        throw runtime_error("rf_pipelines::polynomial_detrender: nt_chunk must be specified (unless axis=AXIS_FREQ)");
}

void polynomial_detrender::_allocate() {
    if (this->kernel.axis == rf_kernels::AXIS_TIME)
        this->coeffs = vector<float>(this->nfreq * (this->kernel.polydeg + 1));
    else if (this->kernel.axis == rf_kernels::AXIS_FREQ)
        this->coeffs = vector<float>(this->nt_chunk * (this->kernel.polydeg + 1));
    else
        throw runtime_error("rf_pipelines::polynomial_detrender: neither TIME nor FREQ axis!");
}

void polynomial_detrender::_deallocate() {
    this->coeffs = vector<float>();
}

void polynomial_detrender::_process_chunk(float *intensity, ssize_t istride, float *weights, ssize_t wstride, ssize_t pos) {
    // Note xdiv(nt_chunk,nds) here
    this->kernel.detrend(nfreq, xdiv(nt_chunk,nds), intensity, istride, weights, wstride, epsilon, coeffs.data());
    this->_handle_coeffs(pos, coeffs.data(), intensity, istride, weights, wstride);
}

void polynomial_detrender::_handle_coeffs(ssize_t pos, float *coeffs, float *intensity, ssize_t istride, float *weights, ssize_t wstride) {
}

Json::Value polynomial_detrender::jsonize() const {
    Json::Value ret;

    ret["class_name"] = "polynomial_detrender";
    ret["nt_chunk"] = int(this->get_prebind_nt_chunk());
    ret["axis"] = rf_kernels::axis_type_to_string(this->kernel.axis);
    ret["polydeg"] = this->kernel.polydeg;
    ret["epsilon"] = this->epsilon;

    return ret;
}

shared_ptr<polynomial_detrender> polynomial_detrender::from_json(const Json::Value &x) {
    if (string_from_json(x,"class_name") != "polynomial_detrender")
        throw runtime_error("rf_pipelines: expected class_name=\"pipeline\" in pipeline json constructor");

    rf_kernels::axis_type axis = axis_type_from_json(x, "axis");
    int nt_chunk = int_from_json(x, "nt_chunk");
    double polydeg = double_from_json(x, "polydeg");
    double epsilon = double_from_json(x, "epsilon");
	
    return make_shared<polynomial_detrender> (axis, nt_chunk, polydeg, epsilon);
}

// Externally callable factory function
shared_ptr<wi_transform> make_polynomial_detrender(int nt_chunk, rf_kernels::axis_type axis, int polydeg, double epsilon)
{
    return make_shared<polynomial_detrender> (axis, nt_chunk, polydeg, epsilon);
}


namespace {
    struct _init {
	_init() {
	    pipeline_object::register_json_deserializer("polynomial_detrender", polynomial_detrender::from_json);
	}
    } init;
}


}  // namespace rf_pipelines
