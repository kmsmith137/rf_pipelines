#include "rf_kernels/spline_detrender.hpp"
#include "rf_pipelines_internals.hpp"
//#include "ch_frb_io/chlog.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif

spline_detrender::spline_detrender(int nt_chunk_, rf_kernels::axis_type axis_, int nbins_, double epsilon_) :
    wi_transform("spline_detrender"),
    nbins(nbins_),
    epsilon(epsilon_),
    axis(axis_)
{
    // Temporary
    if (axis != rf_kernels::AXIS_FREQ)
        throw runtime_error("rf_pipelines::spline_detrender: only AXIS_FREQ is currently implemented");

    // Superfluous for now, but will make sense when AXIS_TIME and/or AXIS_NONE are implemented.
    if ((nt_chunk_ == 0) && (axis != rf_kernels::AXIS_FREQ))
        throw runtime_error("rf_pipelines::spline_detrender: nt_chunk must be specified (unless axis=AXIS_FREQ)");

    stringstream ss;
    ss << "spline_detrender(nt_chunk=" << nt_chunk_ << ", axis=" << rf_kernels::axis_type_to_string(axis) << ", nbins=" << nbins << ", epsilon=" << epsilon << ")";

    this->name = ss.str();
    this->nt_chunk = nt_chunk_;
    this->kernel_chunk_size = 8;
    this->nds = 0;  // allows spline_detrender to run inside a wi_sub_pipeline.
}

// Called after this->nfreq is initialized.
void spline_detrender::_bind_transform(Json::Value &json_attrs)
{
    this->kernel = make_unique<rf_kernels::spline_detrender> (nfreq, nbins, epsilon);
}

void spline_detrender::_bind_transform_rb(ring_buffer_dict &rb_dict) {
    if (this->ringbuf_nhistory) {
        cout << "Spline_detrender: allocating a ring buffer to store " << this->ringbuf_nhistory << " chunks of history!  nt_chunk " << nt_chunk << ", nds " << nds << ", nbins " << nbins << endl;
    }
}

void spline_detrender::_process_chunk(float *intensity, ssize_t istride, float *weights, ssize_t wstride, ssize_t pos)
{
    // Note xdiv(nt_chunk, nds) here.
    rf_assert(kernel.get() != nullptr);
    kernel->detrend(xdiv(nt_chunk,nds), intensity, istride, weights, wstride);
}

void spline_detrender::_unbind_transform()
{
    this->kernel.reset();
}

Json::Value spline_detrender::jsonize() const
{
    Json::Value ret;

    ret["class_name"] = "spline_detrender";
    ret["nt_chunk"] = int(this->get_prebind_nt_chunk());
    ret["axis"] = rf_kernels::axis_type_to_string(rf_kernels::AXIS_FREQ);
    ret["nbins"] = this->nbins;
    ret["epsilon"] = this->epsilon;

    return ret;
}

shared_ptr<spline_detrender> spline_detrender::from_json(const Json::Value &j)
{
    int nbins = int_from_json(j, "nbins");
    ssize_t nt_chunk = ssize_t_from_json(j, "nt_chunk");
    double epsilon = double_from_json(j, "epsilon");
    rf_kernels::axis_type axis = axis_type_from_json(j, "axis");

    return make_shared<spline_detrender> (nt_chunk, axis, nbins, epsilon);
}

void spline_detrender::set_ringbuffer_size(int nhistory) {
    if (this->state != UNBOUND)
	throw runtime_error("spline_detrender::set_ringbuffer_size() called after bind()");
    if (nhistory < 0)
	throw runtime_error("spline_detrender::set_ringbuffer_size(): nhistory was negative");
    this->ringbuf_nhistory = nhistory;
}

namespace {
    struct _init {
	_init() {
	    pipeline_object::register_json_deserializer("spline_detrender", spline_detrender::from_json);
	}
    } init;
}


// Externally callable factory function
shared_ptr<wi_transform> make_spline_detrender(int nt_chunk, rf_kernels::axis_type axis, int nbins, double epsilon)
{
    return make_shared<spline_detrender> (nt_chunk, axis, nbins, epsilon);
}


}  // namespace rf_pipelines
