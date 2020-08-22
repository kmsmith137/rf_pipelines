#include <rf_kernels/upsample.hpp>
#include <rf_kernels/downsample.hpp>
#include "rf_pipelines_internals.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}  // namespace rf_pipelines
#endif


// FIXME: wi_sub_pipeline is currently implemented as a subclass of 'pipeline',
// but would it be cleaner to implement it as a subclass of 'pipeline_object'
// (or 'chunked_pipeline_object') with a member which is a shared_ptr<pipeline>?
//
// Note: if this is done, make sure to update the json-parsing in the web viewer
// and rfp-time.


// -------------------------------------------------------------------------------------------------
//
// Downsampler


struct downsampler : chunked_pipeline_object {
    downsampler(ssize_t Df, ssize_t Dt, ssize_t nt_chunk);

    // Inherits 'nt_chunk' from base class.
    const ssize_t Df;
    const ssize_t Dt;
    rf_kernels::wi_downsampler kernel;

    // Initialized in _bindc().
    shared_ptr<ring_buffer> rb_intensity_in;
    shared_ptr<ring_buffer> rb_weights_in;
    shared_ptr<ring_buffer> rb_intensity_out;
    shared_ptr<ring_buffer> rb_weights_out;
    ssize_t nfreq_out = 0;
    ssize_t nds_out = 0;
    ssize_t nt_out = 0;

    virtual void _bindc(ring_buffer_dict &rb_dict, Json::Value &json_attrs) override;
    virtual bool _process_chunk(ssize_t pos) override;
    virtual void _unbindc() override;
};


downsampler::downsampler(ssize_t Df_, ssize_t Dt_, ssize_t nt_chunk_) :
    chunked_pipeline_object("wi_sub_downsampler", false),  // can_be_first=false
    Df(Df_),
    Dt(Dt_),
    kernel(Df_, Dt_)
{ 
    this->nt_chunk = nt_chunk_;
    rf_assert(nt_chunk > 0);
}


void downsampler::_bindc(ring_buffer_dict &rb_dict, Json::Value &json_attrs)
{
    this->rb_intensity_in = get_buffer(rb_dict, "INTENSITY_HIRES");
    this->rb_weights_in = get_buffer(rb_dict, "WEIGHTS_HIRES");

    // Error checking is imperfect here, since the wi_sub_pipeline constructor
    // contains a complete set of checks.

    if (rb_intensity_in->cdims != rb_weights_in->cdims)
	_throw("'intensity' and 'weights' buffers have different dimensions");
    if (rb_intensity_in->nds != rb_weights_in->nds)
	_throw("'intensity' and 'weights' buffers have different downsampling");
    if (rb_intensity_in->cdims.size() != 1)
	_throw("expected intensity/weights arrays to be two-dimensional");

    this->nfreq_out = xdiv(rb_intensity_in->cdims[0], Df);
    this->nds_out = rb_intensity_in->nds * Dt;

    rf_assert(nt_chunk % (8*nds_out) == 0);
    this->nt_out = xdiv(nt_chunk, nds_out);

    this->rb_intensity_out = create_buffer(rb_dict, "INTENSITY", { nfreq_out }, nds_out);
    this->rb_weights_out = create_buffer(rb_dict, "WEIGHTS", { nfreq_out }, nds_out);
}


bool downsampler::_process_chunk(ssize_t pos)
{
    ring_buffer_subarray i_in(rb_intensity_in, pos, pos + nt_chunk, ring_buffer::ACCESS_READ);
    ring_buffer_subarray w_in(rb_weights_in, pos, pos + nt_chunk, ring_buffer::ACCESS_READ);
    ring_buffer_subarray i_out(rb_intensity_out, pos, pos + nt_chunk, ring_buffer::ACCESS_APPEND);
    ring_buffer_subarray w_out(rb_weights_out, pos, pos + nt_chunk, ring_buffer::ACCESS_APPEND);

    kernel.downsample(this->nfreq_out, this->nt_out, 
		      i_out.data, i_out.stride,
		      w_out.data, w_out.stride,
		      i_in.data, i_in.stride,
		      w_in.data, w_in.stride);

    return true;
}


void downsampler::_unbindc()
{
    this->rb_intensity_in.reset();
    this->rb_weights_in.reset();
    this->rb_intensity_out.reset();
    this->rb_weights_out.reset();
}


// -------------------------------------------------------------------------------------------------
//
// Upsampler


struct upsampler : chunked_pipeline_object {
    upsampler(ssize_t Df, ssize_t Dt, ssize_t nt_chunk, double w_cutoff);

    // Inherits 'nt_chunk' from base class.
    const ssize_t Df;
    const ssize_t Dt;
    const double w_cutoff;
    rf_kernels::weight_upsampler kernel;

    // Initialized in _bindc().
    shared_ptr<ring_buffer> rb_weights_in;
    shared_ptr<ring_buffer> rb_weights_out;
    ssize_t nfreq_in = 0;
    ssize_t nds_in = 0;
    ssize_t nt_in = 0;

    virtual void _bindc(ring_buffer_dict &rb_dict, Json::Value &json_attrs) override;
    virtual bool _process_chunk(ssize_t pos) override;
};


upsampler::upsampler(ssize_t Df_, ssize_t Dt_, ssize_t nt_chunk_, double w_cutoff_) :
    chunked_pipeline_object("wi_sub_upsampler", false),  // can_be_first=false
    Df(Df_),
    Dt(Dt_),
    w_cutoff(w_cutoff_),
    kernel(Df_, Dt_)
{ 
    this->nt_chunk = nt_chunk_;
    rf_assert(nt_chunk > 0);
}


void upsampler::_bindc(ring_buffer_dict &rb_dict, Json::Value &json_attrs)
{
    this->rb_weights_in = get_buffer(rb_dict, "WEIGHTS");
    this->rb_weights_out = get_buffer(rb_dict, "WEIGHTS_HIRES");

    if ((rb_weights_in->cdims.size() != 1) || (rb_weights_out->cdims.size() != 1))
	_throw("expected weights arrays to be two-dimensional");
    if (rb_weights_out->cdims[0] != (rb_weights_in->cdims[0] * Df))
	_throw("expected nfreq_hires == nfreq_lores * Df");
    if (rb_weights_in->nds != (rb_weights_out->nds * Dt))
	_throw("expected nds_lores == nds_hires * Dt");

    this->nfreq_in = rb_weights_in->cdims[0];
    this->nds_in = rb_weights_in->nds;

    rf_assert(nt_chunk % (8*nds_in) == 0);
    this->nt_in = xdiv(nt_chunk, nds_in);
}


bool upsampler::_process_chunk(ssize_t pos)
{
    ring_buffer_subarray w_in(rb_weights_in, pos, pos + nt_chunk, ring_buffer::ACCESS_READ);
    ring_buffer_subarray w_out(rb_weights_out, pos, pos + nt_chunk, ring_buffer::ACCESS_RW);

    kernel.upsample(nfreq_in, nt_in, w_out.data, w_out.stride, w_in.data, w_in.stride, w_cutoff);
    return true;
}


// -------------------------------------------------------------------------------------------------
//
// wi_sub_pipeline


wi_sub_pipeline::wi_sub_pipeline(const shared_ptr<pipeline_object> &sub_pipeline_, const initializer &ini_params_) :
    pipeline("wi_sub_pipeline", ""),
    ini_params(ini_params_),
    sub_pipeline(sub_pipeline_)
{
    // Initial ini_params sanity checking.
    if (ini_params.w_cutoff < 0.0)
	_throw("ini_params.w_cutoff cannot be negative");
    if (ini_params.nt_chunk < 0)
	_throw("ini_params.nt_chunk cannot be negative");
    if (ini_params.nfreq_out < 0)
	_throw("expected ini_params.nfreq_out >= 0");
    if (ini_params.nds_out < 0)
	_throw("expected ini_params.nds_out >= 0");
    if (ini_params.Df < 0)
	_throw("expected ini_params.Df >= 0");
    if (ini_params.Dt < 0)
	_throw("expected ini_params.Dt >= 0");
    if (ini_params.nfreq_out==0 && ini_params.Df==0)
	_throw("either ini_params.nfreq_out or ini_params.Df must be specified");	
    if (ini_params.nds_out==0 && ini_params.Dt==0)
	_throw("either ini_params.nds_out or ini_params.Dt must be specified");	
    if (ini_params.nds_out && ini_params.Dt && xmod(ini_params.nds_out, ini_params.Dt))
	_throw("ini_params.nds_out must be a multiple of ini_params.Dt, if both are specified");

    stringstream ss;
    ss << "wi_sub_pipeline(";
    if (ini_params.Df > 0)
	ss << "Df=" << ini_params.Df;
    if ((ini_params.Df > 0) && (ini_params.nfreq_out > 0))
	ss << ",";
    if (ini_params.nfreq_out > 0)
	ss << "nfreq_out=" << ini_params.nfreq_out;
    if (ini_params.Dt > 0)
	ss << ",Dt=" << ini_params.Dt;
    if (ini_params.nds_out > 0)
	ss << ",nds_out=" << ini_params.nds_out;
    ss << ")";

    this->name = ss.str();
}


void wi_sub_pipeline::_bind(ring_buffer_dict &rb_dict, Json::Value &json_attrs)
{
    if (!has_key(rb_dict, "INTENSITY"))
	_throw("buffer 'INTENSITY' does not exist in pipeline");
    if (!has_key(rb_dict, "WEIGHTS"))
	_throw("buffer 'WEIGHTS' does not exist in pipeline");

    // Don't call pipeline_object::get_buffer() since this has side effects.    
    auto rb_intensity = rb_dict["INTENSITY"];
    auto rb_weights = rb_dict["WEIGHTS"];

    if (rb_intensity->cdims != rb_weights->cdims)
	_throw("'intensity' and 'weights' buffers have different dimensions");
    if (rb_intensity->nds != rb_weights->nds)
	_throw("'intensity' and 'weights' buffers have different downsampling");
    if (rb_intensity->cdims.size() != 1)
	_throw("expected intensity/weights arrays to be two-dimensional");

    ssize_t nfreq_in = rb_intensity->cdims[0];
    ssize_t nds_in = rb_intensity->nds;

    if ((ini_params.Df != 0) && (ini_params.nfreq_out != 0) && (nfreq_in != ini_params.nfreq_out * ini_params.Df))
	_throw("nfreq_in (" + to_string(nfreq_in) + ") does not match expected value (" + to_string(ini_params.nfreq_out * ini_params.Df) + ")");

    if ((ini_params.Dt != 0) && (ini_params.nds_out != 0) && (nds_in != xdiv(ini_params.nds_out, ini_params.Dt)))
	_throw("nds_in (" + to_string(nds_in) + ") does not match expected value (" + to_string(xdiv(ini_params.nds_out, ini_params.Dt)));

    // Compute (Df,Dt)
    ssize_t Df = ini_params.Df;  // can be zero
    ssize_t Dt = ini_params.Dt;  // can be zero
    
    if (Df == 0) {
	// Note: wi_sub_pipeline constructor has already checked (ini_params.nfreq_out > 0).
	if (xmod(nfreq_in, ini_params.nfreq_out))
	    _throw("nfreq_in (" + to_string(nfreq_in) + ") is not a multiple of ini_params.nfreq_out (" + to_string(ini_params.nfreq_out) + ")");
	Df = xdiv(nfreq_in, ini_params.nfreq_out);
    }

    if (Dt == 0) {
	// Note: wi_sub_pipeline constructor has already checked (ini_params.nds_out > 0).
	if (xmod(ini_params.nds_out, nds_in))
	    _throw("nds_in (" + to_string(nds_in) + ") does not divide ini_params.nds_out (" + to_string(nds_in) + ")");
	Dt = xdiv(ini_params.nds_out, nds_in);
    }

    // Choose nt_chunk for downsampler/upsampler.
    // The logic here is similar to chunked_pipeline_object::finalize_nt_chunk().

    ssize_t nt_chunk = ini_params.nt_chunk;       // can be zero
    ssize_t min_nt_chunk = 8 * Dt * nds_in;  // nt_chunk must be a multiple of this

    if (nt_chunk == 0)
	nt_chunk = max(nt_chunk_in/min_nt_chunk, ssize_t(1)) * min_nt_chunk;
    else if ((ini_params.nt_chunk % min_nt_chunk) != 0)
	_throw("ini_params.nt_chunk (" + to_string(ini_params.nt_chunk) + " must be a multiple of " + to_string(min_nt_chunk) + " in this pipeline");

    // Note: can't call pipeline::add() directly (get error message "...add() was called after bind()")
    rf_assert(this->elements.size() == 0);
    this->elements.push_back(make_shared<downsampler> (Df, Dt, nt_chunk));
    this->elements.push_back(this->sub_pipeline);
    this->elements.push_back(make_shared<upsampler> (Df, Dt, nt_chunk, ini_params.w_cutoff));

    ring_buffer_dict rb_dict2;    
    rb_dict2["INTENSITY_HIRES"] = rb_dict["INTENSITY"];
    rb_dict2["WEIGHTS_HIRES"] = rb_dict["WEIGHTS"];

    pipeline::_bind(rb_dict2, json_attrs);
}


// virtual override
void wi_sub_pipeline::_unbind()
{
    pipeline::_unbind();
    this->elements.clear();
}


// virutal override
ssize_t wi_sub_pipeline::get_preferred_chunk_size()
{
    return 0;
}


// virtual override
void wi_sub_pipeline::_visit_pipeline(std::function<void(const std::shared_ptr<pipeline_object>&,int)> f, const std::shared_ptr<pipeline_object> &self, int depth)
{
    // Overriding pipeline::_visit_pipeline() is necessary here, for two reasons.
    //
    // First, pipeline::_visit_pipeline() would also visit the "internal" downsampler/upsampler,
    // and I thought it made most sense to "hide" this implementation detail from the caller of
    // visit_pipeline().  (This is a design decision that could be reversed.)
    //
    // Second, we want _visit_pipeline() to work before bind() is called, and in the current
    // implementation, pipeline::elements is an empty vector before bind() is called.

    rf_assert(self.get() == this);
    f(self, depth);
    
    visit_pipeline(f, sub_pipeline, depth+1);
}


// virtual override
Json::Value wi_sub_pipeline::jsonize() const 
{
    Json::Value ret;

    ret["class_name"] = "wi_sub_pipeline";
    ret["sub_pipeline"] = sub_pipeline->jsonize();
    ret["w_cutoff"] = ini_params.w_cutoff;
    ret["nfreq_out"] = Json::Int64(ini_params.nfreq_out);
    ret["nds_out"] = Json::Int64(ini_params.nds_out);
    ret["Df"] = Json::Int64(ini_params.Df);
    ret["Dt"] = Json::Int64(ini_params.Dt);

    return ret;
}


// static 
shared_ptr<wi_sub_pipeline> wi_sub_pipeline::from_json(const Json::Value &j)
{
    wi_sub_pipeline::initializer ini_params;
    ini_params.w_cutoff = double_from_json(j, "w_cutoff");
    ini_params.nfreq_out = int_from_json(j, "nfreq_out");
    ini_params.nds_out = int_from_json(j, "nds_out");
    ini_params.Df = int_from_json(j, "Df");
    ini_params.Dt = int_from_json(j, "Dt");

    if (!j.isMember("sub_pipeline"))
	throw runtime_error("rf_pipelines::wi_sub_pipeline::from_json(): json member 'sub_pipeline' does not exist");

    shared_ptr<pipeline_object> sub_pipeline = pipeline_object::from_json(j["sub_pipeline"]);

    return make_shared<wi_sub_pipeline> (sub_pipeline, ini_params);
}


namespace {
    struct _init {
	_init() {
	    pipeline_object::register_json_deserializer("wi_sub_pipeline", wi_sub_pipeline::from_json);
	}
    } init;
}


}  // namespace rf_pipelines
