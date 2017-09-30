#include "rf_pipelines_internals.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}  // namespace rf_pipelines
#endif


// FIXME: the downsampler and upsampler use a default nt_chunk, which may not end up
// satisfying the necessary divisibility conditions (must be a multiple of nds_out * S,
// where S=8).


// -------------------------------------------------------------------------------------------------
//
// Downsampler


struct downsampler : chunked_pipeline_object {
    wi_sub_pipeline::initializer ini_params;
    ssize_t nfreq_in = 0;
    ssize_t nds_in = 0;

    downsampler(const wi_sub_pipeline::initializer &ini_params);

    shared_ptr<ring_buffer> rb_intensity_in;
    shared_ptr<ring_buffer> rb_weights_in;

    shared_ptr<ring_buffer> rb_intensity_out;
    shared_ptr<ring_buffer> rb_weights_out;

    virtual void _bind_chunked(ring_buffer_dict &rb_dict, Json::Value &json_data) override;
    virtual bool _process_chunk(ssize_t pos) override;
};


downsampler::downsampler(const wi_sub_pipeline::initializer &ini_params_) :
    chunked_pipeline_object("wi_sub_downsampler", false),  // can_be_first=false
    ini_params(ini_params_)
{ }


void downsampler::_bind_chunked(ring_buffer_dict &rb_dict, Json::Value &json_data)
{
    // Get input buffers, sanity check, compute nfreq_in and nds_in.

    this->rb_intensity_in = get_buffer(rb_dict, "INTENSITY_HIRES");
    this->rb_weights_in = get_buffer(rb_dict, "WEIGHTS_HIRES");

    if (rb_intensity_in->cdims != rb_weights_in->cdims)
	_throw("'intensity' and 'weights' buffers have different dimensions");
    if (rb_intensity_in->nds != rb_weights_in->nds)
	_throw("'intensity' and 'weights' buffers have different downsampling");
    if (rb_intensity_in->cdims.size() != 1)
	_throw("expected intensity/weights arrays to be two-dimensional");

    this->nfreq_in = rb_intensity_in->cdims[0];
    this->nds_in = rb_intensity_in->nds;

    // More sanity checking, compute all ini_params.

    if (ini_params.nfreq_out == 0) {
	if (xmod(nfreq_in, ini_params.Df))
	    _throw("nfreq_in (" + to_string(nfreq_in) + ") is not a multiple of ini_params.Df (" + to_string(ini_params.Df) + ")");
	ini_params.nfreq_out = xdiv(nfreq_in, ini_params.Df);
    }
    
    if (ini_params.Df == 0) {
	if (xmod(nfreq_in, ini_params.nfreq_out))
	    _throw("nfreq_out (" + to_string(nfreq_in) + ") is not a multiple of ini_params.nfreq_out (" + to_string(ini_params.nfreq_out) + ")");
	ini_params.Df = xdiv(nfreq_in, ini_params.nfreq_out);
    }

    if (nfreq_in != ini_params.nfreq_out * ini_params.Df)
	_throw("nfreq_in (" + to_string(nfreq_in) + ") does not match expected value (" + to_string(ini_params.nfreq_out * ini_params.Df) + ")");

    if (ini_params.nds_out == 0)
	ini_params.nds_out = nds_in * ini_params.Dt;

    if (ini_params.Dt == 0) {
	if (xmod(ini_params.nds_out, nds_in))
	    _throw("nds_in (" + to_string(nds_in) + ") does not divide ini_params.nds_out (" + to_string(nds_in) + ")");
	ini_params.Dt = xdiv(ini_params.nds_out, nds_in);
    }

    if (nds_in != xdiv(ini_params.nds_out, ini_params.Dt))
	_throw("nds_in (" + to_string(nds_in) + ") does not match expected value (" + to_string(xdiv(ini_params.nds_out, ini_params.Dt)) + ")");

    // Create input buffers.

    this->rb_intensity_out = create_buffer(rb_dict, "INTENSITY", { ini_params.nfreq_out }, ini_params.nds_out);
    this->rb_weights_out = create_buffer(rb_dict, "WEIGHTS", { ini_params.nfreq_out }, ini_params.nds_out);

    // FIXME: construct kernel object here.
}


bool downsampler::_process_chunk(ssize_t pos)
{
    float *int_in = rb_intensity_in->get(pos, pos + nt_chunk, ring_buffer::ACCESS_READ);
    float *wt_in = rb_weights_in->get(pos, pos + nt_chunk, ring_buffer::ACCESS_READ);
    float *int_out = rb_intensity_out->get(pos, pos + nt_chunk, ring_buffer::ACCESS_APPEND);
    float *wt_out = rb_weights_out->get(pos, pos + nt_chunk, ring_buffer::ACCESS_APPEND);
    
    // FIXME: run downsampling kernel here
    
    rb_intensity_in->put(int_in, pos, pos + nt_chunk, ring_buffer::ACCESS_READ);
    rb_weights_in->put(wt_in, pos, pos + nt_chunk, ring_buffer::ACCESS_READ);
    rb_intensity_out->put(int_out, pos, pos + nt_chunk, ring_buffer::ACCESS_APPEND);
    rb_weights_out->put(wt_out, pos, pos + nt_chunk, ring_buffer::ACCESS_APPEND);

    return true;
}


// -------------------------------------------------------------------------------------------------
//
// Upsampler


struct upsampler : chunked_pipeline_object {
    upsampler();

    shared_ptr<ring_buffer> rb_weights_in;
    shared_ptr<ring_buffer> rb_weights_out;

    virtual void _bind_chunked(ring_buffer_dict &rb_dict, Json::Value &json_data) override;
    virtual bool _process_chunk(ssize_t pos) override;
};


upsampler::upsampler() :
    chunked_pipeline_object("wi_sub_upsampler", false)  // can_be_first=false
{ }


void upsampler::_bind_chunked(ring_buffer_dict &rb_dict, Json::Value &json_data)
{
    this->rb_weights_in = get_buffer(rb_dict, "WEIGHTS");
    this->rb_weights_out = get_buffer(rb_dict, "WEIGHTS_HIRES");

    // FIXME: create kernel object
}


bool upsampler::_process_chunk(ssize_t pos)
{
    float *wt_in = rb_weights_in->get(pos, pos + nt_chunk, ring_buffer::ACCESS_READ);
    float *wt_out = rb_weights_out->get(pos, pos + nt_chunk, ring_buffer::ACCESS_RW);

    // FIXME run upsampling krnel

    rb_weights_in->put(wt_in, pos, pos + nt_chunk, ring_buffer::ACCESS_READ);
    rb_weights_out->put(wt_out, pos, pos + nt_chunk, ring_buffer::ACCESS_RW);    
    
    return true;
}


// -------------------------------------------------------------------------------------------------
//
// wi_sub_pipeline


wi_sub_pipeline::wi_sub_pipeline(const shared_ptr<pipeline_object> &sub_pipeline_, const initializer &ini_params_) :
    pipeline("wi_sub_pipeline"),
    ini_params(ini_params_),
    sub_pipeline(sub_pipeline_)
{
    // A little sanity-checking on ini_params.
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

    // Call pipeline::add() to make the 3-element pipeline.
    this->add(make_shared<downsampler> (ini_params));
    this->add(sub_pipeline);
    this->add(make_shared<upsampler> ());
}


void wi_sub_pipeline::_bind(ring_buffer_dict &rb_dict, Json::Value &json_data)
{
    if (!has_key(rb_dict, "INTENSITY"))
	_throw("buffer 'INTENSITY' does not exist in pipeline");
    if (!has_key(rb_dict, "WEIGHTS"))
	_throw("buffer 'WEIGHTS' does not exist in pipeline");

    ring_buffer_dict rb_dict2;    
    rb_dict2["INTENSITY_HIRES"] = rb_dict["INTENSITY"];
    rb_dict2["WEIGHTS_HIRES"] = rb_dict["WEIGHTS"];

    pipeline::_bind(rb_dict2, json_data);
}


// virtual override
Json::Value wi_sub_pipeline::jsonize() const 
{
    Json::Value ret;

    ret["class_name"] = "wi_sub_pipeline";
    ret["sub_pipeline"] = sub_pipeline->jsonize();
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
	    pipeline_object::register_json_constructor("wi_sub_pipeline", wi_sub_pipeline::from_json);
	}
    } init;
}


}  // namespace rf_pipelines
