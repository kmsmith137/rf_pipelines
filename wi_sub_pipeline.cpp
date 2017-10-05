#include <rf_kernels/upsample.hpp>
#include <rf_kernels/downsample.hpp>
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
    downsampler(ssize_t Df, ssize_t Dt, ssize_t nt_chunk);

    // Inherits 'nt_chunk' from base class.
    const ssize_t Df;
    const ssize_t Dt;
    rf_kernels::wi_downsampler kernel;

    // Initialized in _bind_chunked().
    shared_ptr<ring_buffer> rb_intensity_in;
    shared_ptr<ring_buffer> rb_weights_in;
    shared_ptr<ring_buffer> rb_intensity_out;
    shared_ptr<ring_buffer> rb_weights_out;
    ssize_t nfreq_out = 0;
    ssize_t nds_out = 0;
    ssize_t nt_out = 0;

    virtual void _bind_chunked(ring_buffer_dict &rb_dict, Json::Value &json_attrs) override;
    virtual bool _process_chunk(ssize_t pos) override;
};


downsampler::downsampler(ssize_t Df_, ssize_t Dt_, ssize_t nt_chunk_) :
    chunked_pipeline_object("wi_sub_downsampler", false, nt_chunk_),  // can_be_first=false
    Df(Df_),
    Dt(Dt_),
    kernel(Df_, Dt_)
{ }


void downsampler::_bind_chunked(ring_buffer_dict &rb_dict, Json::Value &json_attrs)
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
    this->nt_out = xdiv(nt_chunk, nds_out);

    this->rb_intensity_out = create_buffer(rb_dict, "INTENSITY", { nfreq_out }, nds_out);
    this->rb_weights_out = create_buffer(rb_dict, "WEIGHTS", { nfreq_out }, nds_out);
}


bool downsampler::_process_chunk(ssize_t pos)
{
    float *int_in = rb_intensity_in->get(pos, pos + nt_chunk, ring_buffer::ACCESS_READ);
    float *wt_in = rb_weights_in->get(pos, pos + nt_chunk, ring_buffer::ACCESS_READ);
    float *int_out = rb_intensity_out->get(pos, pos + nt_chunk, ring_buffer::ACCESS_APPEND);
    float *wt_out = rb_weights_out->get(pos, pos + nt_chunk, ring_buffer::ACCESS_APPEND);

    kernel.downsample(this->nfreq_out, this->nt_out,
		      int_out, rb_intensity_out->get_stride(), 
		      wt_out, rb_weights_out->get_stride(),
		      int_in, rb_intensity_in->get_stride(),
		      wt_in, rb_weights_in->get_stride());
    
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
    upsampler(ssize_t Df, ssize_t Dt, ssize_t nt_chunk, double w_cutoff);

    // Inherits 'nt_chunk' from base class.
    const ssize_t Df;
    const ssize_t Dt;
    const double w_cutoff;
    rf_kernels::weight_upsampler kernel;

    // Initialized in _bind_chunked().
    shared_ptr<ring_buffer> rb_weights_in;
    shared_ptr<ring_buffer> rb_weights_out;
    ssize_t nfreq_in = 0;
    ssize_t nt_in = 0;

    virtual void _bind_chunked(ring_buffer_dict &rb_dict, Json::Value &json_attrs) override;
    virtual bool _process_chunk(ssize_t pos) override;
};


upsampler::upsampler(ssize_t Df_, ssize_t Dt_, ssize_t nt_chunk_, double w_cutoff_) :
    chunked_pipeline_object("wi_sub_upsampler", false, nt_chunk_),  // can_be_first=false
    Df(Df_),
    Dt(Dt_),
    w_cutoff(w_cutoff_),
    kernel(Df_, Dt_)
{ }


void upsampler::_bind_chunked(ring_buffer_dict &rb_dict, Json::Value &json_attrs)
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
    this->nt_in = xdiv(nt_chunk, rb_weights_in->nds);
}


bool upsampler::_process_chunk(ssize_t pos)
{
    float *wt_in = rb_weights_in->get(pos, pos + nt_chunk, ring_buffer::ACCESS_READ);
    float *wt_out = rb_weights_out->get(pos, pos + nt_chunk, ring_buffer::ACCESS_RW);

    kernel.upsample(nfreq_in, nt_in, 
		    wt_out, rb_weights_out->get_stride(), 
		    wt_in, rb_weights_in->get_stride(), 
		    w_cutoff);

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
    ssize_t nt_chunk = ini_params.nt_chunk;  // can be zero
    ssize_t min_nt_chunk = 8 * Dt * nds_in;  // nt_chunk must be a multiple of this

    if (nt_chunk == 0) {
	nt_chunk = (nt_chunk_in / min_nt_chunk) * min_nt_chunk;
	nt_chunk = max(nt_chunk, min_nt_chunk);
    }
    else if ((ini_params.nt_chunk % min_nt_chunk) != 0)
	_throw("ini_params.nt_chunk (" + to_string(ini_params.nt_chunk) + " must be a multiple of " + to_string(min_nt_chunk) + " in this pipeline");

    // Note: can't call pipeline::add() directly (get error message "...add() was called after bind()")
    rf_assert(this->elements.size() == 0);
    this->elements.push_back(make_shared<downsampler> (Df, Dt, nt_chunk));
    this->elements.push_back(sub_pipeline);
    this->elements.push_back(make_shared<upsampler> (Df, Dt, nt_chunk, ini_params.w_cutoff));
    this->_update_name();

    ring_buffer_dict rb_dict2;    
    rb_dict2["INTENSITY_HIRES"] = rb_dict["INTENSITY"];
    rb_dict2["WEIGHTS_HIRES"] = rb_dict["WEIGHTS"];

    pipeline::_bind(rb_dict2, json_attrs);
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
	    pipeline_object::register_json_constructor("wi_sub_pipeline", wi_sub_pipeline::from_json);
	}
    } init;
}


}  // namespace rf_pipelines
