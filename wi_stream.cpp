#include "rf_pipelines_internals.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}  // namespace rf_pipelines
#endif


wi_stream::wi_stream(const string &name, ssize_t nfreq_, ssize_t nt_chunk_) :
    chunked_pipeline_object(name, true, nt_chunk_),   // can_be_first=true
    nfreq(nfreq_)
{
    if (nfreq < 0)
	_throw("expected nfreq >= 0 in wi_stream constructor");
}


void wi_stream::_bind_chunked(ring_buffer_dict &rb_dict, Json::Value &json_attrs)
{
    // Optional subclass-specific initializations.
    this->_bind_stream(json_attrs);

    if (nfreq == 0)
	_throw("nfreq must be initialized in either wi_stream constructor, or in wi_stream::_bind_stream()");

    this->rb_intensity = this->create_buffer(rb_dict, "INTENSITY", {nfreq}, 1);
    this->rb_weights = this->create_buffer(rb_dict, "WEIGHTS", {nfreq}, 1);

    json_attrs["nfreq"] = Json::Int64(nfreq);
}


bool wi_stream::_process_chunk(ssize_t pos)
{
    float *intensity = rb_intensity->get(pos, pos+nt_chunk, ring_buffer::ACCESS_APPEND);
    float *weights = rb_weights->get(pos, pos+nt_chunk, ring_buffer::ACCESS_APPEND);

    bool ret = _fill_chunk(intensity, rb_intensity->get_stride(), weights, rb_intensity->get_stride(), pos);    

    rb_intensity->put(intensity, pos, pos+nt_chunk, ring_buffer::ACCESS_APPEND);
    rb_weights->put(weights, pos, pos+nt_chunk, ring_buffer::ACCESS_APPEND);
    
    return ret;
}


// Default virtual
void wi_stream::_bind_stream(Json::Value &json_attrs)
{
    return;
}


}  // namespace rf_pipelines
