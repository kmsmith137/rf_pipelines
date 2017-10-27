#include "rf_pipelines_internals.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}  // namespace rf_pipelines
#endif


wi_stream::wi_stream(const string &class_name_, const string &name_) :
    chunked_pipeline_object(class_name_, name_, true)
{ }


void wi_stream::_bindc(ring_buffer_dict &rb_dict, Json::Value &json_attrs)
{
    // Optional subclass-specific initializations.
    this->_bind_stream(json_attrs);

    if (nfreq <= 0)
	_throw("'nfreq' must be initialized in either constructor, or in _bind_stream()");
    if (nt_chunk <= 0)
	_throw("'nt_chunk' must be initialized in either constructor, or in _bind_stream()");

    this->rb_intensity = this->create_buffer(rb_dict, "INTENSITY", {nfreq}, 1);
    this->rb_weights = this->create_buffer(rb_dict, "WEIGHTS", {nfreq}, 1);

    json_attrs["nfreq"] = Json::Int64(nfreq);
}


bool wi_stream::_process_chunk(ssize_t pos)
{
    ring_buffer_subarray intensity(rb_intensity, pos, pos+nt_chunk, ring_buffer::ACCESS_APPEND);
    ring_buffer_subarray weights(rb_weights, pos, pos+nt_chunk, ring_buffer::ACCESS_APPEND);

    return _fill_chunk(intensity.data, intensity.stride, weights.data, weights.stride, pos);
}


void wi_stream::_unbindc()
{
    this->_unbind_stream();
    this->rb_intensity.reset();
    this->rb_weights.reset();
}


// Default virtuals
void wi_stream::_bind_stream(Json::Value &json_attrs) { }
void wi_stream::_unbind_stream() { }


}  // namespace rf_pipelines
