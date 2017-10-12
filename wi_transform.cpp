#include "rf_pipelines_internals.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}  // namespace rf_pipelines
#endif


wi_transform::wi_transform(const string &name_) :
    chunked_pipeline_object(name_, false)  // can_be_first=false
{ }


// virtual override
void wi_transform::_bindc(ring_buffer_dict &rb_dict, Json::Value &json_attrs)
{
    this->_prebind_nfreq = nfreq;
    this->_prebind_nds = nds;

    this->rb_intensity = this->get_buffer(rb_dict, "INTENSITY");
    this->rb_weights = this->get_buffer(rb_dict, "WEIGHTS");

    if (rb_intensity->cdims != rb_weights->cdims)
	_throw("'intensity' and 'weights' buffers have different dimensions");
    if (rb_intensity->nds != rb_weights->nds)
	_throw("'intensity' and 'weights' buffers have different downsampling");
    if (rb_intensity->cdims.size() != 1)
	_throw("expected intensity/weights arrays to be two-dimensional");

    ssize_t expected_nfreq = this->nfreq;
    ssize_t expected_nds = this->nds;

    this->nfreq = rb_intensity->cdims[0];
    this->nds = rb_intensity->nds;

    if ((expected_nfreq > 0) && (nfreq != expected_nfreq))
	_throw("pipeline nfreq (" + to_string(nfreq) + ") doesn't match expected value (" + to_string(expected_nfreq) + ")");
    if ((expected_nds > 0) && (nds != expected_nds))
	_throw("pipeline nds (" + to_string(nds) + ") doesn't match expected value (" + to_string(expected_nds) + ")");

    // Optional subclass-specific initializations.
    this->_bind_transform(json_attrs);
}


// virtual override
bool wi_transform::_process_chunk(ssize_t pos)
{
    ring_buffer_subarray intensity(rb_intensity, pos, pos+nt_chunk, ring_buffer::ACCESS_RW);
    ring_buffer_subarray weights(rb_weights, pos, pos+nt_chunk, ring_buffer::ACCESS_RW);

    this->_process_chunk(intensity.data, intensity.stride, weights.data, weights.stride, pos);
    return true;
}


// Default virtual
void wi_transform::_bind_transform(Json::Value &json_attrs)
{
    return;
}


}  // namespace rf_pipelines

