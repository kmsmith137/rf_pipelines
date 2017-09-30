#include "rf_pipelines_internals.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}  // namespace rf_pipelines
#endif


wi_transform::wi_transform(const string &name_, ssize_t nt_chunk_, ssize_t nfreq_, ssize_t nds_) :
    chunked_pipeline_object(name_, false, nt_chunk_),  // can_be_first=false
    nfreq(nfreq_),
    nds(nds_)
{
    if (nfreq < 0)
	_throw("expected nfreq >= 0 in wi_transform constructor");
    if (nds < 0)
	_throw("expected nds >= 0 in wi_transform constructor");
}


// virtual override
void wi_transform::_bind_chunked(ring_buffer_dict &rb_dict, Json::Value &json_data)
{
    this->_save_nfreq = nfreq;
    this->_save_nds = nds;

    this->rb_intensity = this->get_buffer(rb_dict, "INTENSITY");
    this->rb_weights = this->get_buffer(rb_dict, "WEIGHTS");

    if (rb_intensity->cdims != rb_weights->cdims)
	_throw("'intensity' and 'weights' buffers have different dimensions");
    if (rb_intensity->nds != rb_weights->nds)
	_throw("'intensity' and 'weights' buffers have different downsampling");
    if (rb_intensity->cdims.size() != 1)
	_throw ("expected intensity/weights arrays to be two-dimensional");

    ssize_t expected_nfreq = this->nfreq;
    ssize_t expected_nds = this->nds;

    this->nfreq = rb_intensity->cdims[0];
    this->nds = rb_intensity->nds;

    if ((expected_nfreq > 0) && (nfreq != expected_nfreq))
	_throw("pipeline nfreq (" + to_string(nfreq) + ") doesn't match expected value (" + to_string(expected_nfreq) + ")");
    if ((expected_nds > 0) && (nds != expected_nds))
	_throw("pipeline nds (" + to_string(nds) + ") doesn't match expected value (" + to_string(expected_nds) + ")");

    // Optional subclass-specific initializations.
    this->_bind_transform(json_data);
}


// virtual override
bool wi_transform::_process_chunk(ssize_t pos)
{
    float *intensity = rb_intensity->get(pos, pos+nt_chunk, ring_buffer::ACCESS_RW);
    float *weights = rb_weights->get(pos, pos+nt_chunk, ring_buffer::ACCESS_RW);

    this->_process_chunk(intensity, rb_intensity->get_stride(), weights, rb_weights->get_stride(), pos);

    rb_intensity->put(intensity, pos, pos+nt_chunk, ring_buffer::ACCESS_RW);
    rb_weights->put(weights, pos, pos+nt_chunk, ring_buffer::ACCESS_RW);
    return true;
}


// Default virtual
void wi_transform::_bind_transform(Json::Value &json_data)
{
    return;
}


ssize_t wi_transform::get_orig_nfreq() const
{
    return is_bound() ? _save_nfreq : nfreq;
}


ssize_t wi_transform::get_orig_nds() const
{
    return is_bound() ? _save_nds : nds;
}


}  // namespace rf_pipelines

