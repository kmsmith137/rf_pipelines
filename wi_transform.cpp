#include "rf_pipelines_internals.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}  // namespace rf_pipelines
#endif


wi_transform::wi_transform(const string &class_name_, const string &name_) :
    chunked_pipeline_object(class_name_, name_, false)  // can_be_first=false
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
    if ((expected_nds == 1) && (nds != 1))
	_throw("this transform does not (yet?) support downsampling, and cannot be used in a wi_sub_pipeline");
    if ((expected_nds > 0) && (nds != expected_nds))
	_throw("pipeline nds (" + to_string(nds) + ") doesn't match expected value (" + to_string(expected_nds) + ")");

    // Optional subclass-specific initializations.
    this->_bind_transform(json_attrs);
    this->_bind_transform_rb(rb_dict);
}


// virtual override
bool wi_transform::_process_chunk(ssize_t pos)
{
    ring_buffer_subarray intensity(rb_intensity, pos, pos+nt_chunk, ring_buffer::ACCESS_RW);
    ring_buffer_subarray weights(rb_weights, pos, pos+nt_chunk, ring_buffer::ACCESS_RW);

    this->_process_chunk(intensity.data, intensity.stride, weights.data, weights.stride, pos);
    return true;
}


void wi_transform::_unbindc()
{
    this->_unbind_transform();
    
    this->rb_intensity.reset();
    this->rb_weights.reset();

    // We revert 'nfreq' and 'nds' to their "prebind" values.
    this->nfreq = this->get_prebind_nfreq();
    this->nds = this->get_prebind_nds();
}


// We override chunked_pipeline::finalize_nt_chunk() in order to incorporate 'kernel_chunk_size'.
void wi_transform::finalize_nt_chunk()
{
    if (nt_chunk_in <= 0)
	_throw("finalize_nt_chunk(): expected nt_chunk_in > 0.  Note that finalize_nt_chunk() should be called during bind(), after ring buffers are allocated");

    // I don't think there is any way to fail this check, while passing the previous one, but you never know...
    if (nds <= 0)
	_throw("wi_transform::finalize_nt_chunk(): internal error: expected nds > 0");

    if (nt_chunk > 0) {
	this->_check_nt_chunk();
	return;
    }

    ssize_t m = max(nt_chunk_in, ssize_t(512));
    ssize_t n = max(kernel_chunk_size, ssize_t(1)) * nds;

    this->nt_chunk = n * max(m/n, ssize_t(1));
    this->_check_nt_chunk();
}


// We override chunked_pipeline::_check_nt_chunk() in order to incorporate 'kernel_chunk_size'.
void wi_transform::_check_nt_chunk() const
{
    rf_assert(nds > 0);
    rf_assert(nt_chunk > 0);
    rf_assert(nt_chunk_in > 0);

    if (nt_chunk % nds)
	_throw("nt_chunk (=" + to_string(nt_chunk) + ") must be a multiple of nds (=" + to_string(nds) + ")");

    if ((kernel_chunk_size > 0) && (nt_chunk % (kernel_chunk_size * nds)))
	_throw("nt_chunk (=" + to_string(nt_chunk) + ") must be a multiple of kernel_chunk_size*nds (=" + to_string(kernel_chunk_size*nds) + ")");
}


// Default virtuals
void wi_transform::_bind_transform(Json::Value &json_attrs) { }
void wi_transform::_bind_transform_rb(ring_buffer_dict &rb_dict) { }
void wi_transform::_unbind_transform() { }


}  // namespace rf_pipelines
