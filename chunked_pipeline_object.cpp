#include "rf_pipelines_internals.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
} // emacs pacifier
#endif


chunked_pipeline_object::chunked_pipeline_object(const string &name_, bool can_be_first_, ssize_t nt_chunk_, ssize_t nt_chunk_min_) :
    pipeline_object(name_),
    can_be_first(can_be_first_),
    nt_chunk(nt_chunk_),
    nt_chunk_min(nt_chunk_min_)
{
    if (nt_chunk < 0)
	_throw("expected nt_chunk >= 0 in chunked_pipeline_object constructor");
    if (nt_chunk_min < 0)
	_throw("expected nt_chunk_min >= 0 in chunked_pipeline_object constructor");
}


// virtual override
ssize_t chunked_pipeline_object::get_preferred_chunk_size()
{
    if (!can_be_first)
	return 0;
    if (nt_chunk == 0)
	_throw("in chunked_pipeline_objects with can_be_first=true, nt_chunk must be initialized to a nonzero value before bind() is called");
    return nt_chunk;
}


// Helper function, no-ops if nt_chunk has already been initialized to a nonzero value.
// Can be called any time during initialization or bind(), but at latest will be called at the end of bind().
void chunked_pipeline_object::finalize_nt_chunk()
{
    if (nt_chunk > 0) {
	this->_check_nt_chunk();
	return;
    }

    ssize_t m = max(nt_chunk_in, ssize_t(512));
    ssize_t n = (nt_chunk_min > 0) ? nt_chunk_min : 1;

    for (const auto &p: this->all_ring_buffers)
	n = lcm(n, p->nds);

    this->nt_chunk = n * max(m/n, ssize_t(1));
    this->_check_nt_chunk();
}


// Internal helper function, assumes nt_chunk has been initialized.
void chunked_pipeline_object::_check_nt_chunk() const
{
    rf_assert(nt_chunk > 0);

    if ((nt_chunk_min > 0) && (nt_chunk % nt_chunk_min))
	_throw("nt_chunk (=" + to_string(nt_chunk) + ") must be a multiple of nt_chunk_min (=" + to_string(nt_chunk_min) + ")");

    for (const auto &p: this->all_ring_buffers) {
	if (nt_chunk % p->nds)
	    _throw("nt_chunk (=" + to_string(nt_chunk) + ") must be a multiple of all ring buffer downsampling factors (found nds=" + to_string(p->nds) + ")");
    }
}


// virtual override
void chunked_pipeline_object::_bind(ring_buffer_dict &rb_dict, Json::Value &json_attrs)
{
    this->_prebind_nt_chunk = nt_chunk;
    
    // Note: all calls to get_buffer() or create_buffer() are in _bindc(), which is defined by the subclass.
    this->_bindc(rb_dict, json_attrs);
    
    this->finalize_nt_chunk();
    this->nt_chunk_out = (nt_chunk_in % nt_chunk) ? nt_chunk : nt_chunk_in;
    this->nt_maxgap = nt_chunk - gcd(nt_chunk_in, nt_chunk);
    this->nt_contig = nt_chunk;
}


void chunked_pipeline_object::_unbind()
{
    // We revert 'nt_chunk' to its "prebind" value.
    this->nt_chunk = this->get_prebind_nt_chunk();
    this->_unbindc();
}


// virtual override
ssize_t chunked_pipeline_object::_advance()
{
    ssize_t ret = SSIZE_MAX;

    while (pos_lo <= pos_hi - nt_chunk) {
	bool alive = _process_chunk(pos_lo);
	if (!alive)
	    ret = min(ret, pos_hi);
	
	pos_lo += nt_chunk;
    }

    return ret;
}


// default virtual
void chunked_pipeline_object::_unbindc() { }


}  // namespace rf_pipelines
