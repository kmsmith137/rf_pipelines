#include "rf_pipelines_internals.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
} // emacs pacifier
#endif


chunked_pipeline_object::chunked_pipeline_object(const string &name_, bool can_be_first_, ssize_t nt_chunk_) :
    pipeline_object(name_),
    can_be_first(can_be_first_),
    nt_chunk(nt_chunk_)
{
    if (nt_chunk < 0)
	_throw("expected nt_chunk >= 0 in chunked_pipeline_object constructor");
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


// virtual override
void chunked_pipeline_object::_bind(ring_buffer_dict &rb_dict, Json::Value &json_attrs)
{
    this->_prebind_nt_chunk = nt_chunk;

    if (nt_chunk == 0)
	nt_chunk = nt_chunk_in;

    this->nt_chunk_out = (nt_chunk_in % nt_chunk) ? nt_chunk : nt_chunk_in;
    this->nt_maxgap = nt_chunk - gcd(nt_chunk_in, nt_chunk);
    this->nt_contig = nt_chunk;
    
    // Note: all calls to get_buffer() or create_buffer() are in _bindc(), which is defined by the subclass.
    this->_bindc(rb_dict, json_attrs);
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
