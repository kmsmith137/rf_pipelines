#include "rf_pipelines_internals.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif


// FIXME put somewhere more general?
template<typename T>
inline void strided_2d_copy(T *dst, const T *src, ssize_t m, ssize_t n, ssize_t dst_stride, ssize_t src_stride)
{
    for (ssize_t i = 0; i < m; i++)
	memcpy(dst + i*dst_stride, src + i*src_stride, n * sizeof(T));
}


pipeline_spool::pipeline_spool(const vector<string> &bufnames_) :
    pipeline_object("pipeline_spool")
{
    this->bufnames = bufnames_;
}

// Overrides pipeline_object::_bind()
void pipeline_spool::_bind(ring_buffer_dict &rb_dict, Json::Value &json_attrs)
{
    this->nt_chunk_out = nt_chunk_in;
    this->nt_contig = nt_chunk_in;
    this->nt_maxgap = 0;

    this->ring_buffers.resize(bufnames.size());
    this->spooled_buffers.resize(bufnames.size());

    for (unsigned int i = 0; i < bufnames.size(); i++) {
	shared_ptr<ring_buffer> rb = this->get_buffer(rb_dict, bufnames[i]);

	rf_assert(rb->nds > 0);
	rf_assert(nt_chunk_in % rb->nds == 0);

	shared_ptr<spooled_buffer> sb = make_shared<spooled_buffer> ();
	sb->cdims = rb->cdims;
	sb->csize = rb->csize;
	sb->nds = rb->nds;
	sb->nt = 0;

	ring_buffers[i] = rb;
	spooled_buffers[i] = sb;
    }
}


// Overrides pipeline_object::_advance()
ssize_t pipeline_spool::_advance()
{
    rf_assert(pos_lo % nt_chunk_in == 0);
    rf_assert(pos_hi % nt_chunk_in == 0);
    rf_assert(ring_buffers.size() == bufnames.size());
    rf_assert(spooled_buffers.size() == bufnames.size());

    while (pos_lo < pos_hi) {
	for (unsigned int i = 0; i < bufnames.size(); i++) {
	    shared_ptr<ring_buffer> rb = this->ring_buffers[i];
	    shared_ptr<spooled_buffer> sb = this->spooled_buffers[i];
	    
	    ssize_t csize = sb->csize;
	    ssize_t nt_ds = xdiv(nt_chunk_in, sb->nds);

	    shared_ptr<float> dst = make_sptr<float> (csize * nt_ds);
	    ring_buffer_subarray src(rb, pos_lo, pos_lo + nt_chunk_in, ring_buffer::ACCESS_READ);

	    rf_assert(sb->nt == pos_lo);
	    sb->_incremental_data.push_back(dst);
	    sb->nt += nt_chunk_in;

	    strided_2d_copy(dst.get(), src.data, csize, nt_ds, nt_ds, src.stride);
	}

	pos_lo += nt_chunk_in;
    }

    return SSIZE_MAX;
}


// Overrides pipeline_object::_end_pipeline()
void pipeline_spool::_end_pipeline(Json::Value &j)
{
    rf_assert(ring_buffers.size() == bufnames.size());
    rf_assert(spooled_buffers.size() == bufnames.size());

    for (unsigned int i = 0; i < bufnames.size(); i++) {
	shared_ptr<spooled_buffer> sb = this->spooled_buffers[i];
	
	rf_assert(sb->_incremental_data.size() == xdiv(sb->nt, nt_chunk_in));
	rf_assert(sb->data.size() == 0);

	ssize_t csize = sb->csize;
	ssize_t nchunks = sb->_incremental_data.size();
	ssize_t nt_ds = xdiv(nt_chunk_in, sb->nds);

	sb->data.resize(csize * nchunks * nt_ds);

	for (ssize_t ichunk = 0; ichunk < nchunks; ichunk++) {
	    strided_2d_copy(&sb->data[ichunk*nt_ds],               // dst pointer
			    sb->_incremental_data[ichunk].get(),   // src pointer
			    csize, nt_ds,                          // 2-d array dimensions
			    nchunks * nt_ds,                       // dst stride
			    nt_ds);                                // src stride
	}

	sb->_incremental_data.clear();
    }
}


// Overrides pipeline_object::_unbind()
void pipeline_spool::_unbind()
{
    this->ring_buffers.clear();
}


// Overrides pipeline_object::_deallocate()
void pipeline_spool::_deallocate()
{
    this->ring_buffers.clear();
    this->spooled_buffers.clear();
}


shared_ptr<pipeline_spool::spooled_buffer> pipeline_spool::get_spooled_buffer(const string &bufname) const
{
    if (state != DONE)
	throw runtime_error("rf_pipelines::pipeline_spool::get_spooled_buffer() must be called after pipeline exit");

    rf_assert(spooled_buffers.size() == bufnames.size());

    for (unsigned int i = 0; i < bufnames.size(); i++) {
	if (bufnames[i] == bufname) {
	    rf_assert(spooled_buffers[i]->_incremental_data.size() == 0);
	    return spooled_buffers[i];
	}
    }

    throw runtime_error("rf_pipelines::pipeline_spool::get_spooled_buffer(): bufname '" + bufname + "' not found");
}


}  // namespace rf_pipelines
