#include "rf_pipelines_internals.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}  // emacs pacifier
#endif


zoomable_tileset::zoomable_tileset(const vector<vector<ssize_t>> &cdims_, ssize_t ny_arr_, ssize_t nds_min_, ssize_t nds_max_) :
    cdims(cdims_), 
    ny_arr(ny_arr_),
    nds_min(nds_min_),
    nds_max(nds_max_)
{
    if (cdims.size() == 0)
	throw runtime_error("rf_pipelines::zoomable_tileset constructor: expected number of ring buffers to be > 0");
    
    for (size_t i = 0; i < cdims.size(); i++)
	ring_buffer::check_cdims(cdims[i]);

    if (ny_arr <= 0)
	throw runtime_error("rf_pipelines::zoomable_tileset constructor: expected ny_arr > 0");
    if (nds_min < 0)
	throw runtime_error("rf_pipelines::zoomable_tileset constructor: expected nds_min >= 0");
    if (nds_max < 0)
	throw runtime_error("rf_pipelines::zoomable_tileset constructor: expected nds_max >= 0");
    if ((nds_min > 0) && !is_power_of_two(nds_min))
	throw runtime_error("rf_pipelines::zoomable_tileset constructor: expected nds_min to be a power of two");
    if ((nds_max > 0) && !is_power_of_two(nds_max))
	throw runtime_error("rf_pipelines::zoomable_tileset constructor: expected nds_max to be a power of two");
    if ((nds_min > 0) && (nds_max > 0) && (nds_min > nds_max))
	throw runtime_error("rf_pipelines::zoomable_tileset constructor: expected nds_min < nds_max");
}


// Default virtual
void zoomable_tileset::downsample_rbvec(rbvec_t &rb_out, rbvec_t &rb_in, ssize_t pos, ssize_t nt)
{
    for (size_t i = 0; i < rb_out.size(); i++) {
	ssize_t nt_out = xdiv(nt, rb_out[i]->nds);
	ssize_t csize = rb_out[i]->csize;
	
	ring_buffer_subarray a_out(rb_out[i], pos, pos+nt, ring_buffer::ACCESS_APPEND);
	ring_buffer_subarray a_in(rb_in[i], pos, pos+nt, ring_buffer::ACCESS_READ);

	ssize_t ostride = a_out.stride;
	ssize_t istride = a_in.stride;

	// FIXME (low-priority): assembly-language kernel would speed this up.
	// But plotting is unlikely to be a bottleneck!

	for (ssize_t j = 0; j < csize; j++) {
	    float *p_out = a_out.data + j*ostride;
	    float *p_in = a_in.data + j*istride;

	    for (ssize_t k = 0; k < nt_out; k++)
		p_out[k] = 0.5 * (p_in[2*k] + p_in[2*k+1]);
	}
    }
}

// Default virtual
void zoomable_tileset::extend_rbvec(rbvec_t &rb_out, ssize_t pos, ssize_t nt)
{
    for (size_t i = 0; i < rb_out.size(); i++) {
	ssize_t nt_out = xdiv(nt, rb_out[i]->nds);
	ssize_t csize = rb_out[i]->csize;

	ring_buffer_subarray a(rb_out[i], pos, pos+nt, ring_buffer::ACCESS_APPEND);
	ssize_t stride = a.stride;

	for (ssize_t j = 0; j < csize; j++)
	    memset(a.data + j*stride, 0, nt_out * sizeof(float));
    }
}


}  // namespace rf_pipelines
