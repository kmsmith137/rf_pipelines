#include "rf_pipelines_internals.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}  // emacs pacifier
#endif


zoomable_tileset::zoomable_tileset(const vector<vector<ssize_t>> &cdims_, ssize_t ny_arr_, ssize_t nds_arr_) :
    cdims(cdims_), 
    ny_arr(ny_arr_),
    nds_arr(nds_arr_)
{
    if (cdims.size() == 0)
	throw runtime_error("rf_pipelines::zoomable_tileset constructor: expected number of ring buffers to be > 0");
    
    for (size_t i = 0; i < cdims.size(); i++)
	ring_buffer::check_cdims(cdims[i]);

    if (ny_arr <= 0)
	throw runtime_error("rf_pipelines::zoomable_tileset constructor: expected ny_arr > 0");
    if (nds_arr <= 0)
	throw runtime_error("rf_pipelines::zoomable_tileset constructor: expected nds_arr > 0");
    if (!is_power_of_two(nds_arr))
	throw runtime_error("rf_pipelines::zoomable_tileset constructor: expected nds_arr to be a power of two");
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


// -------------------------------------------------------------------------------------------------


zoomable_tileset_state::zoomable_tileset_state(const shared_ptr<zoomable_tileset> &zt_, const shared_ptr<outdir_manager> &mp_, const Json::Value &json_attrs, ssize_t img_ny_) :
    zt(zt_),
    mp(mp_),
    img_nzoom(ssize_t_from_json(json_attrs, "img_nzoom")),
    img_nds(ssize_t_from_json(json_attrs, "img_nds")),
    img_nx(ssize_t_from_json(json_attrs, "img_nx")),
    img_ny(img_ny_)
{
    if (!zt)
	throw runtime_error("zoomable_tileset constructor: expected nonempty shared_ptr<zoomable_tileset>");
    if (!mp)
	throw runtime_error("zoomable_tileset constructor: expected nonempty shared_ptr<outdir_manager>");
    if (img_nzoom <= 0)
	throw runtime_error("zoomable_tileset constructor: expected img_nzoom > 0");
    if (img_nds <= 0)
	throw runtime_error("zoomable_tileset constructor: expected img_nds > 0");
    if (!is_power_of_two(img_nds))
	throw runtime_error("zoomable_tileset constructor: expected img_nds to be a power of two");
    if (img_nx <= 0)
	throw runtime_error("zoomable_tileset constructor: expected img_nx > 0");
    if (img_ny <= 0)
	throw runtime_error("zoomable_tileset constructor: expected img_ny > 0");
    if (xdiv(img_ny, zt->ny_arr) != 0)
	throw runtime_error("zoomable_tileset::ny_arr must be equal to, or a divisor of, img_ny");
    
    this->ds_offset = integer_log2(img_nds) - integer_log2(zt->nds_arr);

    int nouter = img_nzoom + max(ds_offset,0);
    int ninner = zt->cdims.size();

    this->ring_buffers.resize(nouter);

    for (int i = 0; i < nouter; i++) {
	// FIXME I think this slightly overallocates the ring buffers.
	int nds = zt->nds_arr * (1 << i);
	int nt_contig = nds * (2 * img_nx);

	this->ring_buffers[i].resize(ninner);

	for (int j = 0; j < ninner; j++) {
	    this->ring_buffers[i][j] = make_shared<ring_buffer> (zt->cdims[j], nds);
	    this->ring_buffers[i][j]->update_params(nt_contig, nt_contig);   // (nt_contig, nt_maxlag)
	}
    }
}


// The arguments (nt_contig, nt_maxlag) have the same meaning as in ring_buffer::update_params().
void zoomable_tileset_state::update_params(ssize_t nt_contig, ssize_t nt_maxlag)
{
    int ninner = zt->cdims.size();
    int nt0 = zt->nds_arr * (2 * img_nx);

    rf_assert(ring_buffers.size() > 0);
    rf_assert(ring_buffers[0].size() == ninner);
    
    for (int i = 0; i < ninner; i++)
	this->ring_buffers[0][i]->update_params(nt0 + nt_contig, nt0 + nt_maxlag);
} 


void zoomable_tileset_state::allocate()
{
    for (size_t i = 0; i < ring_buffers.size(); i++)
	for (size_t j = 0; j < ring_buffers[i].size(); j++)
	    ring_buffers[i][j]->allocate();
}


void zoomable_tileset_state::deallocate()
{
    for (size_t i = 0; i < ring_buffers.size(); i++)
	for (size_t j = 0; j < ring_buffers[i].size(); j++)
	    ring_buffers[i][j]->deallocate();
}


}  // namespace rf_pipelines
