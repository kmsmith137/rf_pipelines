#include "rf_pipelines_internals.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}  // emacs pacifier
#endif


zoomable_tileset::zoomable_tileset(const vector<vector<ssize_t>> &cdims_, const string &img_prefix_, ssize_t img_ny_, ssize_t ny_arr_, ssize_t nds_arr_) :
    cdims(cdims_), 
    img_prefix(img_prefix_),
    img_ny(img_ny_),
    ny_arr(ny_arr_),
    nds_arr(nds_arr_)
{
    if (cdims.size() == 0)
	throw runtime_error("rf_pipelines::zoomable_tileset constructor: expected number of ring buffers to be > 0");
    
    for (size_t i = 0; i < cdims.size(); i++)
	ring_buffer::check_cdims(cdims[i]);

    if (img_prefix.size() == 0)
	throw runtime_error("rf_pipelines::zoomable_tileset constructor: expected 'img_prefix' to be a nonempty string");
    if (img_ny <= 0)
	throw runtime_error("rf_pipelines::zoomable_tileset constructor: expected img_ny > 0");
    if (ny_arr <= 0)
	throw runtime_error("rf_pipelines::zoomable_tileset constructor: expected ny_arr > 0");
    if (nds_arr <= 0)
	throw runtime_error("rf_pipelines::zoomable_tileset constructor: expected nds_arr > 0");
    if (!is_power_of_two(nds_arr))
	throw runtime_error("rf_pipelines::zoomable_tileset constructor: expected nds_arr to be a power of two");
    if ((img_ny % ny_arr) != 0)
	throw runtime_error("rf_pipelines::zoomable_tileset constructor: ny_arr must be equal to, or a divisor of, img_ny");
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


// This helper ensures that the zoomable_tileset_state constructor doesn't segfault on an empty shared_ptr<zoomable_tileset>.
inline shared_ptr<zoomable_tileset> _check_zt(const shared_ptr<zoomable_tileset> &zt)
{
    if (!zt)
	throw runtime_error("zoomable_tileset constructor: expected nonempty shared_ptr<zoomable_tileset>");
    return zt;
}


zoomable_tileset_state::zoomable_tileset_state(const shared_ptr<zoomable_tileset> &zt_, const pipeline_object &p) :
    zt(_check_zt(zt_)),
    mp(p.out_mp),
    img_prefix(zt_->img_prefix),
    img_nzoom(p._params.img_nzoom),
    img_nds(p._params.img_nds),
    img_nx(p._params.img_nx),
    img_ny(zt_->img_ny),
    nds_arr(zt_->nds_arr),
    ny_arr(zt_->ny_arr),
    debug(p._params.debug)
{
    // These asserts should have been checked previously, either in run_params::check()
    // or in the zoomable_tileset constructor.

    rf_assert(mp);
    rf_assert(mp->outdir.size() > 0);
    rf_assert(img_prefix.size() > 0);
    rf_assert((img_nzoom >= 1) && (img_nzoom <= 10));
    rf_assert((img_nds > 0) && is_power_of_two(img_nds));
    rf_assert((nds_arr > 0) && is_power_of_two(nds_arr));
    rf_assert((img_nx > 0) && (img_nx % 2 == 0));
    rf_assert((img_ny > 0) && (ny_arr > 0) && (img_ny % ny_arr == 0));
    rf_assert(p.state == pipeline_object::BINDING);
    
    this->ds_offset = integer_log2(img_nds) - integer_log2(nds_arr);

    int nouter = max(img_nzoom + ds_offset, ssize_t(1));
    int ninner = zt->cdims.size();

    this->ring_buffers.resize(nouter);

    for (int i = 0; i < nouter; i++) {
	this->ring_buffers[i].resize(ninner);

	// FIXME I think this slightly overallocates the ring buffers.
	int nds = nds_arr * (1 << i);
	int nt_contig = 2 * this->nt_per_block(i);

	for (int j = 0; j < ninner; j++) {
	    this->ring_buffers[i][j] = make_shared<ring_buffer> (zt->cdims[j], nds, debug);
	    this->ring_buffers[i][j]->update_params(nt_contig, nt_contig);   // (nt_contig, nt_maxlag)
	}
    }

    this->nblocks.resize(nouter, ssize_t(0));
    this->rgb_zoom.resize(img_nzoom, nullptr);
    this->_initialize_json();
}


// The arguments (nt_contig, nt_maxlag) have the same meaning as in ring_buffer::update_params().
void zoomable_tileset_state::update_params(ssize_t nt_contig, ssize_t nt_maxlag)
{
    rf_assert(!is_allocated);
    rf_assert(ring_buffers.size() > 0);  // paranoid

    int nt0 = 2 * nt_per_block(0);

    for (size_t i = 0; i < ring_buffers[0].size(); i++)
	ring_buffers[0][i]->update_params(nt0 + nt_contig, nt0 + nt_maxlag);
} 


void zoomable_tileset_state::allocate()
{
    if (is_allocated)
	return;

    for (size_t i = 0; i < ring_buffers.size(); i++)
	for (size_t j = 0; j < ring_buffers[i].size(); j++)
	    ring_buffers[i][j]->allocate();

    // FIXME a little overallocation here.
    for (ssize_t i = 0; i < img_nzoom; i++) {
	uint8_t *p = new uint8_t[ny_arr * img_nx * 3];
	this->rgb_alloc.push_back(unique_ptr<uint8_t[]> (p));
	this->rgb_zoom[i] = p;
    }

    if (ny_arr != img_ny) {
	uint8_t *p = new uint8_t[img_ny * img_nx * 3];
	this->rgb_alloc.push_back(unique_ptr<uint8_t[]> (p));
	this->rgb_us = p;
    }

    this->is_allocated = true;
}


void zoomable_tileset_state::deallocate()
{
    if (!is_allocated)
	return;

    this->reset();

    for (size_t i = 0; i < ring_buffers.size(); i++)
	for (size_t j = 0; j < ring_buffers[i].size(); j++)
	    ring_buffers[i][j]->deallocate();

    for (ssize_t i = 0; i < img_nzoom; i++)
	this->rgb_zoom[i] = nullptr;

    this->rgb_us = nullptr;
    this->rgb_alloc.clear();
    this->is_allocated = false;
}


void zoomable_tileset_state::reset()
{
    if (!is_allocated)
	throw std::runtime_error("rf_pipelines internal error: zoomable_tileset_state::reset() called with is_allocated=false");

    for (size_t i = 0; i < ring_buffers.size(); i++)
	for (size_t j = 0; j < ring_buffers[i].size(); j++)
	    ring_buffers[i][j]->reset();

    for (size_t i = 0; i < ring_buffers.size(); i++)
	nblocks[i] = 0;

    this->curr_pos = 0;
    this->_initialize_json();
    this->is_flushed = false;
}


void zoomable_tileset_state::advance(ssize_t pos)
{
    if (!is_allocated)
	throw std::runtime_error("rf_pipelines internal error: zoomable_tileset_state::advance() called with is_allocated=false");
    if (is_flushed)
	throw std::runtime_error("rf_pipelines internal error: zoomable_tileset_state::advance() called with is_flushed=true");
    
    rf_assert(pos >= curr_pos);

    ssize_t nb_f = pos / this->nt_per_block(0);

    while (nblocks[0] < nb_f) {
	_advance_by_one_block(0);

	for (size_t irb = 1; irb < ring_buffers.size(); irb++) {
	    ssize_t nb0 = nblocks[irb];
	    ssize_t nb1 = nblocks[irb-1];
	    ssize_t nt0 = nt_per_block(irb);

	    if (nb1 == 2*nb0+1)
		break;

	    rf_assert(nb1 == 2*nb0+2);

	    zt->downsample_rbvec(ring_buffers[irb], ring_buffers[irb-1], nb0*nt0, nt0);
	    this->_advance_by_one_block(irb);
	}
    }

    this->curr_pos = pos;
}


void zoomable_tileset_state::flush()
{
    if (!is_allocated)
	throw std::runtime_error("rf_pipelines internal error: zoomable_tileset_state::flush() called with is_allocated=false");
    if (is_flushed)
	return;

    ssize_t pos = this->curr_pos;

    for (size_t irb = 0; irb < ring_buffers.size(); irb++) {
	ssize_t pos0 = (nblocks[irb]) * this->nt_per_block(irb);
	ssize_t pos1 = (nblocks[irb]+1) * this->nt_per_block(irb);

	rf_assert((pos >= pos0) && (pos <= pos1));
	
	if ((pos == pos0) || (pos == pos1))
	    continue;

	if (irb > 0)
	    zt->downsample_rbvec(ring_buffers[irb], ring_buffers[irb-1], pos0, pos-pos0);

	zt->extend_rbvec(ring_buffers[irb], pos, pos1-pos);

	_advance_by_one_block(irb);
	pos = pos1;
    }

    this->is_flushed = true;
}


// _advance_by_one_block(irb)
//
// Called when there is enough data in ring buffer 'irb' to complete a block.
// One "block" consists of (nds_arr * img_nx * (1 << rb)) time samples.

void zoomable_tileset_state::_advance_by_one_block(int irb)
{
    rf_assert((irb >= 0) && (irb < int(ring_buffers.size())));

    if (irb >= ds_offset) {
	int izoom = irb - ds_offset;
	rf_assert(izoom >= 0 && izoom < img_nzoom);

	ssize_t nb = nblocks[irb];
	ssize_t nt = nt_per_block(irb);

	zt->plot_rbvec(rgb_zoom[izoom], ring_buffers[irb], nb*nt, nt);
	_emit_plot(izoom, nb);
    }

    this->nblocks[irb]++;
}


// _emit_plot(izoom, iplot)
// Called when the rgb_zoom[izoom] array is filled.

void zoomable_tileset_state::_emit_plot(int izoom, ssize_t iplot)
{
    uint8_t *rgb = rgb_zoom[izoom];

    // We may need to upsample from shape (ny_arr, img_nx, 3) to shape (img_ny, img_nx, 3).
    if (img_ny != ny_arr) {
	ssize_t Dy = xdiv(img_ny, ny_arr);
	ssize_t stride = 3 * img_nx;

	for (ssize_t i = 0; i < ny_arr; i++)
	    for (ssize_t j = 0; j < Dy; j++)
		memcpy(rgb_us + (i*Dy+j)*stride, rgb + i*stride, stride);

	rgb = rgb_us;
    }

    // Now 'rgb' is a shape (img_ny, img_nx, 3) array.
    // Next step is to write the png file.
    
    stringstream ss;
    ss << "img_" << izoom << "_" << iplot << ".png";

    string basename = ss.str();
    string fullname = mp->add_file(basename);

    // write_rgb8_png(filename, rgb, m, n, ymajor, ytop_to_bottom)
    write_rgb8_png(fullname, rgb, img_ny, img_nx, true, false);

    Json::Value j;
    j["filename"] = fullname;
    j["it0"] = Json::Int64(iplot * img_nds * ssize_t(1 << izoom));
    j["nx"] = Json::Int64(img_nx);

    json_output[izoom]["files"].append(j);
    
    // This part applies in the case where the lowest ring buffer is downsampled 
    // relative to the tileset.  Whenever a block is ready in the lowest ring buffer,
    // we recurisvely emit a "tree" of upsampled tiles.

    if ((izoom > 0) && (izoom <= -ds_offset)) {
	_upsample_rgb(rgb_zoom[izoom-1], rgb_zoom[izoom], false);
	_emit_plot(izoom-1, 2*iplot);
	
	_upsample_rgb(rgb_zoom[izoom-1], rgb_zoom[izoom], true);
	_emit_plot(izoom-1, 2*iplot+1);
    }
}


// Helper called by _emit_plot() to upsample an rgb array along the time axis.
// Both arrays 'rgb_dst', 'rgb_src' have shapes (ny_arr, img_nx, 3).

void zoomable_tileset_state::_upsample_rgb(uint8_t *rgb_dst, const uint8_t *rgb_src, bool second_half)
{
    ssize_t m = ny_arr;
    ssize_t n = img_nx / 2;
    
    if (second_half)
	rgb_src += n;

    for (ssize_t i = 0; i < m; i++) {
	uint8_t *dst_row = rgb_dst + i * (6*n);
	const uint8_t *src_row = rgb_src + i * (6*n);

	for (ssize_t j = 0; j < n; j++) {
	    uint8_t r = src_row[3*j];
	    uint8_t g = src_row[3*j+1];
	    uint8_t b = src_row[3*j+2];

	    dst_row[6*j] = r;
	    dst_row[6*j+1] = g;
	    dst_row[6*j+2] = b;
	    dst_row[6*j+3] = r;
	    dst_row[6*j+4] = g;
	    dst_row[6*j+5] = b;
	}
    }
}


void zoomable_tileset_state::_initialize_json()
{
    json_output = Json::Value(Json::arrayValue);

    for (int i = 0; i < img_nzoom; i++) {
	Json::Value j;
	j["ny"] = Json::Int64(img_ny);
	j["nt_per_pix"] = Json::Int64(img_nds * (1 << i));
	j["files"] = Json::Value(Json::arrayValue);   // initially empty, but will be populated later

	json_output.append(j);
    }
}


}  // namespace rf_pipelines
