#include <cassert>
#include "rf_pipelines_internals.hpp"
#include "rf_pipelines_inventory.hpp"

#ifdef HAVE_CH_FRB_IO
#include <ch_frb_io.hpp>
#else
namespace ch_frb_io { class assembled_chunk; }
#endif

using namespace std;

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif

void chime_wi_transform::set_chime_stream(std::shared_ptr<ch_frb_io::intensity_network_stream> stream,
                                          int beam_id) {
    if (this->state != UNBOUND)
	throw runtime_error("chime_wi_transform::set_chime_stream() called after bind()");
    if (stream && (beam_id < 0))
	throw runtime_error("chime_wi_transform::set_chime_stream(): chime_stream was specified, but chime_beam_id was uninitialized or negative");
#ifndef HAVE_CH_FRB_IO
    if (stream)
	throw runtime_error("chime_wi_transform::set_chime_stream(): chime_stream was specified, but rf_pipelines was compiled without ch_frb_io!");
#endif

    this->chime_stream = stream;
    this->chime_beam_id = beam_id;
}

void chime_wi_transform::_bind_transform(Json::Value &json_attrs) {
    if (!json_attrs.isMember("freq_lo_MHz") || !json_attrs.isMember("freq_hi_MHz"))
	throw runtime_error("chime_wi_transform: expected json_attrs to contain members 'freq_lo_MHz' and 'freq_hi_MHz'");
    if (!json_attrs.isMember("dt_sample"))
	throw runtime_error("chime_wi_transform: expected json_attrs to contain member 'dt_sample'");
    double freq_lo_MHz = json_attrs["freq_lo_MHz"].asDouble();
    double freq_hi_MHz = json_attrs["freq_hi_MHz"].asDouble();
    double dt_sample = json_attrs["dt_sample"].asDouble();
}

void chime_wi_transform::_start_pipeline(Json::Value &json_attrs)
{
#ifdef HAVE_CH_FRB_IO
    if (chime_stream) {
	this->chime_initial_fpga_count = uint64_t_from_json(json_attrs, "initial_fpga_count");
	this->chime_fpga_counts_per_sample = int_from_json(json_attrs, "fpga_counts_per_sample");
	this->chime_fpga_counts_initialized = true;
	
	if (chime_stream->ini_params.fpga_counts_per_sample != chime_fpga_counts_per_sample)
	    throw runtime_error("chime_wi_transform: value of 'fpga_counts_per_sample' in chime_intensity_stream does not match the value in _start_pipeline()");
    }
#endif
}

std::shared_ptr<ch_frb_io::assembled_chunk> chime_wi_transform::assembled_chunk_for_pos(ssize_t pos) {
    shared_ptr<ch_frb_io::assembled_chunk> chunk;
    if (chime_stream) {
	if (!chime_fpga_counts_initialized)
	    throw runtime_error("rf_pipelines::chime_mask_counter internal error: fpga count fields were not initialized as expected");
    
	// The 'pos' argument is the current pipeline position in units of time samples -- convert to FPGA counts
	uint64_t fpga_counts = pos * this->chime_fpga_counts_per_sample + this->chime_initial_fpga_count;

	// The last argument in find_assembled_chunk() is 'toplevel'.
	chunk = chime_stream->find_assembled_chunk(chime_beam_id, fpga_counts, true);
    }
    return chunk;
}

} // namespace




