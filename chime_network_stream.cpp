#include <algorithm>
#include "rf_pipelines_internals.hpp"

#ifdef HAVE_CH_FRB_IO
#include <ch_frb_io.hpp>
#endif

using namespace std;

namespace rf_pipelines {
#if 0
};  // pacify emacs c-mode!
#endif


#ifndef HAVE_CH_FRB_IO

shared_ptr<wi_stream> make_chime_network_stream(const shared_ptr<ch_frb_io::intensity_network_stream> &stream, int assembler_id)
{
    throw runtime_error("rf_pipelines::make_chime_network_stream() was called, but rf_pipelines was compiled without ch_frb_io");
}

shared_ptr<wi_stream> make_chime_network_stream(int udp_port, int beam_id)
{
    throw runtime_error("rf_pipelines::make_chime_network_stream() was called, but rf_pipelines was compiled without ch_frb_io");
}

#else  // HAVE_CH_FRB_IO


struct chime_network_stream : public wi_stream
{
    shared_ptr<ch_frb_io::intensity_network_stream> stream;
    const int beam_id;

    int assembler_id = -1;

    // FIXME these fields will move into a json object soon (I think!)
    double dt_sample = 0.0;
    double freq_lo_MHz = 0.0;
    double freq_hi_MHz = 0.0;

    // FIXME this is a hack that should be removed, see below.
    shared_ptr<ch_frb_io::assembled_chunk> first_chunk;

    chime_network_stream(const shared_ptr<ch_frb_io::intensity_network_stream> &stream_, int beam_id_);
    virtual ~chime_network_stream() { }

    virtual bool _fill_chunk(float *intensity, ssize_t istride, float *weights, ssize_t wstride, ssize_t pos) override;
    virtual void _start_pipeline(Json::Value &j) override;
    virtual void _end_pipeline(Json::Value &j) override;
};


chime_network_stream::chime_network_stream(const shared_ptr<ch_frb_io::intensity_network_stream> &stream_, int beam_id_) :
    wi_stream("chime_network_stream"),
    stream(stream_), 
    beam_id(beam_id_)
{ 
    if (!stream)
	throw runtime_error("rf_pipelines: empty stream pointer passed to chime_network_stream constructor");

    const vector<int> &beam_ids = stream->ini_params.beam_ids;
    
    for (unsigned int i = 0; i < beam_ids.size(); i++) {
	if (beam_ids[i] == beam_id) {
	    this->assembler_id = i;
	    break;
	}
    }

    if (assembler_id < 0) {
	throw runtime_error("chime_network_stream constructor: beam_id=" + stringify(beam_id)
			    + " not found in stream beam_id list " + stringify(beam_ids));
    }

    this->nfreq = ch_frb_io::constants::nfreq_coarse_tot * stream->ini_params.nupfreq;
    this->nt_chunk = ch_frb_io::constants::nt_per_assembled_chunk;
    
    // The frequency band is not currently part of the L0-L1 packet format, so we just use hardcoded values.
    // If we ever want to reuse the packet format in a different experiment, this should be fixed by updating the protocol.
    this->freq_lo_MHz = 400.0;
    this->freq_hi_MHz = 800.0;
    this->dt_sample = chime_file_stream_base::chime_seconds_per_fpga_count * stream->ini_params.fpga_counts_per_sample;
}


void chime_network_stream::_start_pipeline(Json::Value &j)
{
    // tells network thread to start reading packets, returns immediately
    stream->start_stream();

    // FIXME this is awkward: we want to initialize 'initial_fpga_count' in _start_pipeline(), but the
    // only way to do this is by getting the first chunk, which we save until _fill_chunk() is called later.
    // This could be improved by defining a member function of ch_frb_io::intensity_network_stream which
    // returns the initial FPGA count as soon as the first packet is received (rather than waiting until the
    // first chunk is finalized).

    this->first_chunk = stream->get_assembled_chunk(assembler_id);

    if (!first_chunk)
	_throw("no data could be processed");

    uint64_t fpga0 = uint64_t(first_chunk->isample) * uint64_t(first_chunk->fpga_counts_per_sample);
    
    j["initial_fpga_count"] = Json::UInt64(fpga0);
    j["fpga_counts_per_sample"] = first_chunk->fpga_counts_per_sample;
}


bool chime_network_stream::_fill_chunk(float *intensity, ssize_t istride, float *weights, ssize_t wstride, ssize_t pos)
{
    if (first_chunk) {
	first_chunk->decode(intensity, weights, istride, wstride);
	first_chunk.reset();
	return true;
    }

    shared_ptr<ch_frb_io::assembled_chunk> chunk = stream->get_assembled_chunk(assembler_id);

    if (!chunk)
	return false;

    rf_assert(this->nfreq == ch_frb_io::constants::nfreq_coarse_tot * chunk->nupfreq);
    chunk->decode(intensity, weights, istride, wstride);
    return true;
}


void chime_network_stream::_end_pipeline(Json::Value &j)
{
    first_chunk.reset();
    stream->join_threads();
}


// -------------------------------------------------------------------------------------------------


shared_ptr<wi_stream> make_chime_network_stream(const shared_ptr<ch_frb_io::intensity_network_stream> &stream, int beam_id)
{
    return make_shared<chime_network_stream> (stream, beam_id);
}


shared_ptr<wi_stream> make_chime_network_stream(int udp_port, int beam_id)
{
    ch_frb_io::intensity_network_stream::initializer ini_params;
    ini_params.beam_ids = { beam_id };

    if (udp_port > 0)
	ini_params.udp_port = udp_port;

    auto stream = ch_frb_io::intensity_network_stream::make(ini_params);
    return make_chime_network_stream(stream, beam_id);
}


#endif  // HAVE_CH_FRB_IO

}   // namespace rf_pipelines
