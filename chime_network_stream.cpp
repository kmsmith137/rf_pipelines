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

shared_ptr<wi_stream> make_chime_network_stream(int udp_port, int beam_id)
{
    throw runtime_error("rf_pipelines::make_chime_network_stream() was called, but rf_pipelines was compiled without ch_frb_io");
}

#else  // HAVE_CH_FRB_IO


struct chime_network_stream : public wi_stream
{
    shared_ptr<ch_frb_io::intensity_network_stream> stream;

    chime_network_stream(int udp_port, int beam_id);
    virtual ~chime_network_stream() { }

    virtual void stream_start();
    virtual void stream_body(wi_run_state &run_state);
};


chime_network_stream::chime_network_stream(int udp_port, int beam_id)
{
    ch_frb_io::intensity_network_stream::initializer ini_params;
    ini_params.beam_ids = { beam_id };

    if (udp_port > 0)
	ini_params.udp_port = udp_port;

    this->stream = ch_frb_io::intensity_network_stream::make(ini_params);
    
    // The frequency band is not currently part of the L0-L1 packet format, so we just use hardcoded values.
    // If we ever want to reuse the packet format in a different experiment, this should be fixed by updating the protocol.
    this->freq_lo_MHz = 400.0;
    this->freq_hi_MHz = 800.0;
    this->nt_maxwrite = ch_frb_io::constants::nt_per_assembled_chunk;
    
    // Initialization of base class members {nfreq, dt_sample} is deferred to stream_start().
}


void chime_network_stream::stream_start()
{
    // tells network thread to start reading packets, returns immediately
    stream->start_stream();
    
    int nupfreq;
    int nt_per_packet;
    uint64_t fpga_counts_per_sample;
    uint64_t fpga_count0;

    // block until first packet is received
    bool alive = stream->get_first_packet_params(nupfreq, nt_per_packet, fpga_counts_per_sample, fpga_count0);

    if (!alive)
	throw runtime_error("chime_network_stream: no packets received");

    // now we can initialize {nfreq, dt_sample}
    this->nfreq = ch_frb_io::constants::nfreq_coarse_tot * nupfreq;
    this->dt_sample = constants::chime_seconds_per_fpga_count * fpga_counts_per_sample;
}


void chime_network_stream::stream_body(wi_run_state &run_state)
{
    bool startflag = false;

    for (;;) {
	shared_ptr<ch_frb_io::assembled_chunk> chunk = stream->get_assembled_chunk(0);

	if (!chunk)
	    break;

	rf_assert(this->nfreq == ch_frb_io::constants::nfreq_coarse_tot * chunk->nupfreq);

	double t0 = chunk->isample * chunk->fpga_counts_per_sample * constants::chime_seconds_per_fpga_count;

	if (!startflag) {
	    run_state.start_substream(t0);
	    startflag = true;
	}

	float *dst_intensity = nullptr;
	float *dst_weights = nullptr;
	ssize_t dst_stride = 0;
	bool zero_flag = false;

	run_state.setup_write(nt_maxwrite, dst_intensity, dst_weights, dst_stride, zero_flag, t0);
	chunk->decode(dst_intensity, dst_weights, dst_stride);
	run_state.finalize_write(nt_maxwrite);
    }
    
    stream->join_threads();
    run_state.end_substream();
}


shared_ptr<wi_stream> make_chime_network_stream(int udp_port, int beam_id)
{
    shared_ptr<chime_network_stream> ret = make_shared<chime_network_stream> (udp_port, beam_id);
    return ret;
}


#endif  // HAVE_CH_FRB_IO

}   // namespace rf_pipelines
