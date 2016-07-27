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
    shared_ptr<ch_frb_io::intensity_beam_assembler> assembler;

    chime_network_stream(int udp_port, int beam_id);
    virtual ~chime_network_stream() { }

    virtual void stream_start();
    virtual void stream_body(wi_run_state &run_state);
};


chime_network_stream::chime_network_stream(int udp_port, int beam_id)
{
    this->assembler = ch_frb_io::intensity_beam_assembler::make(beam_id);
    this->stream = ch_frb_io::intensity_network_stream::make({assembler}, udp_port);
    
    // The frequency band is not currently part of the L0-L1 packet format, so we just use hardcoded values.
    // If we ever want to reuse the packet format in a different experiment, this should be fixed by updating the protocol.
    this->freq_lo_MHz = 400.0;
    this->freq_hi_MHz = 800.0;
    this->nt_maxwrite = ch_frb_io::assembled_chunk::nt_per_chunk;
    
    // Initialization of base class members {nfreq, dt_sample} is deferred to stream_start().
}


void chime_network_stream::stream_start()
{
    // tells network thread to start reading packets, returns immediately
    stream->start_stream();
    
    // block until first packet is received
    int fpga_counts_per_sample, nupfreq;
    assembler->wait_for_stream_params(fpga_counts_per_sample, nupfreq);

    // now we can initialize {nfreq, dt_sample}
    this->nfreq = 1024 * nupfreq;                       // FIXME hardcoded 1024
    this->dt_sample = 2.5e-6 * fpga_counts_per_sample;  // FIXME hardcoded 2.5e-6
}


void chime_network_stream::stream_body(wi_run_state &run_state)
{
    bool startflag = false;
    shared_ptr<ch_frb_io::assembled_chunk> chunk;

    for (;;) {
	bool alive = assembler->get_assembled_chunk(chunk);
	if (!alive)
	    break;

	rf_assert(nfreq == 1024 * chunk->nupfreq);   // FIXME hardcoded 1024

	double t0 = chunk->chunk_t0 * chunk->fpga_counts_per_sample * 2.5e-6;   // FIXME hardcoded 2.5e-6
	const float *src_intensity = chunk->intensity;
	const float *src_weights = chunk->weights;

	if (!startflag) {
	    run_state.start_substream(t0);
	    startflag = true;
	}

	float *dst_intensity = nullptr;
	float *dst_weights = nullptr;
	ssize_t dst_stride = 0;
	bool zero_flag = false;

	run_state.setup_write(nt_maxwrite, dst_intensity, dst_weights, dst_stride, zero_flag, t0);

	for (int ifreq = 0; ifreq < nfreq; ifreq++) {
	    memcpy(dst_intensity + ifreq*dst_stride, src_intensity + ifreq*nt_maxwrite, nt_maxwrite * sizeof(float));
	    memcpy(dst_weights + ifreq*dst_stride, src_weights + ifreq*nt_maxwrite, nt_maxwrite * sizeof(float));
	}

	run_state.finalize_write(nt_maxwrite);
	chunk.reset();
    }
    
    // "true" joins both the network and assembler threads
    stream->wait_for_end_of_stream(true);

    if (!startflag)
	throw runtime_error("chime_network_stream: no packets received");

    run_state.end_substream();
}


shared_ptr<wi_stream> make_chime_network_stream(int udp_port, int beam_id)
{
    shared_ptr<chime_network_stream> ret = make_shared<chime_network_stream> (udp_port, beam_id);
    return ret;
}


#endif  // HAVE_CH_FRB_IO

}   // namespace rf_pipelines
