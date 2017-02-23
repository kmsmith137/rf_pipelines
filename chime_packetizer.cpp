#include <sstream>
#include "rf_pipelines_internals.hpp"

#ifdef HAVE_CH_FRB_IO
#include <ch_frb_io.hpp>
#endif

using namespace std;

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif


#ifndef HAVE_CH_FRB_IO

shared_ptr<wi_transform> make_chime_packetizer(const string &dstname, int nfreq_per_packet, int nt_per_chunk, int nt_per_packet, float wt_cutoff,
					       int beam_id)
{
    throw runtime_error("rf_pipelines::make_chime_packetizer() was called, but rf_pipelines was compiled without ch_frb_io");
}

#else  // HAVE_CH_FRB_IO


struct chime_packetizer : public wi_transform {
    // Note: inherits { nfreq, nt_chunk, nt_prepad, nt_postpad } from base class wi_transform

    ch_frb_io::intensity_network_ostream::initializer ini_params;
    uint64_t current_fpga_count = 0;

    shared_ptr<ch_frb_io::intensity_network_ostream> ostream;

    chime_packetizer(const string &dstname, int nfreq_per_packet, int nt_per_chunk, int nt_per_packet, float wt_cutoff, double target_gbps, int beam_id=0);

    virtual void set_stream(const wi_stream &stream);
    virtual void start_substream(int isubstream, double t0);
    virtual void process_chunk(double t0, double t1, float *intensity, float *weights, ssize_t stride, float *pp_intensity, float *pp_weights, ssize_t pp_stride);
    virtual void end_substream();
    virtual std::string get_name() const { return "chime_packetizer"; }
};


chime_packetizer::chime_packetizer(const string &dstname, int nfreq_coarse_per_packet, int nt_per_chunk, int nt_per_packet, float wt_cutoff, double target_gbps,
int beam_id)
{
    // Argument checking

    constexpr int nfreq_coarse = ch_frb_io::constants::nfreq_coarse_tot;

    if ((nfreq_coarse_per_packet <= 0) || (nfreq_coarse % nfreq_coarse_per_packet))
	throw runtime_error("chime_packetizer: currently nfreq_coarse_per_packet must be a divisor of " + to_string(nfreq_coarse));
    if (nt_per_chunk <= 0)
	throw runtime_error("chime_packetizer: nt_per_chunk > 0 is required");
    if (nt_per_packet <= 0)
	throw runtime_error("chime_packetizer: nt_per_packet > 0 is required");
    if (nt_per_chunk % nt_per_packet)
	throw runtime_error("chime_packetizer: nt_per_chunk must be a multiple of nt_per_packet");

    // Initialize ini_params (some initializations deferred to set_stream)

    this->ini_params.beam_ids = { beam_id };
    this->ini_params.dstname = dstname;
    this->ini_params.nfreq_coarse_per_packet = nfreq_coarse_per_packet;
    this->ini_params.nt_per_chunk = nt_per_chunk;
    this->ini_params.nt_per_packet = nt_per_packet;
    this->ini_params.wt_cutoff = wt_cutoff;
    this->ini_params.target_gbps = target_gbps;

    this->ini_params.coarse_freq_ids.resize(nfreq_coarse, -1);
    for (int i = 0; i < nfreq_coarse; i++)
	this->ini_params.coarse_freq_ids[i] = i;

    // Initialize base class members
    this->nt_chunk = nt_per_chunk;
    this->name = "chime_packetizer";
}


void chime_packetizer::set_stream(const wi_stream &stream)
{
    constexpr int nfreq_coarse = ch_frb_io::constants::nfreq_coarse_tot;

    if (stream.nfreq % nfreq_coarse)
	throw runtime_error("chime_packetizer: currently stream.nfreq must be a multiple of " + to_string(nfreq_coarse) + " (see comment in rf_pipelines.hpp)");

    this->nfreq = stream.nfreq;
    this->ini_params.nupfreq = stream.nfreq / nfreq_coarse;

    // infer fpga_counts_per_sample from stream.dt_sample
    double f = stream.dt_sample / constants::chime_seconds_per_fpga_count;
    this->ini_params.fpga_counts_per_sample = int(f+0.5);   // round to nearest integer

    if (fabs(f - ini_params.fpga_counts_per_sample) > 0.01) {
	// We use a stringstream here since to_string() gives a weird formatting
	stringstream ss;
	ss << "chime_packetizer: currently stream.dt_sample must be a multiple of " << constants::chime_seconds_per_fpga_count <<  " seconds (see comment in rf_pipelines.hpp)";
	throw runtime_error(ss.str());
    }

    this->ostream = ch_frb_io::intensity_network_ostream::make(ini_params);
}


void chime_packetizer::start_substream(int isubstream, double t0)
{
    if (isubstream > 0)
	throw runtime_error("chime_packetizer: multiple substreams are not currently supported");

    // infer current_fpga_count from t0
    if (t0 > 0.0) {
	double u = t0 / constants::chime_seconds_per_fpga_count / ini_params.fpga_counts_per_sample;
	this->current_fpga_count = int(u+0.5) * ini_params.fpga_counts_per_sample;
    }
}


void chime_packetizer::process_chunk(double t0, double t1, float *intensity, float *weights, ssize_t stride, float *pp_intensity, float *pp_weights, ssize_t pp_stride)
{
    this->ostream->send_chunk(intensity, weights, stride, current_fpga_count);
    this->current_fpga_count += nt_chunk * ini_params.fpga_counts_per_sample;
}


void chime_packetizer::end_substream()
{
    bool join_network_thread = true;
    ostream->end_stream(join_network_thread);
    ostream.reset();
}


shared_ptr<wi_transform> make_chime_packetizer(const string &dstname, int nfreq_per_packet, int nt_per_chunk, int nt_per_packet, float wt_cutoff, double target_gbps, int beam_id)
{
    return make_shared<chime_packetizer> (dstname, nfreq_per_packet, nt_per_chunk, nt_per_packet, wt_cutoff, target_gbps, beam_id);
}


#endif  // HAVE_CH_FRB_IO


}  // namespace rf_pipelines
