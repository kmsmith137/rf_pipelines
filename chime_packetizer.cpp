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

shared_ptr<wi_transform> make_chime_packetizer(const string &dstname, int nfreq_per_packet, int nt_per_chunk, int nt_per_packet, float wt_cutoff)
{
    throw runtime_error("rf_pipelines::make_chime_packetizer() was called, but rf_pipelines was compiled without ch_frb_io");
}

#else  // HAVE_CH_FRB_IO


struct chime_packetizer : public wi_transform {
    // Note: inherits { nfreq, nt_chunk, nt_prepad, nt_postpad } from base class wi_transform

    string dstname;
    int nfreq_per_packet;
    int nt_per_packet;
    float wt_cutoff;
    double target_gbps;

    int nupfreq;
    int fpga_counts_per_sample;
    uint64_t current_fpga_count = 0;

    shared_ptr<ch_frb_io::intensity_network_ostream> stream_obj;

    chime_packetizer(const string &dstname, int nfreq_per_packet, int nt_per_chunk, int nt_per_packet, float wt_cutoff, double target_gbps);

    virtual void set_stream(const wi_stream &stream);
    virtual void start_substream(int isubstream, double t0);
    virtual void process_chunk(double t0, double t1, float *intensity, float *weights, ssize_t stride, float *pp_intensity, float *pp_weights, ssize_t pp_stride);
    virtual void end_substream();
    virtual std::string get_name() const { return "chime_packetizer"; }
};


chime_packetizer::chime_packetizer(const string &dstname_, int nfreq_per_packet_, int nt_per_chunk, int nt_per_packet_, float wt_cutoff_, double target_gbps_)
{
    constexpr int nfreq_coarse = ch_frb_io::constants::nfreq_coarse;

    this->nt_chunk = nt_per_chunk;
    this->name = "chime_packetizer";

    this->dstname = dstname_;
    this->nfreq_per_packet = nt_per_packet_;
    this->nt_per_packet = nt_per_packet_;
    this->wt_cutoff = wt_cutoff_;
    this->target_gbps = target_gbps_;
    
    // initialization of the following members is deferred to set_stream() or start_substream()
    //   nfreq, nupfreq, fpga_counts_per_sample, current_fpga_count

    if ((nfreq_per_packet <= 0) || (nfreq_coarse % nfreq_per_packet))
	throw runtime_error("chime_packetizer: currently nfreq_per_packet must be a divisor of " + to_string(nfreq_coarse));
    if (nt_per_chunk <= 0)
	throw runtime_error("chime_packetizer: nt_per_chunk > 0 is required");
    if (nt_per_packet <= 0)
	throw runtime_error("chime_packetizer: nt_per_packet > 0 is required");
    if (nt_per_chunk % nt_per_packet)
	throw runtime_error("chime_packetizer: nt_per_chunk must be a multiple of nt_per_packet");
}


void chime_packetizer::set_stream(const wi_stream &stream)
{
    constexpr int nfreq_coarse = ch_frb_io::constants::nfreq_coarse;

    if (stream.nfreq % nfreq_coarse)
	throw runtime_error("chime_packetizer: currently stream.nfreq must be a multiple of " + to_string(nfreq_coarse));

    this->nfreq = stream.nfreq;
    this->nupfreq = stream.nfreq / nfreq_coarse;

    // infer fpga_counts_per_sample from stream.dt_sample
    double f = stream.dt_sample / constants::chime_seconds_per_fpga_count;
    this->fpga_counts_per_sample = int(f+0.5);   // round to nearest integer

    if (fabs(f - fpga_counts_per_sample) > 0.01)
	throw runtime_error("chime_packetizer: currently stream.dt_sample must be a multiple of " + to_string(constants::chime_seconds_per_fpga_count) + " seconds");

    vector<int> ibeam = { 0 };

    vector<int> ifreq_chunk(nfreq);
    for (int i = 0; i < nfreq; i++)
	ifreq_chunk[i] = i;

    this->stream_obj = ch_frb_io::intensity_network_ostream::make(dstname, ibeam, ifreq_chunk, nupfreq, nt_chunk, nfreq_per_packet, nt_per_packet, fpga_counts_per_sample, wt_cutoff, target_gbps);
}


void chime_packetizer::start_substream(int isubstream, double t0)
{
    if (isubstream > 0)
	throw runtime_error("chime_packetizer: multiple substreams are not currently supported");

    // infer current_fpga_count from t0
    if (t0 > 0.0) {
	double u = t0 / constants::chime_seconds_per_fpga_count / fpga_counts_per_sample;
	this->current_fpga_count = int(u+0.5) * fpga_counts_per_sample;
    }
}


void chime_packetizer::process_chunk(double t0, double t1, float *intensity, float *weights, ssize_t stride, float *pp_intensity, float *pp_weights, ssize_t pp_stride)
{
    this->stream_obj->send_chunk(intensity, weights, stride, current_fpga_count);
    this->current_fpga_count += nt_chunk * fpga_counts_per_sample;
}


void chime_packetizer::end_substream()
{
    bool join_network_thread = true;
    stream_obj->end_stream(join_network_thread);
    stream_obj.reset();
}


shared_ptr<wi_transform> make_chime_packetizer(const string &dstname, int nfreq_per_packet, int nt_per_chunk, int nt_per_packet, float wt_cutoff, double target_gbps)
{
    return make_shared<chime_packetizer> (dstname, nfreq_per_packet, nt_per_chunk, nt_per_packet, wt_cutoff, target_gbps);
}


#endif  // HAVE_CH_FRB_IO


}  // namespace rf_pipelines
