#ifndef RF_PIPELINE_CHIME_PACKETIZER_H
#define RF_PIPELINE_CHIME_PACKETIZER_H

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif

#ifndef HAVE_CH_FRB_IO
//
#else  // HAVE_CH_FRB_IO

#include <ch_frb_io.hpp>

struct chime_packetizer : public wi_transform {
    // Note: inherits { nfreq, nt_chunk, nt_prepad, nt_postpad } from base class wi_transform
    ch_frb_io::intensity_network_ostream::initializer ini_params;
    uint64_t current_fpga_count = 0;

    std::shared_ptr<ch_frb_io::intensity_network_ostream> ostream;

    chime_packetizer(const std::string &dstname, int nfreq_per_packet, int nt_per_chunk, int nt_per_packet, float wt_cutoff, double target_gbps, int beam_id=0);

    virtual void set_stream(const wi_stream &stream);
    virtual void start_substream(int isubstream, double t0);
    virtual void process_chunk(double t0, double t1, float *intensity, float *weights, ssize_t stride, float *pp_intensity, float *pp_weights, ssize_t pp_stride);
    virtual void end_substream();
    virtual std::string get_name() const { return "chime_packetizer"; }
};

#endif  // HAVE_CH_FRB_IO


}  // namespace rf_pipelines

#endif
