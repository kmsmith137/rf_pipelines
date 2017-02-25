#include "rf_pipelines_internals.hpp"

#ifdef HAVE_CH_FRB_IO
#include <ch_frb_io.hpp>
#endif

namespace rf_pipelines {
#if 0
};  // pacify emacs c-mode!
#endif


class chime_file_stream_base : public wi_stream {
protected:
    // Specified at construction.  For an explanation of the 'noise_source_align' field see rf_pipelines.hpp.
    // Note that the 'nt_chunk' constructor argument is used to initialize the base class member wi_stream::nt_maxwrite.
    std::vector<std::string> filename_list;
    ssize_t noise_source_align;   // if zero, no alignment will be performed

    int curr_ifile = -1;  // index of current file in filename_list

    // Only nonzero if noise source alignment is requested
    ssize_t samples_per_noise_switch = 0;
    ssize_t initial_discard_count = 0;

    double time_lo;
    double time_hi;
    ssize_t nt;
    bool frequencies_are_increasing;

public:
    chime_file_stream_base(const std::vector<std::string> &filename_list_, ssize_t nt_chunk, ssize_t noise_source_align);
    virtual ~chime_file_stream_base() { }

    virtual void stream_start();
    virtual void stream_body(wi_run_state &run_state);

protected:
    virtual void load_file(const std::string& filename);
    virtual void close_file();
    virtual void set_params_from_file();
    virtual void check_file_consistency();
    virtual void read_data(float* dst_int, float* dst_wt, ssize_t it_file, ssize_t n, ssize_t dst_stride);

};














}   // namespace rf_pipelines
