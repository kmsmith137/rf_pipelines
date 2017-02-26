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

    // Note: inherits members { nfreq, freq_lo_MHz, freq_hi_MHz, dt_sample } from base class wi_stream.
    double time_lo = 0.0;
    double time_hi = 0.0;
    ssize_t nt = 0;
    bool frequencies_are_increasing = false;

public:
    chime_file_stream_base(const std::vector<std::string> &filename_list_, ssize_t nt_chunk, ssize_t noise_source_align);
    virtual ~chime_file_stream_base() { }

    virtual void stream_start();
    virtual void stream_body(wi_run_state &run_state);

protected:
    // Reads file from disk into a file-type-specific internal data structure.
    // This will be followed by calls to set_params_from_file() and/or check_file_consistency().
    virtual void load_file(const std::string& filename) = 0;

    // Responsible for initializing wi_stream members: { nfreq, freq_lo_MHz, freq_hi_MHz, dt_sample }
    // and chime_file_stream_base members: { time_lo, time_hi, nt, frequencies_are_increasing }
    virtual void set_params_from_file() = 0;
    
    // Checks consistency between file and data members mentioned above, but doesn't initialize them.
    virtual void check_file_consistency() const = 0;

    // Reads a shape-(nfreq,n) block from current file.
    // The 'dst_int' and 'dst_wt' arrays have shape (nfreq, n) and stride 'dst_stride'.
    // The 'it_file' argument is the initial timestamp of the block, relative to the start of the current file.
    //
    // The subclass can choose whether the 'dst_int' and 'dst_wt" arrays are indexed from lowest frequency 
    // channel to highest, or from highest frequency channel to lowest (the rf_pipelines_default).  This is 
    // done by setting the 'frequencies_are_increasing' flag in set_params_from_file().  The way it is 
    // implemented is via logic in the base class which flips the sign of 'dst_stride' if necessary.

    virtual void read_data(float* dst_int, float* dst_wt, ssize_t it_file, ssize_t n, ssize_t dst_stride) const = 0;

    virtual void close_file() = 0;
};


}   // namespace rf_pipelines
