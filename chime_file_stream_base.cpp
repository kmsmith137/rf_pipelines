#include <algorithm>
#include "rf_pipelines_internals.hpp"
#include "chime_file_stream_base.hpp"

#ifdef HAVE_CH_FRB_IO
#include <ch_frb_io.hpp>
#endif

using namespace std;

namespace rf_pipelines {
#if 0
};  // pacify emacs c-mode!
#endif

    
chime_file_stream_base::chime_file_stream_base(const vector<string> &filename_list_, ssize_t nt_chunk, ssize_t noise_source_align_) :
    filename_list(filename_list_), noise_source_align(noise_source_align_)
{
    rf_assert(filename_list.size() > 0);
    rf_assert(noise_source_align >= 0);

    // We defer initialization of the base class members { nfreq, freq_lo_MHz, freq_hi_MHz, dt_sample } to stream_start().
    this->nt_maxwrite = (nt_chunk > 0) ? nt_chunk : 1024;
}

// virtual
void chime_file_stream_base::stream_start()
{
    load_file(filename_list[0]);
    this->curr_ifile = 0;

    set_params_from_file();

    if (noise_source_align > 0) {
	if (dt_sample < 1.0e-4)
	    throw std::runtime_error("chime_file_stream_base constructor: unexpectedly small dt_sample");

	// The CHIME noise source switch occurs at FPGA counts which are multiples of 2^23.
	double fsamp = (1 << 23) * constants::chime_seconds_per_fpga_count / dt_sample;
	this->samples_per_noise_switch = ssize_t(fsamp + 0.5);

	if (fabs(samples_per_noise_switch - fsamp) > 1.0e-4)
	    throw std::runtime_error("chime_file_stream_base constructor: dt_sample does not appear to evenly divide dt_noise_source");
	if (samples_per_noise_switch % noise_source_align != 0)
	    throw std::runtime_error("chime_file_stream_base constructor: 'noise_source_align' must evenly divide " + to_string(samples_per_noise_switch));

	double fi0 = time_lo / dt_sample;
	ssize_t i0 = ssize_t(fi0 + 0.5);

	if (fabs(fi0 - i0) > 1.0e-4)
	    throw std::runtime_error("chime_file_stream_base constructor: file timestamp does not appear to evenly divide dt_sample");

	i0 = i0 % noise_source_align;
	this->initial_discard_count = (noise_source_align - i0) % noise_source_align;
    }
}

    
// virtual
void chime_file_stream_base::stream_body(wi_run_state &run_state)
{
    ssize_t it_file = 0;                         // time index in current file.  Can be negative!  This means there is a time gap between files.
    ssize_t it_chunk = -initial_discard_count;   // index in current chunk.  Can be negative!  This means we're still discarding initial samples.
    double stream_t0 = time_lo + initial_discard_count * dt_sample;
    int nfiles = filename_list.size();

    float *intensity;
    float *weights;
    ssize_t stride;
    bool zero_flag = true;

    // Reminder: setup_write() sets the 'intensity' and 'weights' buffers to zero.
    // Therefore, if part of the buffer is unmodified here, it will be treated as masked (weight-zero).
    run_state.start_substream(stream_t0);
    run_state.setup_write(nt_maxwrite, intensity, weights, stride, zero_flag, stream_t0);

    for (;;) {
	// At top of loop, the state of the stream is represented by the following variables;
	//   curr_ifile  -> index of current file
	//   it_file     -> time index within current file
	//   it_chunk    -> time index within output chunk

	if (curr_ifile >= nfiles) {
	    // End of stream
	    // FIXME should be able to write fewer than nt_maxwrite samples here
	    run_state.finalize_write(nt_maxwrite);
	    run_state.end_substream();

            close_file();
	    this->curr_ifile = -1;
	    return;
	}
	else if (it_file >= nt) {
	    // End of file.
	    curr_ifile++;
	    it_file = 0;

	    if (curr_ifile >= nfiles) {
                close_file();
		continue;
	    }
	    
	    // Open next file and do consistency tests.
	    double old_t1 = time_hi;
            load_file(filename_list[curr_ifile]);

            check_file_consistency();

            // pick up new time_lo, etc
            set_params_from_file();

	    double new_t0 = time_lo;
	    double gap = (new_t0 - old_t1) / dt_sample;  // time gap between files, as multiple of dt_sample

	    // FIXME another option here would be to start a new substream if a large gap is encountered.
	    if (gap < -0.01)
		throw runtime_error("chime_file_stream_base: acquisition files are not in time order?!");
	    if (gap > 10000.01)
		throw runtime_error("chime_file_stream_base: excessively long gap in acquisition");

	    ssize_t ngap = (ssize_t)(gap + 0.5);  // round to integer
	    rf_assert(ngap >= 0);

	    if (fabs(gap-ngap) > 0.01)
		throw runtime_error("chime_file_stream_base: time gap between files is not an integer multiple of dt_sample");
	    
	    it_file = -ngap;
	}
	else if (it_chunk >= nt_maxwrite) {
	    // Write chunk
	    run_state.finalize_write(nt_maxwrite);
	    run_state.setup_write(nt_maxwrite, intensity, weights, stride, zero_flag, time_lo + it_file * dt_sample);
	    it_chunk = 0;
	}
	else if (it_file < 0) {
	    // Skip gap between files.
	    ssize_t n = min(-it_file, nt_maxwrite - it_chunk);
	    it_file += n;
	    it_chunk += n;
	}
	else if (it_chunk < 0) {
	    // Drop data, by advancing it_chunk.
	    // This only happens if 'noise_source_align' is specified, and we need to drop the first few samples in the stream.
	    ssize_t n = min(-it_chunk, nt - it_file);
	    it_file += n;
	    it_chunk += n;
	}
	else {
	    // Read data from file
	    ssize_t n = min(nt - it_file, nt_maxwrite - it_chunk);
	    
	    // A note on frequency channel ordering.  In rf_pipelines, frequencies must 
	    // be ordered from highest frequency to lowest.  In the ch_frb_io::intensity_hdf5_file,
	    // either ordering can be used depending on the value of the 'frequencies_are_increasing' 
	    // flag.  If this flag is set, then we need to reorder by using a negative stride.

	    bool incflag = frequencies_are_increasing;
	    float *dst_int = incflag ? (intensity + (nfreq-1)*stride + it_chunk) : (intensity + it_chunk);
	    float *dst_wt = incflag ? (weights + (nfreq-1)*stride + it_chunk) : (weights + it_chunk);
	    ssize_t dst_stride = incflag ? (-stride) : stride;

            read_data(dst_int, dst_wt, it_file, n, dst_stride);

	    it_file += n;
	    it_chunk += n;	    
	}
    }
}



}   // namespace rf_pipelines
