#include <algorithm>
#include "rf_pipelines_internals.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
};  // pacify emacs c-mode!
#endif

    
chime_file_stream_base::chime_file_stream_base(const string &class_name, const vector<string> &filename_list_, ssize_t nt_chunk_, ssize_t noise_source_align_) :
    wi_stream(class_name),
    filename_list(filename_list_),
    noise_source_align(noise_source_align_)
{
    rf_assert(filename_list.size() > 0);
    rf_assert(noise_source_align >= 0);

    // We initialize nt_chunk here, but defer initialization of nfreq to chime_file_stream_base::_bind_stream().
    this->nt_chunk = (nt_chunk_ > 0) ? nt_chunk_ : 1024;
}


// Virtual override
void chime_file_stream_base::_bind_stream(Json::Value &json_attrs)
{
    load_file(filename_list[0]);
    this->curr_ifile = 0;

    // Initializes the following fields: nfreq, freq_lo_MHz, freq_hi_MHz, dt_sample, time_lo, time_hi, nt, frequencies_are_increasing.
    this->set_params_from_file();

    double dt = time_lo / dt_sample;
    double df = dt_sample / chime_seconds_per_fpga_count;

    this->initial_sample_index = ssize_t(dt + 0.5);
    this->fpga_counts_per_sample = ssize_t(df + 0.5);

    if (dt_sample < 1.0e-4)
	_throw("chime_file_stream_base constructor: unexpectedly small dt_sample");
    if (dt_sample > 0.01)
	_throw("chime_file_stream_base constructor: unexpectedly large dt_sample");
    if (time_lo < 0.0)
	_throw("chime_file_stream_base constructor: expected time_lo >= 0");
    if (time_lo > 1.0e10 * dt_sample)
	_throw("chime_file_stream_base constructor: unexpectedly large time_lo");
    if (fabs(df - fpga_counts_per_sample) > 1.0e-4)
	_throw("chime_file_stream_base constructor: dt_sample does not appear to be a multiple of 2.56e-6 secconds (one FPGA count)");
    if (fabs(dt - initial_sample_index) > 0.01)
	_throw("chime_file_stream_base constructor: initial timestamp does not appear to be a multiple of dt_sample");

    // Initialize initial_discard_count.
    if (noise_source_align > 0) {
	// The CHIME noise source switch occurs at FPGA counts which are multiples of 2^23.
	ssize_t fpga_counts_per_noise_switch = 1 << 23;
	if (fpga_counts_per_noise_switch % fpga_counts_per_sample)
	    _throw("chime_file_stream_base_constructor: expected fpga_counts_per_sample to evenly divide 2^23");

	ssize_t samples_per_noise_switch = xdiv(fpga_counts_per_noise_switch, fpga_counts_per_sample);
	if (samples_per_noise_switch % noise_source_align)
	    _throw("chime_file_stream_base constructor: 'noise_source_align' must evenly divide " + to_string(samples_per_noise_switch));

	ssize_t i0 = initial_sample_index % noise_source_align;
	this->initial_discard_count = (noise_source_align - i0) % noise_source_align;
	this->initial_sample_index += initial_discard_count;
    }

    json_attrs["freq_lo_MHz"] = this->freq_lo_MHz;
    json_attrs["freq_hi_MHz"] = this->freq_hi_MHz;
    json_attrs["dt_sample"] = this->dt_sample;
    json_attrs["t_initial"] = this->time_lo + initial_discard_count * dt_sample;
}


// Virtual override
void chime_file_stream_base::_start_pipeline(Json::Value &json_attrs)
{
    // FIXME I plan to switch to a single set of json attrs soon.
    // Then we can move these two lines of code to _bind_stream(), and get rid of _start_pipeline().

    json_attrs["initial_fpga_count"] = Json::Int64(this->initial_sample_index * this->fpga_counts_per_sample);
    json_attrs["fpga_counts_per_sample"] = Json::Int64(this->fpga_counts_per_sample);
}
    
// virtual
bool chime_file_stream_base::_fill_chunk(float *intensity, ssize_t istride, float *weights, ssize_t wstride, ssize_t pos)
{
    int nfiles = filename_list.size();
    ssize_t it_chunk = 0;  // Index in current chunk.

    // We set the 'intensity' and 'weights' buffers to zero here.
    // Therefore, if the buffers (or a subset thereof) are unmodified, they will be treated as masked (weight-zero).

    for (int ifreq = 0; ifreq < nfreq; ifreq++) {
	memset(intensity + ifreq*istride, 0, nt_chunk * sizeof(float));
	memset(weights + ifreq*wstride, 0, nt_chunk * sizeof(float));
    }

    for (;;) {

	// At top of loop, the state of the stream is represented by the following variables:
	//   curr_ifile  -> index of current file
	//   it_file     -> time index within current file
	//   nt_file     -> number of time samples in current file
	//   it_chunk    -> time index within output chunk
	
	if (curr_ifile >= nfiles)
	    return false;

	if (it_file >= nt_file) {
	    // End of file.
	    curr_ifile++;
	    close_file();

	    if (curr_ifile >= nfiles)
		return false;
	    
	    // Open next file and do consistency tests.
	    double old_t1 = time_hi;
            load_file(filename_list[curr_ifile]);

            check_file_consistency();

	    // Read new time_lo, time_hi, nt_file.
            set_params_from_file();

	    // Infer gap length from timestamps.
	    double new_t0 = time_lo;
	    double gap = (new_t0 - old_t1) / dt_sample;  // time gap between files, as multiple of dt_sample

	    if (gap < -0.01)
		throw runtime_error("chime_file_stream_base: acquisition files are not in time order?!");
	    if (gap > 10000.01)
		throw runtime_error("chime_file_stream_base: excessively long gap in acquisition");

	    ssize_t ngap = (ssize_t)(gap + 0.5);  // round to integer
	    rf_assert(ngap >= 0);

	    if (fabs(gap-ngap) > 0.01)
		throw runtime_error("chime_file_stream_base: time gap between files is not an integer multiple of dt_sample");

	    // Note that it_file is negative if there is a gap between files.
	    it_file = -ngap;
	    continue;
	}

	if (it_chunk >= nt_chunk)
	    return true;  // End of chunk, but not end of file.

	if (initial_discard_count > 0) {
	    // Still discarding initial samples (only happens if 'noise_source_align' is specified)
	    ssize_t n = min(initial_discard_count, nt_file - it_file);
	    it_file += n;
	    initial_discard_count -= n;
	    continue;
	}

	if (it_file < 0) {
	    // Skip gap between files.
	    ssize_t n = min(-it_file, nt_chunk - it_chunk);
	    it_file += n;
	    it_chunk += n;
	    continue;
	}

	// If we get here, we will read data from file into the chunk.

	ssize_t n = min(nt_file - it_file, nt_chunk - it_chunk);
	    
	// A note on frequency channel ordering.  In rf_pipelines, frequencies must 
	// be ordered from highest frequency to lowest.  In the ch_frb_io::intensity_hdf5_file,
	// either ordering can be used depending on the value of the 'frequencies_are_increasing' 
	// flag.  If this flag is set, then we need to reorder by using a negative stride.
	
	bool incflag = frequencies_are_increasing;
	float *dst_int = incflag ? (intensity + (nfreq-1)*istride + it_chunk) : (intensity + it_chunk);
	float *dst_wt = incflag ? (weights + (nfreq-1)*wstride + it_chunk) : (weights + it_chunk);
	ssize_t dst_istride = incflag ? (-istride) : istride;
	ssize_t dst_wstride = incflag ? (-wstride) : wstride;
	
	// Call read_data() method in subclass.
	this->read_data(dst_int, dst_wt, it_file, n, dst_istride, dst_wstride);
	
	it_file += n;
	it_chunk += n;	    
    }
}


void chime_file_stream_base::_end_pipeline(Json::Value &json_output)
{
    if ((curr_ifile >= 0) && (curr_ifile < int(filename_list.size()))) {
	close_file();
	curr_ifile = -1;
    }
}


void chime_file_stream_base::_unbind_stream()
{
    if ((curr_ifile >= 0) && (curr_ifile < int(filename_list.size())))
	close_file();

    this->curr_ifile = -1;
    this->initial_discard_count = 0;
}
    

}   // namespace rf_pipelines
