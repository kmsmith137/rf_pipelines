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

shared_ptr<wi_stream> make_chime_stream_from_acqdir(const string &filename, ssize_t nt_chunk)
{
    throw runtime_error("rf_pipelines::make_chime_stream_from_acqdir() was called, but rf_pipelines was compiled without ch_frb_io");
}

shared_ptr<wi_stream> make_chime_stream_from_filename(const string &filename, ssize_t nt_chunk)
{
    throw runtime_error("rf_pipelines::make_chime_stream_from_filename() was called, but rf_pipelines was compiled without ch_frb_io");
}

shared_ptr<wi_stream> make_chime_stream_from_filename_list(const vector<string> &filename_list, ssize_t nt_chunk)
{
    throw runtime_error("rf_pipelines::make_chime_stream_from_filename_list() was called, but rf_pipelines was compiled without ch_frb_io");
}

#else  // HAVE_CH_FRB_IO


class chime_file_stream : public wi_stream
{
protected:
    // Specified at construction.  Note that the 'nt_chunk' argument becomes 
    vector<string> filename_list;

    shared_ptr<ch_frb_io::intensity_hdf5_file> curr_file;
    int curr_ifile;  // index of current file in filename_list

public:
    chime_file_stream(const vector<string> &filename_list_, ssize_t nt_chunk);
    virtual ~chime_file_stream() { }

    virtual void stream_body(wi_run_state &run_state);
};

    
chime_file_stream::chime_file_stream(const vector<string> &filename_list_, ssize_t nt_chunk) :
    filename_list(filename_list_)
{
    rf_assert(filename_list.size() > 0);

    this->curr_file = make_shared<ch_frb_io::intensity_hdf5_file> (filename_list[0]);
    this->curr_ifile = 0;

    this->nfreq = curr_file->nfreq;
    this->freq_lo_MHz = curr_file->freq_lo_MHz;
    this->freq_hi_MHz = curr_file->freq_hi_MHz;
    this->dt_sample = curr_file->dt_sample;
    this->nt_maxwrite = (nt_chunk > 0) ? nt_chunk : 1024;
}

    
// virtual
void chime_file_stream::stream_body(wi_run_state &run_state)
{
    ssize_t it_file = 0;    // time index in current file.  Can be negative!  This means there is a time gap between files.
    ssize_t it_chunk = 0;
    int nfiles = filename_list.size();

    float *intensity;
    float *weights;
    ssize_t stride;
    bool zero_flag = true;

    run_state.start_substream(curr_file->time_lo);
    run_state.setup_write(nt_maxwrite, intensity, weights, stride, zero_flag, curr_file->time_lo);

    for (;;) {
	//
	// At top of loop, the state of the stream is represented by the following variables;
	//   curr_file   -> pointer to currently open file object
	//   curr_ifile  -> index of current file
	//   it_file     -> time index within current file
	//   it_chunk    -> time index within output chunk
	//

	if (curr_ifile >= nfiles) {
	    // End of stream
	    // FIXME should be able to write fewer than nt_maxwrite samples here
	    run_state.finalize_write(nt_maxwrite);
	    run_state.end_substream();
	    return;
	}
	else if (it_file >= curr_file->nt_logical) {
	    // End of file.
	    curr_ifile++;
	    it_file = 0;

	    if (curr_ifile >= nfiles) {
		curr_file = shared_ptr<ch_frb_io::intensity_hdf5_file> ();
		continue;
	    }
	    
	    // Open next file and do consistency tests.
	    double old_t1 = curr_file->time_hi;
	    curr_file = make_shared<ch_frb_io::intensity_hdf5_file> (filename_list[curr_ifile]);

	    if (curr_file->nfreq != this->nfreq)
		throw runtime_error("chime_file_stream: not every .h5 file has the same number of frequency channels?!");
	    if (fabs(curr_file->freq_lo_MHz - this->freq_lo_MHz) > 1.0e-4 * this->freq_lo_MHz)
		throw runtime_error("chime_file_stream: not every .h5 file has the same frequency range?!");
	    if (fabs(curr_file->freq_hi_MHz - this->freq_hi_MHz) > 1.0e-4 * this->freq_hi_MHz)
		throw runtime_error("chime_file_stream: not every .h5 file has the same time sampling rate?!");
	    if (fabs(curr_file->dt_sample - this->dt_sample) > 1.0e-4 * this->dt_sample)
		throw runtime_error("chime_file_stream: not every .h5 file has the same time sampling rate?!");

	    double new_t0 = curr_file->time_lo;
	    double gap = (new_t0 - old_t1) / dt_sample;  // time gap between files, as multiple of dt_sample

	    // FIXME another option here would be to start a new substream if a large gap is encountered.
	    if (gap < -0.01)
		throw runtime_error("chime_file_stream: acquisition files are not in time order?!");
	    if (gap > 10000.01)
		throw runtime_error("chime_file_stream: excessively long gap in acquisition");

	    ssize_t ngap = (ssize_t)(gap + 0.5);  // round to integer
	    rf_assert(ngap >= 0);

	    if (fabs(gap-ngap) > 0.01)
		throw runtime_error("chime_file_stream: time gap between files is not an integer multiple of dt_sample");
	    
	    it_file = -ngap;
	}
	else if (it_chunk >= nt_maxwrite) {
	    // Write chunk
	    run_state.finalize_write(nt_maxwrite);
	    run_state.setup_write(nt_maxwrite, intensity, weights, stride, zero_flag, curr_file->time_lo + it_file * dt_sample); 
	    it_chunk = 0;
	}
	else if (it_file < 0) {
	    // Skip gap between files.
	    ssize_t n = min(-it_file, nt_maxwrite-it_chunk);
	    it_file += n;
	    it_chunk += n;
	}
	else {
	    // Read data from file
	    ssize_t n = min(curr_file->nt_logical - it_file, nt_maxwrite - it_chunk);
	    
	    //
	    // A note on frequency channel ordering.  In rf_pipelines, frequencies must 
	    // be ordered from highest frequency to lowest.  In the ch_frb_io::intensity_hdf5_file,
	    // either ordering can be used depending on the value of the 'frequencies_are_increasing' 
	    // flag.  If this flag is set, then we need to reorder by using a negative stride.
	    //
	    bool incflag = curr_file->frequencies_are_increasing;
	    float *dst_int = incflag ? (intensity + (nfreq-1)*stride + it_chunk) : (intensity + it_chunk);
	    float *dst_wt = incflag ? (weights + (nfreq-1)*stride + it_chunk) : (weights + it_chunk);
	    ssize_t dst_stride = incflag ? (-stride) : stride;

	    curr_file->get_unpolarized_intensity(dst_int, dst_wt, it_file, n, dst_stride);
	    it_file += n;
	    it_chunk += n;	    
	}
    }
}


// Returns true if filename is of the form NNNNNNNN.h5, where N=[0,9]
static bool is_chime_file(const string &basename)
{
    if (basename.size() != 11)
	return false;
    if (!endswith(basename, ".h5"))
	return false;

    for (int i = 0; i < 8; i++)
	if (!isdigit(basename.at(i)))
	    return false;

    return true;
}


// Lists all files of the form ${dirname}/NNNNNNNN.h5, where N=[0,9]
static void list_chime_acqdir(vector<string> &chime_files, const string &dirname, bool allow_empty=false)
{
    vector<string> all_files;
    listdir(all_files, dirname);

    bool wflag = false;
    chime_files.resize(0);

    for (unsigned int i = 0; i < all_files.size(); i++) {
	string basename = all_files[i];

	if (is_chime_file(basename))
	    chime_files.push_back(dirname + "/" + basename);
	else if (!wflag && (endswith(basename,".h5") || endswith(basename,".hdf5"))) {
	    cerr << "warning: directory '" << dirname << "' contains hdf5 files which are not of the 'CHIME' form NNNNNNNN.h5\n";
	    wflag = true;
	}
    }

    if (!allow_empty && chime_files.size()==0)
	throw runtime_error("No CHIME files found in acquisition directory '" + dirname + "'");

    std::sort(chime_files.begin(), chime_files.end());    
}


shared_ptr<wi_stream> make_chime_stream_from_acqdir(const string &dirname, ssize_t nt_chunk)
{
    bool allow_empty = false;
    vector<string> filename_list;
    list_chime_acqdir(filename_list, dirname, allow_empty);
    cerr << dirname << ": " << filename_list.size() << " data files found in acqdir";

    return make_chime_stream_from_filename_list(filename_list, nt_chunk);
}

shared_ptr<wi_stream> make_chime_stream_from_filename(const string &filename, ssize_t nt_chunk)
{
    vector<string> filename_list;
    filename_list.push_back(filename);

    return make_chime_stream_from_filename_list(filename_list, nt_chunk);    
}

shared_ptr<wi_stream> make_chime_stream_from_filename_list(const vector<string> &filename_list, ssize_t nt_chunk)
{
    if (filename_list.size() == 0)
	throw runtime_error("empty filename_list in make_chime_stream()");

    return make_shared<chime_file_stream> (filename_list, nt_chunk);
}


#endif  // HAVE_CH_FRB_IO


}   // namespace rf_pipelines
