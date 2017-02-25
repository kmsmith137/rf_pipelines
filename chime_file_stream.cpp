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


#ifndef HAVE_CH_FRB_IO

shared_ptr<wi_stream> make_chime_stream_from_acqdir(const string &filename, ssize_t nt_chunk, ssize_t noise_source_align)
{
    throw runtime_error("rf_pipelines::make_chime_stream_from_acqdir() was called, but rf_pipelines was compiled without ch_frb_io");
}

shared_ptr<wi_stream> make_chime_stream_from_filename(const string &filename, ssize_t nt_chunk, ssize_t noise_source_align)
{
    throw runtime_error("rf_pipelines::make_chime_stream_from_filename() was called, but rf_pipelines was compiled without ch_frb_io");
}

shared_ptr<wi_stream> make_chime_stream_from_filename_list(const vector<string> &filename_list, ssize_t nt_chunk, ssize_t noise_source_align)
{
    throw runtime_error("rf_pipelines::make_chime_stream_from_filename_list() was called, but rf_pipelines was compiled without ch_frb_io");
}

#else  // HAVE_CH_FRB_IO


class chime_file_stream : public chime_file_stream_base
{
protected:
    shared_ptr<ch_frb_io::intensity_hdf5_file> curr_file;

public:
    chime_file_stream(const vector<string> &filename_list_, ssize_t nt_chunk, ssize_t noise_source_align);
    virtual ~chime_file_stream() { }

protected:
    virtual void load_file(const std::string& filename);
    virtual void close_file();
    virtual void set_params_from_file();
    virtual void check_file_consistency();
    virtual void read_data(float* dst_int, float* dst_wt, ssize_t it_file, ssize_t n, ssize_t dst_stride);

};

    
chime_file_stream::chime_file_stream(const vector<string> &filename_list_, ssize_t nt_chunk, ssize_t noise_source_align_) :
    chime_file_stream_base(filename_list_, nt_chunk, noise_source_align_)
{
}

//virtual
void chime_file_stream::load_file(const string &fn) {
    this->curr_file = make_shared<ch_frb_io::intensity_hdf5_file>(fn);
}

//virtual
void chime_file_stream::close_file() {
     this->curr_file.reset();
}

//virtual
void chime_file_stream::set_params_from_file() {
    this->nfreq = curr_file->nfreq;
    this->freq_lo_MHz = curr_file->freq_lo_MHz;
    this->freq_hi_MHz = curr_file->freq_hi_MHz;
    this->dt_sample = curr_file->dt_sample;
    this->time_lo = curr_file->time_lo;
    this->time_hi = curr_file->time_hi;
    this->nt = curr_file->nt_logical;
    this->frequencies_are_increasing = curr_file->frequencies_are_increasing;
}

//virtual
void chime_file_stream::check_file_consistency() {
    if (curr_file->nfreq != this->nfreq)
        throw runtime_error("chime_file_stream: not every .h5 file has the same number of frequency channels?!");
    if (fabs(curr_file->freq_lo_MHz - this->freq_lo_MHz) > 1.0e-4 * this->freq_lo_MHz)
        throw runtime_error("chime_file_stream: not every .h5 file has the same frequency range?!");
    if (fabs(curr_file->freq_hi_MHz - this->freq_hi_MHz) > 1.0e-4 * this->freq_hi_MHz)
        throw runtime_error("chime_file_stream: not every .h5 file has the same time sampling rate?!");
    if (fabs(curr_file->dt_sample - this->dt_sample) > 1.0e-4 * this->dt_sample)
        throw runtime_error("chime_file_stream: not every .h5 file has the same time sampling rate?!");
}

//virtual
void chime_file_stream::read_data(float* dst_int, float* dst_wt, ssize_t it_file, ssize_t n, ssize_t dst_stride) {
	    curr_file->get_unpolarized_intensity(dst_int, dst_wt, it_file, n, dst_stride);
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
    chime_files.resize(0);
    bool wflag = false;

    for (const string &basename: listdir(dirname)) {
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


shared_ptr<wi_stream> make_chime_stream_from_acqdir(const string &dirname, ssize_t nt_chunk, ssize_t noise_source_align)
{
    bool allow_empty = false;
    vector<string> filename_list;
    list_chime_acqdir(filename_list, dirname, allow_empty);
    cerr << dirname << ": " << filename_list.size() << " data files found in acqdir\n";

    return make_chime_stream_from_filename_list(filename_list, nt_chunk, noise_source_align);
}

shared_ptr<wi_stream> make_chime_stream_from_filename(const string &filename, ssize_t nt_chunk, ssize_t noise_source_align)
{
    vector<string> filename_list;
    filename_list.push_back(filename);

    return make_chime_stream_from_filename_list(filename_list, nt_chunk, noise_source_align);
}

shared_ptr<wi_stream> make_chime_stream_from_filename_list(const vector<string> &filename_list, ssize_t nt_chunk, ssize_t noise_source_align)
{
    if (filename_list.size() == 0)
	throw runtime_error("empty filename_list in make_chime_stream()");

    return make_shared<chime_file_stream> (filename_list, nt_chunk, noise_source_align);
}


#endif  // HAVE_CH_FRB_IO


}   // namespace rf_pipelines
