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

shared_ptr<wi_stream> make_chime_stream_from_acqdir(const string &filename, ssize_t nt_chunk, ssize_t noise_source_align, ssize_t nfiles)
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

    virtual Json::Value jsonize() const override;
    static shared_ptr<chime_file_stream> from_json(const Json::Value &j);

protected:
    virtual void load_file(const std::string& filename) override;
    virtual void close_file() override;
    virtual void set_params_from_file() override;
    virtual void check_file_consistency() const override;
    virtual void read_data(float* dst_int, float* dst_wt, ssize_t it_file, ssize_t n, ssize_t dst_istride, ssize_t dst_wstride) const override;
};

    
chime_file_stream::chime_file_stream(const vector<string> &filename_list_, ssize_t nt_chunk, ssize_t noise_source_align_) :
    chime_file_stream_base("chime_file_stream", filename_list_, nt_chunk, noise_source_align_)
{
}

// virtual override
void chime_file_stream::load_file(const string &fn) {
    this->curr_file = make_shared<ch_frb_io::intensity_hdf5_file>(fn);
}

// virtual override
void chime_file_stream::close_file() {
     this->curr_file.reset();
}

// virtual override
void chime_file_stream::set_params_from_file() {
    this->nfreq = curr_file->nfreq;
    this->freq_lo_MHz = curr_file->freq_lo_MHz;
    this->freq_hi_MHz = curr_file->freq_hi_MHz;
    this->dt_sample = curr_file->dt_sample;
    this->time_lo = curr_file->time_lo;
    this->time_hi = curr_file->time_hi;
    this->nt_file = curr_file->nt_logical;
    this->frequencies_are_increasing = curr_file->frequencies_are_increasing;
}

// virtual override
void chime_file_stream::check_file_consistency() const {
    if (curr_file->nfreq != this->nfreq)
        throw runtime_error("chime_file_stream: not every .h5 file has the same number of frequency channels?!");
    if (fabs(curr_file->freq_lo_MHz - this->freq_lo_MHz) > 1.0e-4 * this->freq_lo_MHz)
        throw runtime_error("chime_file_stream: not every .h5 file has the same frequency range?!");
    if (fabs(curr_file->freq_hi_MHz - this->freq_hi_MHz) > 1.0e-4 * this->freq_hi_MHz)
        throw runtime_error("chime_file_stream: not every .h5 file has the same time sampling rate?!");
    if (fabs(curr_file->dt_sample - this->dt_sample) > 1.0e-4 * this->dt_sample)
        throw runtime_error("chime_file_stream: not every .h5 file has the same time sampling rate?!");
}

// virtual override
void chime_file_stream::read_data(float* dst_int, float* dst_wt, ssize_t it_file, ssize_t n, ssize_t dst_istride, ssize_t dst_wstride) const
{
    curr_file->get_unpolarized_intensity(dst_int, dst_wt, it_file, n, dst_istride, dst_wstride);
}

// virtual override
Json::Value chime_file_stream::jsonize() const
{
    Json::Value ret;
    Json::Value &jf = ret["filename_list"];

    ret["class_name"] = "chime_file_stream";
    ret["noise_source_align"] = int(noise_source_align);
    ret["nt_chunk"] = int(this->get_orig_nt_chunk());

    for (const string &f: filename_list)
	jf.append(f);

    return ret;
}

// static
shared_ptr<chime_file_stream> chime_file_stream::from_json(const Json::Value &j)
{
    Json::Value a = array_from_json(j, "filename_list");
    vector<string> vs;

    for (const Json::Value &v: a) {
	if (!v.isString())
	    throw runtime_error("chime_file_stream::from_json(): filename was not a string as expected");
	vs.push_back(v.asString());
    }

    int nt_chunk = int_from_json(j, "nt_chunk");
    int noise_source_align = int_from_json(j, "noise_source_align");
    return make_shared<chime_file_stream> (vs, nt_chunk, noise_source_align);
}


namespace {
    struct _init {
	_init() {
	    pipeline_object::register_json_constructor("chime_file_stream", chime_file_stream::from_json);
	}
    } init;
}


// -------------------------------------------------------------------------------------------------


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


shared_ptr<wi_stream> make_chime_stream_from_acqdir(const string &dirname, ssize_t nt_chunk, ssize_t noise_source_align, ssize_t nfiles)
{
    if (nfiles < 0)
	throw runtime_error("rf_pipelines::chime_stream_from_acqdir: 'nfiles' argument was negative");

    bool allow_empty = false;
    vector<string> filename_list;
    list_chime_acqdir(filename_list, dirname, allow_empty);
    cerr << dirname << ": " << filename_list.size() << " data files found in acqdir\n";

    if ((nfiles > 0) && (nfiles < (ssize_t)filename_list.size()))
	filename_list.resize(nfiles);

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
