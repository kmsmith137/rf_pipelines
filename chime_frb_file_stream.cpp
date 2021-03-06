#include <glob.h>
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

shared_ptr<wi_stream> make_chime_frb_stream_from_glob(const string &glob_pattern, ssize_t nt_chunk, ssize_t noise_source_align)
{
    throw runtime_error("rf_pipelines::make_chime_frb_stream_from_glob() was called, but rf_pipelines was compiled without ch_frb_io");
}

shared_ptr<wi_stream> make_chime_frb_stream_from_filename(const string &filename, ssize_t nt_chunk, ssize_t noise_source_align)
{
    throw runtime_error("rf_pipelines::make_chime_frb_stream_from_filename() was called, but rf_pipelines was compiled without ch_frb_io");
}

shared_ptr<wi_stream> make_chime_frb_stream_from_filename_list(const vector<string> &filename_list, ssize_t nt_chunk, ssize_t noise_source_align)
{
    throw runtime_error("rf_pipelines::make_chime_frb_stream_from_filename_list() was called, but rf_pipelines was compiled without ch_frb_io");
}

#else  // HAVE_CH_FRB_IO

using namespace ch_frb_io;

class chime_frb_file_stream : public chime_file_stream_base
{
protected:
    shared_ptr<ch_frb_io::assembled_chunk> chunk;

public:
    chime_frb_file_stream(const vector<string> &filename_list_, ssize_t nt_chunk, ssize_t noise_source_align);
    virtual ~chime_frb_file_stream() { }

    virtual Json::Value jsonize() const override;
    static shared_ptr<chime_frb_file_stream> from_json(const Json::Value &j);

protected:
    virtual void load_file(const std::string &filename) override;
    virtual void close_file() override;
    virtual void set_params_from_file() override;
    virtual void check_file_consistency() const override;
    virtual void read_data(float* dst_int, float* dst_wt, ssize_t it_file, ssize_t n, ssize_t dst_istride, ssize_t dst_wstride) const override;
};

    
chime_frb_file_stream::chime_frb_file_stream(const vector<string> &filename_list_, ssize_t nt_chunk, ssize_t noise_source_align_) :
    chime_file_stream_base("chime_frb_file_stream", filename_list_, nt_chunk, noise_source_align_)
{
}

// virtual override
void chime_frb_file_stream::load_file(const string &fn) {
    this->chunk = assembled_chunk::read_msgpack_file(fn);
}

// virtual override
void chime_frb_file_stream::close_file() {
     this->chunk.reset();
}

// virtual override
void chime_frb_file_stream::set_params_from_file() {
    this->nfreq = chunk->nupfreq * ch_frb_io::constants::nfreq_coarse_tot;
    this->freq_lo_MHz = 400.0;
    this->freq_hi_MHz = 800.0;
    double fpga_t = chime_file_stream_base::chime_seconds_per_fpga_count;
    this->dt_sample = fpga_t * chunk->fpga_counts_per_sample;
    this->time_lo = fpga_t * chunk->fpga_begin;
    this->time_hi = fpga_t * chunk->fpga_end;
    this->nt_file = ch_frb_io::constants::nt_per_assembled_chunk;
    this->frequencies_are_increasing = false;
}

// virtual override
void chime_frb_file_stream::check_file_consistency() const {
    if (this->nfreq != chunk->nupfreq * ch_frb_io::constants::nfreq_coarse_tot)
        throw runtime_error("chime_frb_file_stream: nfreq mismatch!");
    double fpga_t = chime_file_stream_base::chime_seconds_per_fpga_count;
    if (this->dt_sample != fpga_t * chunk->fpga_counts_per_sample)
        throw runtime_error("chime_frb_file_stream: dt_sample mismatch!");
    if (this->nt_file != ch_frb_io::constants::nt_per_assembled_chunk)
        throw runtime_error("chime_frb_file_stream: nt_file mismatch!");
}

// virtual override
void chime_frb_file_stream::read_data(float* dst_int, float* dst_wt, ssize_t it_file, ssize_t n, ssize_t dst_istride, ssize_t dst_wstride) const {
    chunk->decode_subset(dst_int, dst_wt, it_file, n, dst_istride, dst_wstride);
}


// virtual override
Json::Value chime_frb_file_stream::jsonize() const
{
    Json::Value ret;
    Json::Value &jf = ret["filename_list"];

    ret["class_name"] = "chime_frb_file_stream";
    ret["nt_chunk"] = int(this->get_prebind_nt_chunk());
    ret["noise_source_align"] = int(noise_source_align);

    for (const string &f: filename_list)
	jf.append(f);

    return ret;
}

// static 
shared_ptr<chime_frb_file_stream> chime_frb_file_stream::from_json(const Json::Value &j)
{
    Json::Value a = array_from_json(j, "filename_list");
    vector<string> vs;

    for (const Json::Value &v: a) {
	if (!v.isString())
	    throw runtime_error("chime_frb_file_stream::from_json(): filename was not a string as expected");
	vs.push_back(v.asString());
    }

    int nt_chunk = int_from_json(j, "nt_chunk");
    int noise_source_align = int_from_json(j, "noise_source_align");
    return make_shared<chime_frb_file_stream> (vs, nt_chunk, noise_source_align);
}


namespace {
    struct _init {
	_init() {
	    pipeline_object::register_json_deserializer("chime_frb_file_stream", chime_frb_file_stream::from_json);
	}
    } init;
}


// -------------------------------------------------------------------------------------------------

// Lists all files matching given glob pattern.
static void list_glob(vector<string> &filenames, const string &glob_pattern, bool allow_empty=false) {
    filenames.resize(0);

    glob_t theglob;
    memset(&theglob, 0, sizeof(glob_t));

    if (glob(glob_pattern.c_str(), GLOB_BRACE | GLOB_TILDE, NULL, &theglob) == -1) {
        throw runtime_error("glob() failed: " + string(strerror(errno)));
    }

    for (size_t i=0; i<theglob.gl_pathc; i++) {
        filenames.push_back(string(theglob.gl_pathv[i]));
    }

    globfree(&theglob);
    if (!allow_empty && filenames.size()==0)
	throw runtime_error("No CHIME files found for glob " + glob_pattern);
}


shared_ptr<wi_stream> make_chime_frb_stream_from_glob(const string &glob_pattern, ssize_t nt_chunk, ssize_t noise_source_align)
{
    bool allow_empty = false;
    vector<string> filename_list;
    list_glob(filename_list, glob_pattern, allow_empty);
    cout << glob_pattern << ": " << filename_list.size() << " data files found\n";
    return make_chime_frb_stream_from_filename_list(filename_list, nt_chunk, noise_source_align);
}

shared_ptr<wi_stream> make_chime_frb_stream_from_filename(const string &filename, ssize_t nt_chunk, ssize_t noise_source_align)
{
    vector<string> filename_list;
    filename_list.push_back(filename);
    return make_chime_frb_stream_from_filename_list(filename_list, nt_chunk, noise_source_align);
}

shared_ptr<wi_stream> make_chime_frb_stream_from_filename_list(const vector<string> &filename_list, ssize_t nt_chunk, ssize_t noise_source_align)
{
    if (filename_list.size() == 0)
	throw runtime_error("empty filename_list in make_chime_frb_stream()");
    return make_shared<chime_frb_file_stream> (filename_list, nt_chunk, noise_source_align);
}


#endif  // HAVE_CH_FRB_IO


}   // namespace rf_pipelines
