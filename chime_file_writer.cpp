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

shared_ptr<wi_transform> make_chime_file_writer(const string &filename, bool clobber, int bitshuffle, ssize_t nt_chunk)
{
    throw runtime_error("rf_pipelines::make_chime_file_writer() was called, but rf_pipelines was compiled without ch_frb_io");
}

#else  // HAVE_CH_FRB_IO


struct chime_file_writer : public wi_transform {
    // Constructor args
    const string filename;
    const bool clobber;
    const int bitshuffle;

    // Stream params (not available until set_stream() gets called)
    double freq_lo_MHz = 0.0;
    double freq_hi_MHz = 0.0;
    double dt_sample = 0.0;

    std::unique_ptr<ch_frb_io::intensity_hdf5_ofile> ofile;
    int ichunk = 0;
    
    // Used to repack 'intensity' and 'weights' into contiguous arrays
    // FIXME these can go away when ch_frb_io::intensity_hdf5_ofile::append_chunk() supports a 'stride' argument
    uptr<float> intensity_contig_buf;
    uptr<float> weights_contig_buf;


    chime_file_writer(const string &filename_, bool clobber_, int bitshuffle_, ssize_t nt_chunk_) :
	wi_transform("chime_file_writer"),
	filename(filename_),
	clobber(clobber_), 
	bitshuffle(bitshuffle_)
    {
	this->name = "chime_file_writer(" + filename + ")";
	this->nt_chunk = nt_chunk_;
    }

    
    virtual void _bind_transform(Json::Value &json_attrs) override
    {
	if (!json_attrs.isMember("freq_lo_MHz") || !json_attrs.isMember("freq_hi_MHz"))
	    throw runtime_error("chime_file_writer: expected json_attrs to contain members 'freq_lo_MHz' and 'freq_hi_MHz'");
	
	if (!json_attrs.isMember("dt_sample"))
	    throw runtime_error("chime_file_writer: expected json_attrs to contain member 'dt_sample'");

	this->freq_lo_MHz = json_attrs["freq_lo_MHz"].asDouble();
	this->freq_hi_MHz = json_attrs["freq_hi_MHz"].asDouble();
	this->dt_sample = json_attrs["dt_sample"].asDouble();
    }


    virtual void _allocate() override
    {
	this->intensity_contig_buf = make_uptr<float> (nfreq * nt_chunk);
	this->weights_contig_buf = make_uptr<float> (nfreq * nt_chunk);
    }


    virtual void _start_pipeline(Json::Value &j) override
    {
	if (!clobber && file_exists(filename))
	    throw runtime_error(filename + ": file already exists and clobber=false was specified in the the chime_file_writer constructor");

	// Not really correct but that's OK
	vector<string> pol = { "XX" };

	// Note swapped ordering of freq_hi_MHz and freq_lo_MHz.  This is because rf_pipelines always orders frequency channels from highest to lowest.
	this->ofile = make_unique<ch_frb_io::intensity_hdf5_ofile> (filename, nfreq, pol, freq_hi_MHz, freq_lo_MHz, dt_sample, 0, 0, bitshuffle, nt_chunk);
	this->ichunk = 0;
    }


    virtual void _process_chunk(float *intensity, ssize_t istride, float *weights, ssize_t wstride, ssize_t pos) override
    {
	// Repack to contiguous arrays
	for (int ifreq = 0; ifreq < nfreq; ifreq++) {
	    memcpy(&intensity_contig_buf[0] + ifreq*nt_chunk, intensity + ifreq*istride, nt_chunk * sizeof(float));
	    memcpy(&weights_contig_buf[0] + ifreq*nt_chunk, weights + ifreq*wstride, nt_chunk * sizeof(float));
	}

	this->ofile->append_chunk(nt_chunk, &intensity_contig_buf[0], &weights_contig_buf[0], ichunk * nt_chunk, pos * dt_sample);
	this->ichunk++;
    }


    virtual void _end_pipeline(Json::Value &j) override
    {
	// Resetting this pointer will close file
	this->ofile.reset();
    }


    virtual void _reset() override
    {
	// Resetting this pointer will close file
	this->ofile.reset();
	this->ichunk = 0;
    }


    virtual void _deallocate() override
    {
	this->intensity_contig_buf.reset();
	this->weights_contig_buf.reset();
    }


    virtual Json::Value jsonize() const override
    {
	Json::Value ret;
	ret["class_name"] = "chime_file_writer";
	ret["filename"] = filename;
	ret["clobber"] = clobber;
	ret["bitshuffle"] = bitshuffle;
	ret["nt_chunk"] = int(nt_chunk);
	return ret;
    }

    static shared_ptr<chime_file_writer> from_json(const Json::Value &j)
    {
	string filename = string_from_json(j, "filename");
	bool clobber = bool_from_json(j, "clobber");
	int bitshuffle = int_from_json(j, "bitshuffle");
	int nt_chunk = int_from_json(j, "nt_chunk");
	return make_shared<chime_file_writer> (filename, clobber, bitshuffle, nt_chunk);
    }
};


namespace {
    struct _init {
	_init() {
	    pipeline_object::register_json_deserializer("chime_file_writer", chime_file_writer::from_json);
	}
    } init;
}


// See rf_pipelines.hpp for an explanation of the arguments
shared_ptr<wi_transform> make_chime_file_writer(const string &filename, bool clobber, int bitshuffle, ssize_t nt_chunk)
{
    return make_shared<chime_file_writer> (filename, clobber, bitshuffle, nt_chunk);
}

#endif  // HAVE_CH_FRB_IO


}  // namespace rf_pipelines
