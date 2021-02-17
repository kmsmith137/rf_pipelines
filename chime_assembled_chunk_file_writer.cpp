#include <algorithm>
#include "rf_pipelines_internals.hpp"

#ifdef HAVE_CH_FRB_IO
#include <ch_frb_io.hpp>
#endif

using namespace std;
using namespace ch_frb_io;

namespace rf_pipelines {
#if 0
};  // pacify emacs c-mode!
#endif


#ifndef HAVE_CH_FRB_IO

shared_ptr<wi_transform> make_chime_assembled_chunk_file_writer(const string &filename, bool clobber)
{
    throw runtime_error("rf_pipelines::make_chime_assembled_chunk_file_writer() was called, but rf_pipelines was compiled without ch_frb_io");
}

#else  // HAVE_CH_FRB_IO

struct chime_assembled_chunk_file_writer : public wi_transform {
    // Constructor args
    const string filename;
    const bool clobber;

    // Stream params (not available until set_stream() gets called)
    assembled_chunk::initializer chunk_ini;
    int ichunk_offset = 0;

    chime_assembled_chunk_file_writer(const string &filename_, bool clobber_) :
	wi_transform("chime_assembled_chunk_file_writer"),
	filename(filename_),
	clobber(clobber_)
    {
	this->name = "chime_assembled_chunk_file_writer(" + filename + ")";
	this->nt_chunk = constants::nt_per_assembled_chunk;
    }
    
    virtual void _start_pipeline(Json::Value &j) override
    {
        // FIXME --
        chunk_ini.beam_id = 0;
        chunk_ini.nupfreq = this->nfreq / constants::nfreq_coarse_tot;
        chunk_ini.nrfifreq = 0; // ??
        chunk_ini.nt_per_packet = 16; // ??
        chunk_ini.fpga_counts_per_sample = j["fpga_counts_per_sample"].asInt();
        chunk_ini.frame0_nano = j["frame0_nano"].asUInt64();
        // This converts ch_frb_io chunk number (which scales directly
        // to FPGA counts) to rf_pipelines chunk counts (which start
        // from 0 at the beginning of the stream).
        this->ichunk_offset = j["initial_fpga_count"].asUInt64() / (uint64_t)(chunk_ini.fpga_counts_per_sample * constants::nt_per_assembled_chunk);
    }

    virtual void _process_chunk(float *intensity, ssize_t istride, float *weights, ssize_t wstride, ssize_t pos) override
    {
        chunk_ini.ichunk = this->ichunk_offset + pos / constants::nt_per_assembled_chunk;
        shared_ptr<assembled_chunk> ch = assembled_chunk::make(chunk_ini);
        float ilo = 1e9;
        float ihi = -1e9;
        int nbad = 0;
        for (int f=0; f<this->nfreq; f++) {
            for (int t=0; t<this->nt_chunk; t++) {
                if (weights[f*istride + t] == 0) {
                    nbad++;
                    continue;
                }
                float ival = intensity[f*istride + t];
                if (ival < ilo)
                    ilo = ival;
                if (ival > ihi)
                    ihi = ival;
            }
        }
        // Set scale and offset to fill 2 to 252.
        // intensity = scale * (8-bit value) + offset
        float scale = (ihi - ilo) / 250.;
        float offset = ilo - 2.*scale;
        for (int i=0; i<ch->nscales; i++) {
            ch->scales[i] = scale;
            ch->offsets[i] = offset;
        }
        int i = 0;
        for (int f=0; f<this->nfreq; f++) {
            for (int t=0; t<this->nt_chunk; t++) {
                uint8_t val;
                if (weights[f*istride + t] == 0)
                    val = 0;
                else
                    val = uint8_t((intensity[f*istride + t] - offset) / scale);
                ch->data[i] = val;
                i++;
            }
        }
        string thisfn = ch->format_filename(filename);
	if (!clobber && file_exists(thisfn))
	    throw runtime_error(thisfn + ": file already exists and clobber=false was specified in the the chime_assembled_chunk_file_writer constructor");
        ch->write_msgpack_file(thisfn, false);
    }

    virtual Json::Value jsonize() const override
    {
	Json::Value ret;
	ret["class_name"] = "chime_assembled_chunk_file_writer";
	ret["filename"] = filename;
	ret["clobber"] = clobber;
	return ret;
    }

    static shared_ptr<chime_assembled_chunk_file_writer> from_json(const Json::Value &j)
    {
	string filename = string_from_json(j, "filename");
	bool clobber = bool_from_json(j, "clobber");
	return make_shared<chime_assembled_chunk_file_writer> (filename, clobber);
    }
};


namespace {
    struct _init {
	_init() {
	    pipeline_object::register_json_deserializer("chime_assembled_chunk_file_writer", chime_assembled_chunk_file_writer::from_json);
	}
    } init;
}


// See rf_pipelines.hpp for an explanation of the arguments
shared_ptr<wi_transform> make_chime_assembled_chunk_file_writer(const string &filename, bool clobber)
{
    return make_shared<chime_assembled_chunk_file_writer> (filename, clobber);
}

#endif  // HAVE_CH_FRB_IO


}  // namespace rf_pipelines
