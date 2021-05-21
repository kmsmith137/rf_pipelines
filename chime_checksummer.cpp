#include <algorithm>
#include "rf_pipelines_internals.hpp"

#ifdef HAVE_CH_FRB_IO
#include <ch_frb_io.hpp>
#endif

extern "C" {
#include "md5.h"
}

using namespace std;
using namespace ch_frb_io;

namespace rf_pipelines {
#if 0
};  // pacify emacs c-mode!
#endif


#ifndef HAVE_CH_FRB_IO

shared_ptr<wi_transform> make_chime_checksummer(const std::vector<int> &beams)
{
    throw runtime_error("rf_pipelines::make_chime_checksummer() was called, but rf_pipelines was compiled without ch_frb_io");
}

#else  // HAVE_CH_FRB_IO

struct chime_checksummer : public wi_transform {
    // Constructor args
    const std::vector<int> beams;

    // Is this stream's beam ID in the list of active beams?
    bool active = false;

    chime_checksummer(const std::vector<int> &beams_) :
	wi_transform("chime_checksummer"),
        beams(beams_)
    {
	this->name = "chime_checksummer";
	this->nt_chunk = constants::nt_per_assembled_chunk;
    }
    
    virtual void _start_pipeline(Json::Value &j) override
    {
        int beam = j["beam_id"].asInt();
        if (beams.size() == 0) {
            active = true;
        } else {
            for (int b : beams) {
                if (beam == b) {
                    active = true;
                    break;
                }
            }
        }
    }

    virtual void _process_chunk(float *intensity, ssize_t istride, float *weights, ssize_t wstride, ssize_t pos) override
    {
        if (!active)
            return;

        md5_state_t md5;
        md5_byte_t md5_digest[17];
        md5_init(&md5);
        for (int f=0; f<this->nfreq; f++)
            md5_append(&md5, reinterpret_cast<md5_byte_t*>(intensity + f*istride),
                       sizeof(float) * this->nt_chunk);
        md5_finish(&md5, md5_digest);
        cout << "Chunk " << std::dec << pos << ": MD5 intensity ";
        for (int i=0; i<16; i++) {
            cout << (i>0 ? ":" : "") << std::hex << ((md5_digest[i] >> 4) & 0xf) << (md5_digest[i] & 0xf);
        }

        md5_init(&md5);
        for (int f=0; f<this->nfreq; f++)
            md5_append(&md5, reinterpret_cast<md5_byte_t*>(weights + f*wstride),
                       sizeof(float) * this->nt_chunk);
        md5_finish(&md5, md5_digest);
        cout << ", weight ";
        for (int i=0; i<16; i++) {
            cout << (i>0 ? ":" : "") << std::hex << ((md5_digest[i] >> 4) & 0xf) << (md5_digest[i] & 0xf);
        }
        cout << endl;
    }

    virtual Json::Value jsonize() const override
    {
	Json::Value ret;
	ret["class_name"] = "chime_checksummer";

        Json::Value jbeams;
        for (int b : beams)
            jbeams.append(Json::Value(b));
        ret["beams"] = jbeams;
	return ret;
    }

    static shared_ptr<chime_checksummer> from_json(const Json::Value &j)
    {
        std::vector<int> beams;
        if (j.isMember("beams")) {
            Json::Value a = array_from_json(j, "beams");
            for (const Json::Value &v: a) {
                if (!v.isIntegral())
                    throw runtime_error("parsing chime_checksummer JSON: expected 'beams' to be an array of ints");
                int b = v.asInt();
                beams.push_back(b);
            }
        }
	return make_shared<chime_checksummer> (beams);
    }
};


namespace {
    struct _init {
	_init() {
	    pipeline_object::register_json_deserializer("chime_checksummer", chime_checksummer::from_json);
	}
    } init;
}


// See rf_pipelines.hpp for an explanation of the arguments
shared_ptr<wi_transform> make_chime_checksummer(const std::vector<int> &beams)
{
    return make_shared<chime_checksummer> (beams);
}

#endif  // HAVE_CH_FRB_IO


}  // namespace rf_pipelines
