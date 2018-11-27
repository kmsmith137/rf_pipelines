#include <algorithm>

#include "rf_pipelines_internals.hpp"
#include "rf_pipelines_inventory.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif

typedef lock_guard<mutex> ulock;

intensity_injector::intensity_injector(int nt_chunk) :
    wi_transform("intensity_injector") {
    this->nt_chunk = nt_chunk;
}

void intensity_injector::_bind_transform(Json::Value &json_attrs)
{
    // Should be redundant with asserts elsewhere in rf_pipelines, but just being paranoid!
    rf_assert(this->nds == 1);
}

void intensity_injector::_start_pipeline(Json::Value &json_attrs) {}

void intensity_injector::inject(shared_ptr<inject_data> data) {
    ulock u(mutex);
    to_inject.push_back(data);
}

void intensity_injector::_process_chunk(float *intensity, ssize_t istride, float *weights, ssize_t wstride, ssize_t pos)
{
    // Reminder: previous asserts have already checked that
    //   this->nds == 1
    
    // Recall that the 'pos' argument is the current pipeline position
    // in units of time samples (starts from 0 for first chunk
    // processed by this rf_pipeline)

    vector<shared_ptr<inject_data> > to_inject_now;
    {
        ulock u(mutex);
        //cout << "Intensity_Injector transform: " << to_inject.size() << " chunks of data" << endl;
        for (int idata=0; idata<to_inject.size(); idata++) {
            auto data = to_inject[idata];
            // Index in the current chunk of data of "fpga0" of this injected data entry
            if ((data->sample0 + data->max_offset) < pos) {
                // This inject_data request's time has passed.
                to_inject.erase(to_inject.begin() + idata);
                idata--;
            }
            if ((data->sample0 + data->min_offset) >= (pos + this->nt_chunk)) {
                // This inject_data request is in the future.
                continue;
            }
            to_inject_now.push_back(data);
        }
    }
    for (const auto &data : to_inject_now) {
        // About int sizes here: data->sample0 may be large, as may
        // pos (if we run for a long time).

        // Index in this chunk of offset zero in the inject_data request;
        // may be positive or negative.
        ssize_t sample0 = data->sample0 - pos;
        int nf = 0;
        int ntotal = 0;
        int nbefore = 0;
        int nafter = 0;
        // offset into the "data" array
        ssize_t data_offset = 0;
        for (int i=0; i<data->sample_offset.size(); i++) {
            // save this index in case we need it below
            ssize_t this_data_offset = data_offset;
            // skip this frequency's data regardless of whether we use them
            data_offset += data->ndata[i];
            // sample start,end for this frequency
            ssize_t f0 = sample0 + data->sample_offset[i];
            ssize_t f1 = f0 + data->ndata[i];
            if (f0 >= nt_chunk) {
                // This frequency's data is after this chunk
                nafter++;
                continue;
            }
            if (f1 <= 0) {
                // This frequency's data is before this chunk
                nbefore++;
                continue;
            }
            // Now select the subset of this frequency's data that overlaps
            // this chunk.  f0 and f1 are indices into this chunk of data.
            int inj_offset = 0;
            int ncopy = data->ndata[i];
            if (f0 < 0) {
                // we're injecting the tail end of this frequency's array
                ncopy -= (-f0);
                inj_offset += (-f0);
                // it will start at the beginning of this chunk
                f0 = 0;
            }
            if (f0 + ncopy > nt_chunk) {
                ncopy = nt_chunk - f0;
                // a tail of data remains
            }
            float* indata = intensity + i*istride + f0;
            for (int j=0; j<ncopy; j++)
                indata[j] += data->data[this_data_offset + inj_offset + j];
            nf += 1;
            ntotal += ncopy;
        }
        //cout << "Injected " << nf << " frequency bins, total of " << ntotal << " samples.  N freq before: " << nbefore << ", after: " << nafter << endl;
    }
}

Json::Value intensity_injector::jsonize() const
{
    Json::Value ret;
    ret["class_name"] = "intensity_injector";
    ret["nt_chunk"] = int(nt_chunk);
    return ret;
}

shared_ptr<intensity_injector>
intensity_injector::from_json(const Json::Value &j)
{
    ssize_t nt_chunk = ssize_t_from_json(j, "nt_chunk");
    return make_shared<intensity_injector>(nt_chunk);
}

namespace {
    struct _init {
        _init() {
            pipeline_object::register_json_deserializer("intensity_injector", intensity_injector::from_json);
        }
    } init;
}

// Externally callable
shared_ptr<intensity_injector> make_intensity_injector(int nt_chunk) {
    return make_shared<intensity_injector>(nt_chunk);
}

}  // namespace rf_pipelines

