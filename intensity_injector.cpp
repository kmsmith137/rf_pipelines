#include <algorithm>

#include "rf_pipelines_internals.hpp"
#include "rf_pipelines_inventory.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif

typedef lock_guard<mutex> ulock;

injector::injector(int nt_chunk) :
    wi_transform("injector") {
    this->nt_chunk = nt_chunk;
}

void injector::_bind_transform(Json::Value &json_attrs)
{
    // Should be redundant with asserts elsewhere in rf_pipelines, but just being paranoid!
    rf_assert(this->nds == 1);
}

void injector::inject(shared_ptr<inject_data> data) {
    ulock u(mutex);
    to_inject.push_back(data);
}

uint64_t injector::get_last_fpgacount_seen() {
    uint64_t rtn;
    {
        ulock u(mutex);
        rtn = last_fpgacount_processed;
    }
    return rtn;
}


void injector::_start_pipeline(Json::Value &j)
{
    this->initial_fpga_count = uint64_t_from_json(j, "initial_fpga_count");
    this->fpga_counts_per_sample = int_from_json(j, "fpga_counts_per_sample");
    this->fpga_counts_initialized = true;
}

void injector::_process_chunk(float *intensity, ssize_t istride, float *weights, ssize_t wstride, ssize_t pos)
{
    // Reminder: previous asserts have already checked that
    //   this->nds == 1
    if (!fpga_counts_initialized)
	throw runtime_error("rf_pipelines::injector internal error: fpga count fields were not initialized as expected");
    
    // Recall that the 'pos' argument is the current pipeline position in units of time samples (not FPGA counts)

    vector<shared_ptr<inject_data> > to_inject_now;
    {
        ulock u(mutex);
        //cout << "Injector transform: " << to_inject.size() << " chunks of data" << endl;
        for (int idata=0; idata<to_inject.size(); idata++) {
            auto data = to_inject[idata];

            // Index in the current chunk of data of "fpga0" of this injected data entry
            int sample0 = (data->fpga0 - this->initial_fpga_count) / this->fpga_counts_per_sample - pos;
            if (sample0 + data->max_offset < 0) {
                // This inject_data request's time has passed.
                to_inject.erase(to_inject.begin() + idata);
                idata--;
            }
            if (sample0 + data->min_offset >= this->nt_chunk)
                // This inject_data request is in the future.
                continue;
            //cout << "Data to inject overlaps this chunk!!" << endl;
            to_inject_now.push_back(data);
        }
    }
    for (const auto &data : to_inject_now) {
        int sample0 = (data->fpga0 - this->initial_fpga_count) / this->fpga_counts_per_sample - pos;
        int nf = 0;
        int ntotal = 0;
        int nbefore = 0;
        int nafter = 0;
        // offset into the "data" array
        int data_offset = 0;
        for (int i=0; i<data->sample_offset.size(); i++) {
            int this_data_offset = data_offset;
            data_offset += data->ndata[i];
            // sample start,end for this frequency
            int f0 = sample0 + data->sample_offset[i];
            int f1 = f0 + data->ndata[i];
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

Json::Value injector::jsonize() const
{
    Json::Value ret;
    ret["class_name"] = "injector";
    ret["nt_chunk"] = int(nt_chunk);
    return ret;
}

shared_ptr<injector>
injector::from_json(const Json::Value &j)
{
    ssize_t nt_chunk = ssize_t_from_json(j, "nt_chunk");
    return make_shared<injector>(nt_chunk);
}

namespace {
    struct _init {
        _init() {
            pipeline_object::register_json_deserializer("injector", injector::from_json);
        }
    } init;
}

// Externally callable
shared_ptr<injector> make_injector(int nt_chunk) {
    return make_shared<injector>(nt_chunk);
}

}  // namespace rf_pipelines

