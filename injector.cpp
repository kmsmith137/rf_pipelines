//#include <ch_frb_io.hpp>

#include "rf_pipelines_internals.hpp"
#include "rf_pipelines_inventory.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif

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
    to_inject.push_back(data);
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
    
    // The 'pos' argument is the current pipeline position in units of time samples (not FPGA counts)
    uint64_t fpga_counts_start = pos * this->fpga_counts_per_sample + this->initial_fpga_count;
    uint64_t fpga_counts_end = (pos + nt_chunk) * this->fpga_counts_per_sample + this->initial_fpga_count;

    for (auto data : to_inject) {
        // check for no-overlap
        if ((data->fpga0 >= fpga_counts_end) ||
            (data->fpga_max < fpga_counts_start))
            continue;

        // OVERLAP!!  FIXME
        cout << "Data to inject overlaps this chunk!!" << endl;
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

