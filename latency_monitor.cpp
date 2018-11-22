#include "rf_pipelines_internals.hpp"
#include "rf_pipelines_inventory.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif

latency_monitor::latency_monitor(int nt_chunk_, string where_) :
    wi_transform("latency_monitor"),
    where(where_)
{	
    stringstream ss;
    ss << "latency_monitor(nt_chunk=" << nt_chunk_ << ", where=" << where << ")";

    this->name = ss.str();
    this->nt_chunk = nt_chunk_;
    this->nds = 0;  // allows us to run in a wi_sub_pipeline

    if (nt_chunk == 0)
        throw runtime_error("rf_pipelines::latency_monitor: nt_chunk must be specified");
}

latency_monitor::~latency_monitor() {}

void latency_monitor::_process_chunk(float *intensity, ssize_t istride, float *weights, ssize_t wstride, ssize_t pos) {}

shared_ptr<latency_monitor>
latency_monitor::from_json(const Json::Value &j)
{
    ssize_t nt_chunk = ssize_t_from_json(j, "nt_chunk");
    string where = string_from_json(j, "where");
    return make_shared<latency_monitor> (nt_chunk, where);
}

Json::Value latency_monitor::jsonize() const
{
    Json::Value ret;
    ret["class_name"] = "latency_monitor";
    ret["nt_chunk"] = int(this->get_prebind_nt_chunk());
    ret["where"] = where;
    return ret;
}

namespace {
    struct _init {
        _init() {
            pipeline_object::register_json_deserializer("latency_monitor", latency_monitor::from_json);
        }
    } init;
}


// Externally callable
shared_ptr<latency_monitor> make_latency_monitor(int nt_chunk, const string &where) {
    return make_shared<latency_monitor>(nt_chunk, where);
}

}  // namespace rf_pipelines
