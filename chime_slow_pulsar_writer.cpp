#include "rf_pipelines_internals.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif


chime_slow_pulsar_writer::chime_slow_pulsar_writer(ssize_t nt_chunk_) :
    wi_transform("chime_slow_pulsar_writer")
{
    // Note: nt_chunk is defined in wi_transform base class.
    this->nt_chunk = nt_chunk_;
}


void chime_slow_pulsar_writer::init_real_time_state(const real_time_state &rt_state_)
{
    if (rt_state_.beam_id < 0)
	throw runtime_error("rf_pipelines::chime_slow_pulsar:writer::init_real_time_state(): 'beam_id' is negative, or uninitialized");
    if (!rt_state_.memory_pool)
	throw runtime_error("rf_pipelines::chime_slow_pulsar:writer::init_real_time_state(): 'memory_pool' is an empty pointer, or uninitialized");
    if (!rt_state_.output_devices)
	throw runtime_error("rf_pipelines::chime_slow_pulsar:writer::init_real_time_state(): 'output_devices' is an empty pointer, or uninitialized");
    
    this->rt_state = rt_state_;    
}


void chime_slow_pulsar_writer::set_output_file_params(const output_file_params &of_params_)
{
    // FIXME: probably want some sanity-checking on 'of_params_' here.
    // FIXME: need to make this thread-safe!

    this->of_params = of_params_;
}


// virtual override
void chime_slow_pulsar_writer::_process_chunk(float *intensity, ssize_t istride, float *weights, ssize_t wstride, ssize_t pos)
{
    cout << "!!! chime_slow_pulsar_writer::_process_chunk() not implemented yet !!!\n" << flush;
    throw runtime_error("chime_slow_pulsar_writer::_process_chunk() not implemented yet");
}


// virtual override
void chime_slow_pulsar_writer::_end_pipeline(Json::Value &json_output)
{
    // Empty for now, but will probably want to add code to flush incomplete file to disk.
}


// virtual override
Json::Value chime_slow_pulsar_writer::jsonize() const
{
    Json::Value ret;
    
    ret["class_name"] = "chime_slow_pulsar_writer";
    ret["nt_chunk"] = int(this->nt_chunk);

    return ret;
}


// static member function
shared_ptr<chime_slow_pulsar_writer> chime_slow_pulsar_writer::from_json(const Json::Value &j)
{
    ssize_t nt_chunk = ssize_t_from_json(j, "nt_chunk");
    return make_shared<chime_slow_pulsar_writer> (nt_chunk);
}


// This boilerplate is necessary to "deserialize" a chime_slow_pulsar_writer object from a json file.
namespace {
    struct _init {
	_init() {
	    pipeline_object::register_json_deserializer("chime_slow_pulsar_writer", chime_slow_pulsar_writer::from_json);
	}
    } init;
}


}  // namespace rf_pipelines
