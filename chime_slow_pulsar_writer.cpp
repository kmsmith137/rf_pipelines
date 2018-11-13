#include "rf_pipelines_internals.hpp"
#include <rf_kernels/downsample.hpp>
#include <fstream>

using namespace std;
using namespace rf_kernels;

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif


chime_slow_pulsar_writer::chime_slow_pulsar_writer(ssize_t nt_chunk_) :
    wi_transform("chime_slow_pulsar_writer")
{
    // Note: nt_chunk is defined in wi_transform base class.
    this->nt_chunk = nt_chunk_;
    // this->of_params = new of_params
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

    std::lock_guard<std::mutex> lock(this->of_mutex);
    this->of_params = of_params_;
    this->of_params.nt_out = this->nt_chunk / this->of_params.nds_out / this->nds;
    this->tmp_i = std::shared_ptr<std::vector<float>>(new std::vector<float>(this->of_params.nfreq_out * this->of_params.nt_out));
    this->tmp_w = std::shared_ptr<std::vector<float>>(new std::vector<float>(this->of_params.nfreq_out * this->of_params.nt_out));

    // const int nfreq_out = this->of_params.nfreq_out;
    // const int nt_out = this->nt_chunk / this->of_params.nds_out;
    this->downsampler = std::shared_ptr<rf_kernels::wi_downsampler>(new rf_kernels::wi_downsampler(
                                                             this->nfreq/this->of_params.nfreq_out, this->of_params.nds_out));
}

template<typename T>
T* get_ptr(std::shared_ptr<std::vector<T>> vect)
{
    return &((*vect)[0]);
};

// virtual override
void chime_slow_pulsar_writer::_process_chunk(float *intensity, ssize_t istride, float *weights, ssize_t wstride, ssize_t pos)
{
    // we just need to extract the parameters locally
    this->of_mutex.lock();
    int nfreq_out = this->of_params.nfreq_out;
    int nt_out = this->of_params.nt_out;
    int nds = this->nds;
    auto tmp_i = this->tmp_i;
    auto tmp_w = this->tmp_w;

    this->of_mutex.unlock();

    // std::cout << "current downsample setting: (nfreq_out) " << this->of_params.nfreq_out << " (nt_out) " << this->of_params.nt_out << std::endl;
    if( (nfreq_out > 0) && (nt_out > 0) ){
        // do downsample   
        // int nfreq_out, int nt_out,
        //         float *out_i, int out_istride,
        //         float *out_w, int out_wstride,
        //         const float *in_i, int in_istride,
        //         const float *in_w, int in_wstride);
        std::cout << "downsampling" << std::endl;

        this->downsampler->downsample(nfreq_out, nt_out, 
                                        get_ptr<float>(tmp_i), nt_out,
                                        get_ptr<float>(tmp_w), nt_out,
                                        intensity, istride,
                                        weights, wstride);
        
        std::cout << "writing" << std::endl;
        // basic file io for testing
        std::string fname = "test_out.spdat";
        {
            std::ofstream of(fname, std::ios::binary);
            of.seekp(std::ios::end);
            of.write(reinterpret_cast<char*>(&(this->of_params)), 
                sizeof(chime_slow_pulsar_writer::output_file_params));
            of.write(reinterpret_cast<char*>(&((*this->tmp_i)[0])),
                sizeof(float) * nfreq_out * nt_out);
            of.write(reinterpret_cast<char*>(&((*this->tmp_w)[0])),
                sizeof(float) * nfreq_out * nt_out);
            of.close()
        }
    }
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
