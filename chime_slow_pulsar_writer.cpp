#include "rf_pipelines_internals.hpp"
#include <rf_kernels/downsample.hpp>
#include <fstream>

#include <ch_frb_io.hpp>

using namespace std;
using namespace rf_kernels;
using namespace ch_frb_io;

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif

struct chime_slow_pulsar_context
{

};

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
    ch_chunk::initializer ini_params;
    ini_params.pool = this->rt_state.memory_pool;
    chunk = slow_pulsar_chunk::make_slow_pulsar_chunk(ini_params);
}


void chime_slow_pulsar_writer::set_output_file_params(const output_file_params &of_params_)
{
    // FIXME: probably want some sanity-checking on 'of_params_' here.

    std::lock_guard<std::mutex> lock(this->of_mutex);
    this->of_params = of_params_;
    this->of_params.nt_out = this->nt_chunk / this->of_params.nds_out / this->nds;
    // I consider this an acceptable use of "new"
    this->tmp_i = std::shared_ptr<std::vector<float>>(new std::vector<float>(this->of_params.nfreq_out * this->of_params.nt_out));
    this->tmp_w = std::shared_ptr<std::vector<float>>(new std::vector<float>(this->of_params.nfreq_out * this->of_params.nt_out));

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
    int nbits_out = this->of_params.nbits_out;
    auto tmp_i = this->tmp_i;
    auto tmp_w = this->tmp_w;
    uint64_t frame0 = this->frame0_nano;
    uint64_t fpga_counts_per_sample = this->fpga_counts_per_sample;
    uint64_t nsamp_chunk = this->nt_chunk;
    this->of_mutex.unlock();

    // std::cout << "current downsample setting: (nfreq_out) " << this->of_params.nfreq_out << " (nt_out) " << this->of_params.nt_out << std::endl;
    if( (nfreq_out > 0) && (nt_out > 0) ){

        this->downsampler->downsample(nfreq_out, nt_out, 
                                        get_ptr<float>(tmp_i), nt_out,
                                        get_ptr<float>(tmp_w), nt_out,
                                        intensity, istride,
                                        weights, wstride);
        this->quantize_store(tmp_i, nt_out, tmp_w, nt_out, nbits_out);
        // basic file io for testing
        // std::string fname = "test_out.spdat";
        // {
        //     std::ofstream of(fname, std::ios::binary | std::ios::app);
        //     of.seekp(std::ios::end);
        //     of.write(reinterpret_cast<char*>(&(this->of_params)), 
        //         sizeof(chime_slow_pulsar_writer::output_file_params));
        //     of.write(reinterpret_cast<char*>(&((*this->tmp_i)[0])),
        //         sizeof(float) * nfreq_out * nt_out);
        //     of.write(reinterpret_cast<char*>(&((*this->tmp_w)[0])),
        //         sizeof(float) * nfreq_out * nt_out);
        // }
    }

    this->nchunk += 1;
}

void chime_slow_pulsar_writer::quantize_store(fvec_t in, const ssize_t istride, fvec_t weights, const ssize_t wstride, const int nbits)
{
    //no-op
}

// // make sure there's a slab allocated
// bool chime_slow_pulsar_writer::verify_slab()
// {
//     // must set null
//     // need threadsafe?
//     if(!working_slab){
//         working_slab = memory_slab_t(this->rt_state.memory_pool->get_slab(false, true));
//     }
// }

// // virtual override
// void chime_slow_pulsar_writer::_bind_transform(Json::Value &json_attrs)
// {
//     if (!json_attrs.isMember("frame0_nano") || !json_attrs.isMember("fpga_counts_per_sample"))
//         throw runtime_error("chime_slow_pulsar_writer: expected json_attrs to contain members 'frame0_nano' and 'fpga_counts_per_sample'");
    
//     this->frame0_nano = json_attrs['fpga_counts_per_sample'].asUInt64();
//     this->fpga_counts_per_sample = json_attrs['fpga_counts_per_sample'].asUInt64();

//     std::cout << "frame0: " << this->frame0_nano << " fpga_counts_per_sample: " << this->fpga_counts_per_sample << std::endl;
// }

// virtual override
void chime_slow_pulsar_writer::_start_pipeline(Json::Value &json_attrs)
{
    if (!json_attrs.isMember("frame0_nano") || !json_attrs.isMember("fpga_counts_per_sample"))
        throw runtime_error("chime_slow_pulsar_writer: expected json_attrs to contain members 'frame0_nano' and 'fpga_counts_per_sample'");
    
    // acquire the outfile lock
    std::lock_guard<std::mutex> lg(this->of_mutex);
    this->frame0_nano = json_attrs["frame0_nano"].asUInt64();
    this->fpga_counts_per_sample = json_attrs["fpga_counts_per_sample"].asUInt64();

    // std::cout << "frame0: " << this->frame0_nano << " fpga_counts_per_sample: " << this->fpga_counts_per_sample << std::endl;
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
