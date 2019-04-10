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


void chime_slow_pulsar_writer::get_new_chunk()
{
    // NOTE: complete the initializer!
    std::lock_guard<std::mutex> lg(this->chunk_mutex);
    std::shared_ptr<ch_chunk_initializer> ini_params(new ch_chunk_initializer(this->rt_state.memory_pool));
    ini_params->fpga_counts_per_sample = this->fpga_counts_per_sample;
    ini_params->frame0_nano = this->frame0_nano;
    this->chunk = slow_pulsar_chunk::make_slow_pulsar_chunk(ini_params);
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
    this->get_new_chunk();
}


void chime_slow_pulsar_writer::set_output_file_params(const output_file_params &of_params_)
{
    // FIXME: probably want some sanity-checking on 'of_params_' here.

    std::lock_guard<std::mutex> lg(this->of_mutex);
    this->of_params = of_params_;
    this->of_params.nt_out = this->nt_chunk / this->of_params.nds_out / this->nds;
    // I consider this an acceptable use of "new"
    this->tmp_i = std::shared_ptr<std::vector<float>>(new std::vector<float>(this->of_params.nfreq_out * this->of_params.nt_out));
    this->tmp_w = std::shared_ptr<std::vector<float>>(new std::vector<float>(this->of_params.nfreq_out * this->of_params.nt_out));
    const ssize_t nbytes_charbuf = ch_frb_io::byte_ceil(2 * this->of_params.nbits_out * this->of_params.nfreq_out * this->of_params.nt_out);
    this->nbytes_charbuf = nbytes_charbuf;
    this->tmp_buf = std::shared_ptr<std::vector<char>>(new std::vector<char>(nbytes_charbuf));

    this->downsampler = std::shared_ptr<rf_kernels::wi_downsampler>(new rf_kernels::wi_downsampler(
                                                             this->nfreq/this->of_params.nfreq_out, this->of_params.nds_out));
}


template<typename T>
T* get_ptr(std::shared_ptr<std::vector<T>> vect)
{
    return &((*vect)[0]);
};


// this is a hack to get around sp_header forward declaration in the hpp file
void quantize_store(chime_slow_pulsar_writer* spw, fvec_t in, const ssize_t istride, 
                    fvec_t weights, const ssize_t wstride, const int nbits, sp_header& sph)
{
    // quantize logic here

    // attempt to commit the chunk to slab
    const bool success = spw->chunk->commit_chunk(sph, spw->tmp_buf);

    if(!success){
        // std::cout<< "writing output file" << std::endl;
        // enqueue write request
        std::shared_ptr<write_chunk_request> req(new write_chunk_request());
        req->filename = "test/test.dat";
        req->chunk = spw->chunk;
        spw->rt_state.output_devices->enqueue_write_request(req);
        spw->get_new_chunk();
        quantize_store(spw, in, istride, weights, wstride, nbits, sph);
    }
}


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

    sp_header sph;
    sph.ichunk = this->ichunk;
    sph.nt = nt_out;
    sph.ntds = nds;
    sph.nfreq = nfreq_out;
    sph.nbits = nbits_out;

    // std::cout << "current downsample setting: (nfreq_out) " << this->of_params.nfreq_out << " (nt_out) " << this->of_params.nt_out << std::endl;
    if( (nfreq_out > 0) && (nt_out > 0) ){

        this->downsampler->downsample(nfreq_out, nt_out, 
                                        get_ptr<float>(tmp_i), nt_out,
                                        get_ptr<float>(tmp_w), nt_out,
                                        intensity, istride,
                                        weights, wstride);
        
        //this is inelegant
        quantize_store(this, tmp_i, nt_out, tmp_w, nt_out, nbits_out, sph);
    }

    this->ichunk += 1;
}


// void chime_slow_pulsar_writer::quantize_store(fvec_t in, const ssize_t istride, fvec_t weights, const ssize_t wstride, const int nbits,
//                                                 sp_header& sph)
// {

//     //quantize_kernel

//     const bool success = this->chunk->commit_to_chunk(sph, this->tmp_buf);

//     if(!success){
//         // enqueue write request
//         std::shared_ptr<ch_frb_io::write_chunk_request>> wreq();
//         wreq->filename = "test/test.dat";
//         wreq->chunk = this->chunk;
//         this->output_devices->enqueue_write_request(wreq);
//         this->get_new_chunk();
//         this->quantize_store(in, istride, weights, wstride, nbits, sph);
//     }
// }


// virtual override
void chime_slow_pulsar_writer::_start_pipeline(Json::Value &json_attrs)
{
    if (!json_attrs.isMember("frame0_nano") || !json_attrs.isMember("fpga_counts_per_sample"))
        throw runtime_error("chime_slow_pulsar_writer: expected json_attrs to contain members 'frame0_nano' and 'fpga_counts_per_sample'");
    
    // acquire the outfile lock
    std::lock_guard<std::mutex> lg(this->of_mutex);
    std::lock_guard<std::mutex> lg2(this->chunk_mutex);
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
