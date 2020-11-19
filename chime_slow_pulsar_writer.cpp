#include "rf_pipelines_internals.hpp"
#include <rf_kernels/downsample.hpp>
#include <fstream>
#include <sstream>
#include <chrono>

#include <spshuff.hpp>

#include <ch_frb_io.hpp>

using namespace std;
using namespace rf_kernels;
using namespace ch_frb_io;
using namespace spshuff;

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
    const ssize_t ntime_max = 4096;
    const ssize_t nfreq_max = 16384;
    const ssize_t nsamp_max = nfreq_max * ntime_max;
    const ssize_t nbytes_charbuf_max = encode_ceil_bound(nsamp_max);

    if (nsamp_max % 8 != 0){
        throw runtime_error("rf_pipelines::chime_slow_pulsar_writer::chime_slow_pulsar_writer(): 8 must divide nsamp_max");
    }

    // one-time max size allocation of tmp buffers
    this->tmp_i = std::shared_ptr<std::vector<float>>(new std::vector<float>(nsamp_max));
    this->tmp_w = std::shared_ptr<std::vector<float>>(new std::vector<float>(nsamp_max));
    this->tmp_inorm = std::shared_ptr<std::vector<float>>(new std::vector<float>(ntime_max));
    this->tmp_qbuf = std::shared_ptr<std::vector<uint8_t>>(new std::vector<uint8_t>(ntime_max));
    this->tmp_var = std::shared_ptr<std::vector<float>>(new std::vector<float>(nfreq_max));
    this->tmp_mean = std::shared_ptr<std::vector<float>>(new std::vector<float>(nfreq_max));
    this->tmp_mask = std::shared_ptr<std::vector<uint8_t>>(new std::vector<uint8_t>(nsamp_max / 8));

    // a buffer to receive the huffman encoded data
    this->tmp_ibuf = std::shared_ptr<std::vector<uint32_t>>(new std::vector<uint32_t>(nbytes_charbuf_max/sizeof(uint32_t)));
}


// Caller must hold chunk lock!
void chime_slow_pulsar_writer::_get_new_chunk_with_lock()
{
    ssize_t beam_id, nbins;
    {
        std::lock_guard<std::mutex> lg(this->param_mutex);

        beam_id = this->beam_id;
        nbins = this->nbins;
    }

    // TODO: complete the initializer!
    std::shared_ptr<ch_chunk_initializer> ini_params(new ch_chunk_initializer(this->rt_state.memory_pool));
    ini_params->fpga_counts_per_sample = this->fpga_counts_per_sample;
    ini_params->frame0_nano = this->frame0_nano;
    this->wrote_start = false;
    this->chunk = slow_pulsar_chunk::make_slow_pulsar_chunk(ini_params);
    this->chunk->file_header.nbins = nbins;
    this->chunk->file_header.beam_id = beam_id;
    this->chunk->file_header.version = 4; // TODO shift this hard-coded constant elsewhere (e.g. to the chunk itself)

    // just to be sure...
    this->i0 = 0;
    this->bit0 = 0;
}


void chime_slow_pulsar_writer::init_real_time_state(const real_time_state &rt_state_)
{
    if (!rt_state_.memory_pool)
    throw runtime_error("rf_pipelines::chime_slow_pulsar_writer::init_real_time_state(): 'memory_pool' is an empty pointer, or uninitialized");
    if (!rt_state_.output_devices)
    throw runtime_error("rf_pipelines::chime_slow_pulsar_writer::init_real_time_state(): 'output_devices' is an empty pointer, or uninitialized");
    
    std::lock_guard<std::mutex> lg(this->chunk_mutex);
    this->rt_state = rt_state_;
    this->_get_new_chunk_with_lock();
}


void chime_slow_pulsar_writer::set_params(const ssize_t beam_id, const ssize_t nfreq, const ssize_t ntime, const ssize_t nbins, std::shared_ptr<std::string> base_path)
{   
    std::lock_guard<std::mutex> lg(this->param_mutex);

    // forgo checks on path validity for now
    this->base_path = base_path;

    this->beam_id = beam_id;

    this->nfreq_out = nfreq;
    this->ntime_out = ntime;
    this->nds_freq = this->nfreq/this->nfreq_out;
    this->nds_time = nt_chunk/this->ntime_out;

    this->nsamp = this->nfreq_out * this->ntime_out;

    // TODO check parameter consistency with existing file_header?
    this->nbins = nbins;

    if (this->nbins != -1 && this->nbins != nbins){
        throw runtime_error("rf_pipelines::chime_slow_pulsar_writer::set_params(): nbins must not change over the course of a run");
    }
    if (nbins != 5){
        throw runtime_error("rf_pipelines::chime_slow_pulsar_writer::set_params(): nbins must equal 5; email aroman@perimeterinstitute.ca if you desire a different binning");
    }
    if (this->ntime_out % 8 != 0){
        throw runtime_error("rf_pipelines::chime_slow_pulsar_writer::set_params(): 8 must divide ntime_out");
    }
    if (this->ntime_out > 4096){
        throw runtime_error("rf_pipelines::chime_slow_pulsar_writer::set_params(): ntime_out must be <= 4096");
    }

    this->downsampler = std::shared_ptr<rf_kernels::wi_downsampler>(new rf_kernels::wi_downsampler(
                                                         this->nds_freq, this->nds_time));
    // TODO: force a flush of the existing chunk/get new chunk?? Resolve first chunk header problems
}


template<typename T>
T* get_ptr(std::shared_ptr<std::vector<T>> vect)
{
    return &((*vect)[0]);
};


// this is a hack to get around sp_header forward declaration in the hpp file
// NOTE: this is only ever called from a process that current holds the of_mutex lock!
// thus, no additional locking is necessary. It would be dangerous to call this otherwise
void store_with_lock(chime_slow_pulsar_writer* spw, std::shared_ptr<std::string> base_path,
                     std::shared_ptr<sp_chunk_header> spch)
{
    // attempt to commit the chunk to slab
    const int result = spw->chunk->commit_chunk(spch, spw->tmp_ibuf, spw->compressed_data_len, 
                                    spw->tmp_mask, spw->tmp_mean, spw->tmp_var);

    if(result == 1){
        // the slab is full; enqueue write request
        std::shared_ptr<write_chunk_request> req(new write_chunk_request());
        // req->filename = "/home/aroman/tmp/test.dat";
        std::stringstream str;
        // TODO: add file sep end check?
        // TODO: implement directory structure
        str << *base_path << "/" << spw->ichunk << ".dat";
        req->filename = str.str();
        // std::cout << "writing output file " << req->filename <<std::endl;
        req->chunk = spw->chunk;
        spw->rt_state.output_devices->enqueue_write_request(req);
        double tend = spw->chunk->file_header.end;

        // TODO: address hackey constant below; not in the spirit of rf_pipelines
        double tstart = tend - 2560 * 1024 * spw->fpga_counts_per_sample * 1e-9;
        spw->_get_new_chunk_with_lock();
        // TODO: investigate suspicious nbins in file header
        spw->chunk->file_header.start = tstart;
        spw->chunk->file_header.end = tend;
        store_with_lock(spw, base_path, spch);
    }
    else if(result == 2){
        throw new runtime_error("chime_slow_pulsar_writer: byte size of single chunk exceeds memory slab capacity");
    }
}


// virtual override
void chime_slow_pulsar_writer::_process_chunk(float *intensity, ssize_t istride, float *weights, ssize_t wstride, ssize_t pos)
{
    int nfreq_out, ntime_out;
    {
        std::lock_guard<std::mutex> lg0(this->param_mutex);
        nfreq_out = this->nfreq_out;
        ntime_out = this->ntime_out;
    }


    if( (nfreq_out > 0) && (ntime_out > 0) ){
        // TODO: protect with separate lock relavant to file parameters
        // i.e. move back to a two-mutex approach; one for the chunk,
        // one for writer parameters

        int nds_time;
        std::shared_ptr<rf_kernels::wi_downsampler> downsampler;
        std::shared_ptr<sp_chunk_header> sph(new sp_chunk_header());
        std::shared_ptr<std::string> base_path;
        {
            std::lock_guard<std::mutex> lg1(this->param_mutex);

            nds_time = this->nds_time;        
            sph->ntime = ntime_out;
            sph->nfreq = nfreq_out;
            downsampler = this->downsampler;
            base_path = this->base_path;
        }

        sph->fpgaN = fpga_counts_per_sample;
        sph->fpga0 = this->initial_fpga_count + pos * this->fpga_counts_per_sample;
        sph->frame0_nano = this->frame0_nano;
        uint64_t fpga_counts_per_sample = this->fpga_counts_per_sample;
        uint64_t nsamp_chunk = this->nt_chunk;

        std::lock_guard<std::mutex> lg2(this->chunk_mutex);

        const ssize_t fpga_nano = 2560;
        const double tnow = (sph->fpga0 * fpga_nano + frame0_nano) * 1e-9;
        const double tend = tnow + 1e-9 * fpga_nano * (1024 * fpga_counts_per_sample);
        if(nt_chunk != 1024){
            throw new runtime_error("chime_slow_pulsar_writer: expects nt_chunk = 1024");
        }

        if(!wrote_start){
            chunk->file_header.start = tnow;
            wrote_start = true;
        }

        chunk->file_header.end = tend;

        bit0 = 0;
        i0 = 0;

        // TODO: check for no downsampling, rewrite to be cache-local
        downsampler->downsample(nfreq_out, ntime_out, 
                                        get_ptr<float>(tmp_i), ntime_out,
                                        get_ptr<float>(tmp_w), ntime_out,
                                        intensity, istride,
                                        weights, wstride);
        
        const ssize_t nrow_mask = ntime_out / 8;
        uint8_t mask_byte = 0;
        // estimate channel mean and var, compute mask
        for(ssize_t ifreq = 0; ifreq < nfreq_out; ifreq++){
            float s1 = 0.;
            float s2 = 0.;
            ssize_t ibyte_w = 0;
            ssize_t ibit_w = 0;
            for(ssize_t itime = 0; itime < ntime_out; itime++){
                float v = (*tmp_i)[ifreq * ntime_out + itime];
                s1 += v;
                s2 += v*v;

                // this is a particularly inelegant solution
                float w = (*tmp_w)[ifreq * ntime_out + itime] == 1.;
                if(w == 1.){
                    mask_byte += pow(2, ibit_w);
                }
                else if(w == 0.){
                    // pass
                }
                else{
                    throw new runtime_error("chime_slow_pulsar_writer: invalid value encountered in weight array");
                }

                ibit_w++;
                if(ibit_w == 8){
                    ibit_w = 0;
                    ibyte_w += 1;
                    (*tmp_mask)[ifreq * nrow_mask + ibyte_w] = mask_byte;
                    mask_byte = 0;
                }
            }
            
            // compute the relevant statistics
            float fmean = s1 / float(ntime_out);
            float f2mean = s2 / float(ntime_out);
            float fvar = (float(ntime_out) / float(ntime_out -1)) * (f2mean - fmean * fmean);

            (*tmp_mean)[ifreq]= fmean;
            (*tmp_var)[ifreq] = fvar;
            float stdev = sqrt(fvar);

            // second time pass; normalize samples
            for(ssize_t itime = 0; itime < ntime_out; itime++){
                (*tmp_inorm)[itime] = ((*tmp_i)[ifreq * ntime_out + itime] - fmean) / stdev;
            }

            // quantize just one row
            quantize_naive5(get_ptr<float>(tmp_inorm), get_ptr<uint8_t>(tmp_qbuf), ntime_out);

            // compress just one row, add to buffer
            // this should track the intra-bit state perfectly and yield a contiguous nsamp huffman coded array
            huff_encode_kernel(get_ptr<uint8_t>(tmp_qbuf), 
                               get_ptr<uint32_t>(tmp_ibuf), ntime_out, i0, bit0);
        }

        // set compressed data length
        compressed_data_len = i0;
        if(bit0 > 0){
            compressed_data_len++;
        }

        // TODO: give quantize_store a more hands-off role; quantization and compression happens in loop
        store_with_lock(this, base_path, sph);
        ichunk += 1;
    }
}


// virtual override
void chime_slow_pulsar_writer::_start_pipeline(Json::Value &json_attrs)
{
    // if (!json_attrs.isMember("frame0_nano") || !json_attrs.isMember("fpga_counts_per_sample")
    //     || !json_attrs.isMember("initial_fpga_count"))
    if (!json_attrs.isMember("initial_fpga_count") || !json_attrs.isMember("fpga_counts_per_sample"))
        throw runtime_error("chime_slow_pulsar_writer: expected json_attrs to contain members 'frame0_nano' and 'fpga_counts_per_sample'");
    
    // this->frame0_nano = json_attrs["frame0_nano"].asUInt64();
    this->frame0_nano = 0;
    this->fpga_counts_per_sample = json_attrs["fpga_counts_per_sample"].asUInt64();
    this->initial_fpga_count = json_attrs["initial_fpga_count"].asUInt64();
}


// virtual override
void chime_slow_pulsar_writer::_end_pipeline(Json::Value &json_output)
{
    // TODO: ensure flush of incomplete chunk (?)
}


// virtual override
Json::Value chime_slow_pulsar_writer::jsonize() const
{
    Json::Value ret;
    
    // TODO: expand on this?
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
