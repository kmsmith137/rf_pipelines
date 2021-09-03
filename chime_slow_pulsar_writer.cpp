#include "rf_pipelines_internals.hpp"
#include <rf_kernels/downsample.hpp>
#include <fstream>
#include <sstream>
#include <chrono>

#include <spshuff.hpp>
#include <ch_frb_io.hpp>

#include <immintrin.h>

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
    const ssize_t nbytes_n5_encoder = n5_encoder::min_nbytes(nsamp_max);

    if (nsamp_max % 8 != 0){
        throw runtime_error("rf_pipelines::chime_slow_pulsar_writer::chime_slow_pulsar_writer(): 8 must divide nsamp_max");
    }

    if (nsamp_max % 32 != 0){
        throw runtime_error("rf_pipelines::chime_slow_pulsar_writer::chime_slow_pulsar_writer(): 32 must divide nsamp_max");
    }

    // one-time max size allocation of tmp buffers
    this->tmp_intrin = std::shared_ptr<uint32_t>(aligned_alloc<uint32_t>(8, 32, false), free);
    this->tmp_intrinf1 = std::shared_ptr<float>(aligned_alloc<float>(8, 32, false), free);
    this->tmp_intrinf2 = std::shared_ptr<float>(aligned_alloc<float>(8, 32, false), free);

    this->tmp_i = std::shared_ptr<float>(aligned_alloc<float>(nsamp_max, 32, false), free);
    this->tmp_w = std::shared_ptr<float>(aligned_alloc<float>(nsamp_max, 32, false), free);

    this->tmp_inorm = std::shared_ptr<float>(aligned_alloc<float>(ntime_max, 32, false), free);
    this->tmp_qbuf = std::shared_ptr<uint8_t>(aligned_alloc<uint8_t>(ntime_max, 32, false), free);

    this->tmp_var = std::shared_ptr<std::vector<float>>(new std::vector<float>(nfreq_max));
    this->tmp_mean = std::shared_ptr<std::vector<float>>(new std::vector<float>(nfreq_max));
    this->tmp_mask = std::shared_ptr<std::vector<uint8_t>>(new std::vector<uint8_t>(nsamp_max / 8));

    // a buffer to receive the huffman encoded data
    // this->tmp_ibuf = std::shared_ptr<std::vector<uint32_t>>(new std::vector<uint32_t>(nbytes_charbuf_max/sizeof(uint32_t)));
    this->tmp_ibuf = make_shared<std::vector<uint8_t>>(nbytes_n5_encoder);
}


// Caller must hold chunk lock!
void chime_slow_pulsar_writer::_get_new_chunk_with_locks(const ssize_t beam_id, const ssize_t nbins, const uint64_t frame0_nano)
{
    // TODO: complete the initializer!
    std::shared_ptr<ch_chunk_initializer> ini_params = make_shared<ch_chunk_initializer>(this->rt_state.memory_pool);
    ini_params->fpga_counts_per_sample = this->fpga_counts_per_sample;
    ini_params->frame0_nano = frame0_nano;
    this->wrote_start = false;
    std::cout << "making new chunk" << std::endl;
    this->chunk = slow_pulsar_chunk::make_slow_pulsar_chunk(ini_params);
    this->chunk->file_header.nbins = nbins;
    this->chunk->file_header.beam_id = beam_id;
    this->chunk->file_header.version = 5; // TODO shift this hard-coded constant elsewhere (e.g. to the chunk itself)
}


void chime_slow_pulsar_writer::init_real_time_state(const real_time_state &rt_state_)
{
    if (!rt_state_.memory_pool)
    throw runtime_error("rf_pipelines::chime_slow_pulsar_writer::init_real_time_state(): 'memory_pool' is an empty pointer, or uninitialized");
    if (!rt_state_.output_devices)
    throw runtime_error("rf_pipelines::chime_slow_pulsar_writer::init_real_time_state(): 'output_devices' is an empty pointer, or uninitialized");

    this->rt_state = rt_state_;


}


void chime_slow_pulsar_writer::set_params(const ssize_t beam_id, const ssize_t nfreq_out, 
                const ssize_t ntime_out, const ssize_t nbins, std::shared_ptr<std::string> base_path,
                const uint64_t frame0_nano)
{ 
    std::lock_guard<std::mutex> lg_param(this->param_mutex);

    // this is interpreted as an off switch
    bool is_null = (nfreq_out == 0) || (ntime_out == 0);

    if(is_null){
        pstate = nullptr;
    }
    else{
        pstate = make_shared<chime_slow_pulsar_writer::param_state>();
        // see if we require a new chunk
        if(!chunk || chunk->file_header.nbins != nbins || chunk->file_header.beam_id != beam_id){
            std::lock_guard<std::mutex> lg_chunk(this->chunk_mutex);
            this->_get_new_chunk_with_locks(beam_id, nbins, frame0_nano);
        }
        
        // forgo checks on path validity for now
        pstate->base_path = base_path;

        pstate->beam_id = beam_id;

        // std::shared_ptr<ch_frb_io::assembled_chunk> achunk = stream->get_assembled_chunk();
        pstate->frame0_nano = frame0_nano;

        pstate->nfreq_out = nfreq_out;
        pstate->ntime_out = ntime_out;
        // TODO: remove hard-coded upsampling factor (16)
        pstate->nds_freq = nfreq / pstate->nfreq_out;
        pstate->nds_time = nt_chunk / pstate->ntime_out;

        cout << "Creating downsampler.  nfreq " << nfreq << ", nfreq_out " << pstate->nfreq_out
             << ", nt_chunk " << nt_chunk << ", ntime_out " << pstate->ntime_out
             << ", nds_freq " << pstate->nds_freq << ", nds_time " << pstate->nds_time << endl;

        pstate->downsampler = std::shared_ptr<rf_kernels::wi_downsampler>(new rf_kernels::wi_downsampler(
                                                   pstate->nds_freq, pstate->nds_time));

        pstate->nsamp = pstate->nfreq_out * pstate->ntime_out;

        // TODO check parameter consistency with existing file_header?
        pstate->nbins = nbins;

        if (pstate->nbins != -1 && pstate->nbins != nbins){
            throw runtime_error("rf_pipelines::chime_slow_pulsar_writer::set_params(): nbins must not change over the course of a run");
        }
        if (nbins != 5){
            throw runtime_error("rf_pipelines::chime_slow_pulsar_writer::set_params(): nbins must equal 5; email aroman@perimeterinstitute.ca if you desire a different binning");
        }
        if (pstate->nfreq_out % 2 != 0){
            throw runtime_error("rf_pipelines::chime_slow_pulsar_writer::set_params(): nfreq_out must be a power of 2");
        }   
        if (pstate->ntime_out % 2 != 0){
            throw runtime_error("rf_pipelines::chime_slow_pulsar_writer::set_params(): ntime_out must be a power of 2");
        }   
        if (pstate->ntime_out % 8 != 0){
            throw runtime_error("rf_pipelines::chime_slow_pulsar_writer::set_params(): 8 must divide ntime_out");
        }    
        if (pstate->ntime_out % 64 != 0){
            throw runtime_error("rf_pipelines::chime_slow_pulsar_writer::set_params(): 64 must divide ntime_out");
        }
        if (pstate->ntime_out > 4096){
            throw runtime_error("rf_pipelines::chime_slow_pulsar_writer::set_params(): ntime_out must be <= 4096");
        }
        if (pstate->nsamp % 32 != 0){
            throw runtime_error("rf_pipelines::chime_slow_pulsar_writer::set_params(): 32 must divide nsamp (simd alignment)");
        }
        // TODO: force a flush of the existing chunk/get new chunk?? Resolve first chunk header problems
    }
}


template<typename T>
T* get_ptr(std::shared_ptr<std::vector<T>> vect)
{
    return &((*vect)[0]);
};

template<typename T>
T* get_ptr(std::shared_ptr<T> in_ptr)
{
    return &(*in_ptr);
};


std::string pad_month(int tm_mon){
    std::string ret("");
    int disp_mon = tm_mon + 1;

    if(disp_mon < 10){
        ret += "0";
    }

    ret += std::to_string(disp_mon);
    return ret;
}

std::string pad_day(int tm_mday){
    std::string ret("");

    if(tm_mday < 10){
        ret += "0";
    }

    ret += std::to_string(tm_mday);
    return ret;
}

std::string pad_beam(int ibeam){
    std::string beam_str("");

    int trunc_log = max((int) log10(ibeam), 0);
    for(int i = 0; i < 3 - trunc_log; i++){
        beam_str += "0";
    }
    beam_str += std::to_string(ibeam);
    return beam_str;
}

// this is a hack to get around sp_header forward declaration in the hpp
// NOTE: this is only ever called from a process that current holds the chunk_mutex lock!
// thus, no additional locking is necessary; it would be dangerous to call this otherwise
void store_with_lock(chime_slow_pulsar_writer* spw, 
                     std::shared_ptr<chime_slow_pulsar_writer::param_state> pstate,
                     std::shared_ptr<sp_chunk_header> spch)
{
    // attempt to commit the chunk to slab
    const int result = spw->chunk->commit_chunk(spch, spw->tmp_ibuf, spw->compressed_data_len, 
                                    spw->tmp_mask, spw->tmp_mean, spw->tmp_var);

    if(result == 1){
        // the slab is full; enqueue write request
        std::shared_ptr<write_chunk_request> req = make_shared<write_chunk_request>();
        // req->filename = "/home/aroman/tmp/test.dat";
        std::shared_ptr<std::stringstream> pathstr = make_shared<std::stringstream>();
        // TODO: add file sep end check?
        // // TODO: implement directory structure
        std::time_t utc_start = (std::time_t) spw->chunk->file_header.start;
        std::tm* tnow = std::gmtime(&utc_start);

        *pathstr << *(pstate->base_path) << "/" << (tnow->tm_year + 1900) << "/" << pad_month(tnow->tm_mon) 
        << "/" << pad_day(tnow->tm_mday) << "/" << pad_beam(pstate->beam_id) << "/" <<
        ((ssize_t) spw->chunk->file_header.start) << "_" << ((ssize_t) spw->chunk->file_header.end) << ".dat";

        // *pathstr << *(pstate->base_path) << "/" << spw->ichunk << ".dat";

        req->filename = pathstr->str();
        // std::cout << "writing output file " << req->filename <<std::endl;
        req->chunk = spw->chunk;
        spw->rt_state.output_devices->enqueue_write_request(req);

        // double set times by accident??
        // double tend = spw->chunk->file_header.end;
        // // TODO: address hackey constant below; not in the spirit of rf_pipelines
        // double tstart = tend - 2560 * 1024 * spw->fpga_counts_per_sample * 1e-9;
        spw->_get_new_chunk_with_locks(pstate->beam_id, pstate->nbins, pstate->frame0_nano);
        // TODO: investigate suspicious nbins in file header
        // spw->chunk->file_header.start = tstart;
        // spw->chunk->file_header.end = tend;
        store_with_lock(spw, pstate, spch);
    }
    else if(result == 2){
        throw new runtime_error("chime_slow_pulsar_writer: byte size of single chunk exceeds memory slab capacity");
    }
}

// useful for timing
void dummy_compress(const uint8_t* qbuf, uint32_t* cbuf, const ssize_t nsamp, 
                    ssize_t& i0, ssize_t& bit0)
{
    // based on max write
    const ssize_t rbit = (3 * nsamp) % 32;
    ssize_t write_len = (3 * nsamp) / 32;
    if(rbit > 0){
        write_len += 1;
    }

    uint8_t s = 0;
    for(ssize_t i = 0; i < nsamp; i++){
        s += qbuf[i];
    }

    for(ssize_t i = i0; i < write_len; i++){
        cbuf[i] = (uint32_t) s;
    }

    i0 += write_len;
}


const double get_time_sec()
{
    //const
    auto tnow = std::chrono::system_clock::now();
    return std::chrono::duration_cast<std::chrono::seconds>(tnow.time_since_epoch()).count();
}

inline struct timeval xgettimeofday()
{
    struct timeval tv;

    int err = gettimeofday(&tv, NULL);
    if (_unlikely(err))
	throw std::runtime_error("gettimeofday failed");

    return tv;
}
inline int64_t usec_between(const struct timeval &tv1, const struct timeval &tv2)
{
    return 1000000 * int64_t(tv2.tv_sec - tv1.tv_sec) + int64_t(tv2.tv_usec - tv1.tv_usec);
}

// virtual override
void chime_slow_pulsar_writer::_process_chunk(float *intensity, ssize_t istride, float *weights, ssize_t wstride, ssize_t pos)
{
    //auto dstn_t0 = std::chrono::steady_clock::now();
    struct timeval dstn_t0 = xgettimeofday();

    std::shared_ptr<sp_chunk_header> sph = make_shared<sp_chunk_header>();
    std::shared_ptr<chime_slow_pulsar_writer::param_state> pstate;

    {        
        // TODO: protect with separate lock relavant to file parameters
        // i.e. move back to a two-mutex approach; one for the chunk,
        // one for writer parameters
        std::lock_guard<std::mutex> lg0(this->param_mutex);
        // copy local parameter state
        pstate = this->pstate;
    }

    // check whether the pstate is null (acquisition off)
    // do not proceed or block until chunk is initialized
    if(pstate && chunk){

        // redundant copy to the chunk header; old way of tracking state
        sph->ntime = pstate->ntime_out;
        sph->nfreq = pstate->nfreq_out;

        n5_encoder enc(get_ptr<uint8_t>(this->tmp_ibuf));

        sph->fpgaN = fpga_counts_per_sample;
        sph->fpga0 = this->initial_fpga_count + pos * this->fpga_counts_per_sample;
        sph->frame0_nano = pstate->frame0_nano;
        uint64_t fpga_counts_per_sample = this->fpga_counts_per_sample;
        uint64_t nsamp_chunk = this->nt_chunk;

        const ssize_t fpga_nano = 2560;

        // const double tnow = get_time_sec();
        // std::cout << "new frame0_nano ";
        const double tnow = (sph->fpga0 * fpga_nano + pstate->frame0_nano) * 1e-9;
        // std::cout << tnow << std::endl;
        const double tend = tnow + 1e-9 * fpga_nano * (nt_chunk * fpga_counts_per_sample);

        // make less rigid?
        if(nt_chunk != 1024){
            throw new runtime_error("chime_slow_pulsar_writer: expects nt_chunk = 1024");
        }

        {
            std::lock_guard<std::mutex> lg_chunk0(this->chunk_mutex);
            if(!wrote_start){
                chunk->file_header.start = tnow;
                wrote_start = true;
            }

            chunk->file_header.end = tend;
        }

        bit0 = 0;
        i0 = 0;
        uint32_t tmp = 0; // tmp buffer for huffman encoding

        const bool no_ds = (pstate->nds_time == 1) && (pstate->nds_freq == 1);

        auto t0 = std::chrono::high_resolution_clock::now();

        float* ds_i;
        float* ds_w;
        // ssize_t my_istride;
        // ssize_t my_wstride;
        // if(no_ds){
        //     // no-op
        //     ds_i = intensity;
        //     ds_w = weights;
        //     my_istride = istride;
        //     my_wstride = wstride;
        // }
        // else{
        //     pstate->downsampler->downsample(pstate->nfreq_out, pstate->ntime_out, 
        //                         get_ptr<float>(tmp_i), pstate->ntime_out,
        //                         get_ptr<float>(tmp_w), pstate->ntime_out,
        //                         intensity, istride,
        //                         weights, wstride);

        //     ds_i = get_ptr<float>(tmp_i);
        //     ds_w = get_ptr<float>(tmp_w);
        //     my_istride = pstate->ntime_out;
        //     my_wstride = pstate->ntime_out;
        // }

        // const float* ds_ic = ds_i;
        // const float* ds_wc = ds_w;
        // const ssize_t istridec = my_istride;
        // const ssize_t wstridec = my_wstride;


        // TODO: check for no downsampling, rewrite to be cache-local

        auto t1 = std::chrono::high_resolution_clock::now();

        const ssize_t nrow_mask = pstate->ntime_out / 8;
        const ssize_t nds_freq = pstate->nds_freq;
        const ssize_t nds_time = pstate->nds_time;
        const ssize_t nds_tot = nds_freq * nds_time;
        const ssize_t log_nds_tot = log2l(nds_tot);
        const ssize_t nfreq_out = pstate->nfreq_out;
        const ssize_t ntime_out = pstate->ntime_out;

        auto t2 = std::chrono::high_resolution_clock::now();

        double wdur = 0.;
        double ndur = 0.;
        double qdur = 0.;
        double cdur = 0.;

        const float* weightc = weights;
        const float* intensityc = intensity;
        uint32_t* tmp0 = get_ptr<uint32_t>(tmp_intrin);
        float* tmp1 = get_ptr<float>(tmp_intrinf1);
        float* tmp2 = get_ptr<float>(tmp_intrinf2);
        float fnorm;

        // estimate channel mean and var, compute mask
        for(ssize_t ifreq = 0; ifreq < pstate->nfreq_out; ifreq++){

            if(!no_ds){
                // TODO check separately for downsampling in each axis
                pstate->downsampler->downsample(1, pstate->ntime_out, 
                    get_ptr<float>(tmp_i), pstate->ntime_out,
                    get_ptr<float>(tmp_w), pstate->ntime_out,
                    intensity + nds_freq * ifreq * istride, istride,
                    weights + nds_freq * ifreq * wstride, wstride);

                ds_i = get_ptr<float>(tmp_i);
                ds_w = get_ptr<float>(tmp_w);
            }
            else{
                ds_i = intensity + ifreq * istride;
                ds_w = weights + ifreq * wstride;
            }

            const float* ds_ic = ds_i;
            const float* ds_wc = ds_w;

            float s1 = 0.;
            float s2 = 0.;
            ssize_t ibyte_w = 0;
            ssize_t ibit_w = 0;
            uint8_t* mask_tmp_ar = get_ptr<uint8_t>(tmp_mask);

            // auto t21 = std::chrono::high_resolution_clock::now();
            // for(ssize_t itime = 0; itime < pstate->ntime_out; itime++){
            //     const float v = ds_ic[ifreq * istridec + itime];
            //     s1 += v;
            //     s2 += v*v;

            //     // this is a particularly inelegant solution
            //     // TODO: FIX
            //     float w = ds_wc[ifreq * wstridec + itime];
            //     // std::cout << w << std::endl;
            //     mask_byte += ((uint8_t) w) << ibit_w;

            //     ibit_w++;
            //     if(ibit_w == 8){
            //         ibit_w = 0;
            //         ibyte_w += 1;
            //         (*tmp_mask)[ifreq * nrow_mask + ibyte_w] = mask_byte;
            //         mask_byte = 0;
            //     }
            // }

            // for(ssize_t itime_o = 0; itime_o < pstate->ntime_out; itime_o+=8){
            //     uint8_t mask_byte = 0;
            //     const ssize_t iind = ifreq * istridec + itime_o;
            //     const ssize_t wind = ifreq * wstridec + itime_o;
            //     for(ssize_t itime_i = 0; itime_i < 8; itime_i++){
            //         const float v = ds_ic[iind + itime_i];
            //         s1 += v;
            //         s2 += v*v;

            //         // this is a particularly inelegant solution
            //         // TODO: FIX
            //         const float w = ds_wc[wind + itime_i];
            //         // std::cout << w << std::endl;
            //         mask_byte += ((uint8_t) w) << itime_i;
            //     }
            //     mask_tmp_ar[ifreq * nrow_mask + itime_o] = mask_byte;
            // }
            __m256 ms1 = _mm256_set1_ps(0.);
            __m256 ms2 = _mm256_set1_ps(0.);
            // const __m128i shift0 = _mm_set_epi32(1,1,1,1);
            const __m256i shift_ds = _mm256_set_epi32(log_nds_tot, log_nds_tot, log_nds_tot, log_nds_tot,
                                                      log_nds_tot, log_nds_tot, log_nds_tot, log_nds_tot);
            const __m256i shift0 = _mm256_set_epi32(1,1,1,1,1,1,1,1);
            const __m256i shift1 = _mm256_set_epi32(2,2,2,2,2,2,2,2);
            const __m256i shift2 = _mm256_set_epi32(4,4,4,4,4,4,4,4);
            for(ssize_t iframe_o = 0; iframe_o < ntime_out/8; iframe_o+=1){
                const ssize_t itime_o = iframe_o * 8;
                const __m256 mvari = _mm256_load_ps(ds_ic + itime_o);
                ms1 = _mm256_add_ps(mvari, ms1);
                ms2 = _mm256_fmadd_ps(mvari, mvari, ms2);

                // compute mask
                __m256i mvarw = _mm256_cvtps_epi32(_mm256_load_ps(ds_wc + itime_o));
                mvarw = _mm256_srlv_epi32(mvarw, shift_ds); // divide by downsampling factor
                mvarw = _mm256_add_epi32(mvarw, _mm256_shuffle_epi32(_mm256_sllv_epi32(mvarw, shift0), 177));
                mvarw = _mm256_add_epi32(mvarw, _mm256_shuffle_epi32(_mm256_sllv_epi32(mvarw, shift1), 2));
                mvarw = _mm256_add_epi32(mvarw, _mm256_sllv_epi32(_mm256_permute2f128_si256(mvarw, mvarw, 1), shift2));

                _mm256_store_si256((__m256i*) tmp0, mvarw);
                // uint8_t mask_byte = 0;
                // for(ssize_t itime_i = 0; itime_i < 8; itime_i++){
                    // const uint8_t w = ((uint8_t) ds_wc[itime_o + itime_i]) / (nds_tot);
                    // mask_byte += ((uint8_t) w) << itime_i;
                // }
                // std::cout << std::endl;

                mask_tmp_ar[ifreq * nrow_mask + iframe_o] = (uint8_t) tmp0[0];
                // mask_tmp_ar[ifreq * nrow_mask + itime_o] = mask_byte;
                // std::bitset<8> x((uint8_t) tmp_intrin[0]);
                // std::bitset<8> y((uint8_t) mask_byte);
                // // std::cout << "wut2" << endl;
                // std::cout << (((uint8_t) tmp_intrin[0]) - mask_byte) << " " << x << " " << y << " " << ((uint32_t) ((uint8_t) tmp_intrin[0])) << " " << ((uint32_t) mask_byte) << std::endl;
            }

            _mm256_store_ps(tmp1, ms1);
            _mm256_store_ps(tmp2, ms2);

            for(ssize_t i = 0; i < 8; i++){
                s1 += tmp1[i];
                s2 += tmp2[i];
            }

            // auto t22 = std::chrono::high_resolution_clock::now();

            // std::chrono::duration<double> tmp1 = t22 - t21;
            // wdur += tmp1.count();

            // compute the relevant statistics
            const float fmean = s1 / float(pstate->ntime_out);
            const float f2mean = s2 / float(pstate->ntime_out);
            const float fvar = (float(pstate->ntime_out) / float(pstate->ntime_out -1)) * (f2mean - fmean * fmean);

            (*tmp_mean)[ifreq]= fmean;
            (*tmp_var)[ifreq] = fvar;
            const float stdev = sqrt(fvar);
            // const float wdev = 1./stdev;
            // const float meanc = fmean/stdev;

            enc.set_mean_and_rms(fmean, stdev);

            for(ssize_t itime = 0; itime < pstate->ntime_out; itime += 64){
                enc.encode64_aligned(ds_ic + itime);
            }

            // auto t221 = std::chrono::high_resolution_clock::now();
            // float* tmp_inorm_ptr = get_ptr<float>(tmp_inorm);
            // // second time pass; normalize samples
            // // for(ssize_t itime = 0; itime < pstate->ntime_out; itime++){
            // //    tmp_inorm_ptr[itime] = (ds_ic[ifreq * istridec + itime] - fmean) * wdev;
            // // }
            // const __m256 mwdev = _mm256_broadcast_ss(&wdev);
            // const __m256 mmeanc = _mm256_broadcast_ss(&meanc);
            // const float* this_ic = ds_ic;

            // for(ssize_t itime = 0; itime < pstate->ntime_out; itime+=8){
            //     const __m256 mvar = _mm256_load_ps(this_ic + itime);
            //     _mm256_store_ps(tmp_inorm_ptr + itime, _mm256_fmsub_ps(mvar, mwdev, mmeanc));
            // }

            // auto t222 = std::chrono::high_resolution_clock::now();
            // std::chrono::duration<double> tmp21 = t222 - t221;
            // ndur += tmp21.count();

            // auto t23 = std::chrono::high_resolution_clock::now();
            // quantize just one row
            // quantize_naive5_simd4(get_ptr<float>(tmp_inorm), get_ptr<uint8_t>(tmp_qbuf), pstate->ntime_out);
            // auto t24 = std::chrono::high_resolution_clock::now();

            // std::chrono::duration<double> tmp2 = t24 - t23;
            // qdur += tmp2.count();
            // compress just one row, add to buffer
            // this should track the intra-bit state perfectly and yield a contiguous nsamp huffman coded array
            // huff_encode_kernel(get_ptr<uint8_t>(tmp_qbuf), 
            //                    get_ptr<uint32_t>(tmp_ibuf), pstate->ntime_out, i0, bit0, tmp);
            // dummy_compress(get_ptr<uint8_t>(tmp_qbuf), 
            //                    get_ptr<uint32_t>(tmp_ibuf), pstate->ntime_out, i0, bit0);
            // auto t25 = std::chrono::high_resolution_clock::now();

            // std::chrono::duration<double> tmp3 = t25 - t24;
            // cdur += tmp3.count();
        }

        // flush n5_encoder and set the length in 32 bit words
        compressed_data_len = (enc.flush() + 31) / 32;

        auto t3 = std::chrono::high_resolution_clock::now();

        // set compressed data length
        // compressed_data_len = i0;
        // if(bit0 > 0){
        //     compressed_data_len++;
        // }

        // TODO: give quantize_store a more hands-off role; quantization and compression happens in loop
        {
            std::lock_guard<std::mutex> lg_chunk1(this->chunk_mutex);
            store_with_lock(this, pstate, sph);
        }

        auto t4 = std::chrono::high_resolution_clock::now();
        ichunk += 1;
        std::chrono::duration<double> tdownsample = t1 - t0;
        std::chrono::duration<double> tloop = t3 - t2;
        std::chrono::duration<double> twrite = t4 - t3;
        std::cout << "time per exec: " << this->time_spent_in_transform / ichunk << std::endl;
        std::cout << "\tdownsample: " << tdownsample.count() << std::endl;
        std::cout << "\tloop: " << tloop.count() << std::endl;
        // std::cout << "\t\tweight_loop: " << wdur << std::endl;
        // std::cout << "\t\tnorm_loop: " << ndur << std::endl;
        // std::cout << "\t\tquantize_loop: " << qdur << std::endl;
        // std::cout << "\t\tcompress_loop: " << cdur << std::endl;
        std::cout << "\twrite: " << twrite.count() << std::endl;
    }


    //auto dstn_t1 = std::chrono::steady_clock::now();
    //double dt = (std::chrono::duration<double>(dstn_t1 - dstn_t0).count());
    struct timeval dstn_t1 = xgettimeofday();
    double dt = usec_between(dstn_t0, dstn_t1) * 1e-6;
    cout << "SPS writer took " << dt << " seconds" << endl;
}


// virtual override
void chime_slow_pulsar_writer::_start_pipeline(Json::Value &json_attrs)
{
    if (!json_attrs.isMember("initial_fpga_count") || !json_attrs.isMember("fpga_counts_per_sample"))
        throw runtime_error("chime_slow_pulsar_writer: expected json_attrs to contain members 'frame0_nano' and 'fpga_counts_per_sample'");
    
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

    //ssize_t beam_id;
    ssize_t nfreq_out = 0;
    ssize_t ntime_out = 0;
    ssize_t nbins = 0;
    string base_path = "";
    string key = "nfreq_out";
    if (j.isMember(key))
        nfreq_out = j[key].asInt();
    key = "ntime_out";
    if (j.isMember(key))
        ntime_out = j[key].asInt();
    key = "nbins";
    if (j.isMember(key))
        nbins = j[key].asInt();
    key = "base_path";
    if (j.isMember(key))
        base_path = j[key].asString();

    auto sps = make_shared<chime_slow_pulsar_writer> (nt_chunk);

    if ((base_path.size() > 0) && (nfreq_out > 0) && (ntime_out > 0) && (nbins > 0))
	cout << "FIXME: ignoring initial_params for now" << endl;

    return sps;
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
