#include <algorithm>
#include "rf_pipelines_internals.hpp"

#ifdef HAVE_CH_FRB_IO
#include <ch_frb_io.hpp>
#endif

using namespace std;

namespace rf_pipelines {
#if 0
};  // pacify emacs c-mode!
#endif


#ifndef HAVE_CH_FRB_IO

shared_ptr<wi_stream> make_chime_network_stream(const shared_ptr<ch_frb_io::intensity_network_stream> &stream, int assembler_index, float prescale)
{
    throw runtime_error("rf_pipelines::make_chime_network_stream() was called, but rf_pipelines was compiled without ch_frb_io");
}

shared_ptr<wi_stream> make_chime_network_stream(int udp_port, int assembler_index, float prescale)
{
    throw runtime_error("rf_pipelines::make_chime_network_stream() was called, but rf_pipelines was compiled without ch_frb_io");
}

shared_ptr<wi_stream> make_dummy_chime_network_stream(ssize_t nt_tot, int nupfreq, int nt_per_packet, int fpga_counts_per_sample, double pool_gb)
{
    throw runtime_error("rf_pipelines::make_dummy_chime_network_stream() was called, but rf_pipelines was compiled without ch_frb_io");
}

#else  // HAVE_CH_FRB_IO


// -------------------------------------------------------------------------------------------------
//
// chime_network_stream


struct chime_network_stream : public wi_stream
{
    shared_ptr<ch_frb_io::intensity_network_stream> stream;
    const float prescale;
    const int assembler_id;

    chime_network_stream(const shared_ptr<ch_frb_io::intensity_network_stream> &stream_, int assembler_index_, float prescale_);
    virtual ~chime_network_stream() { }

    virtual bool _fill_chunk(float *intensity, ssize_t istride, float *weights, ssize_t wstride, ssize_t pos) override;
    virtual void _bind_stream(Json::Value &j) override;
    virtual void _start_pipeline(Json::Value &j) override;
    virtual void _end_pipeline(Json::Value &j) override;
};


chime_network_stream::chime_network_stream(const shared_ptr<ch_frb_io::intensity_network_stream> &stream_, int assembler_index_, float prescale_) :
    wi_stream("chime_network_stream"),
    stream(stream_), 
    prescale(prescale_),
    assembler_id(assembler_index_)
{ 
    if (!stream)
	throw runtime_error("rf_pipelines: empty stream pointer passed to chime_network_stream constructor");
    if (assembler_id < 0)
	throw runtime_error("chime_network_stream constructor: assembler_index must be >= 0");

    this->nfreq = ch_frb_io::constants::nfreq_coarse_tot * stream->ini_params.nupfreq;
    this->nt_chunk = ch_frb_io::constants::nt_per_assembled_chunk;
}


void chime_network_stream::_bind_stream(Json::Value &json_attrs)
{
    json_attrs["freq_lo_MHz"] = 400.0;
    json_attrs["freq_hi_MHz"] = 800.0;
    json_attrs["dt_sample"] = chime_file_stream_base::chime_seconds_per_fpga_count * stream->ini_params.fpga_counts_per_sample;
}


void chime_network_stream::_start_pipeline(Json::Value &j)
{
    // tells network thread to start reading packets, returns immediately
    stream->start_stream();

    stream->wait_for_first_packet();
    uint64_t fpga0 = stream->get_first_fpgacount();
    uint64_t frame0_nano = stream->get_frame0_nano();
    
    j["initial_fpga_count"] = Json::UInt64(fpga0);
    j["frame0_nano"] = Json::UInt64(frame0_nano);
    j["fpga_counts_per_sample"] = stream->ini_params.fpga_counts_per_sample;
}


bool chime_network_stream::_fill_chunk(float *intensity, ssize_t istride, float *weights, ssize_t wstride, ssize_t pos)
{
    shared_ptr<ch_frb_io::assembled_chunk> chunk = stream->get_assembled_chunk(assembler_id);
    if (!chunk)
	return false;

    rf_assert(this->nfreq == ch_frb_io::constants::nfreq_coarse_tot * chunk->nupfreq);
    chunk->decode(intensity, weights, istride, wstride, prescale);
    return true;
}


void chime_network_stream::_end_pipeline(Json::Value &j)
{
    stream->join_threads();
}


// -------------------------------------------------------------------------------------------------
//
// chime_dummy_network_stream


struct chime_dummy_network_stream : public wi_stream
{
    ssize_t nt_tot;
    int nupfreq;
    int nt_per_packet;
    int fpga_counts_per_sample;
    double pool_gb;

    vector<shared_ptr<ch_frb_io::assembled_chunk>> chunk_pool;
    ssize_t ichunk = 0;
    ssize_t nchunks;
    
    chime_dummy_network_stream(ssize_t nt_tot_, int nupfreq_, int nt_per_packet_, int fpga_counts_per_sample_, double pool_gb_) :
	wi_stream("chime_dummy_network_stream"),
	nt_tot(nt_tot_),
	nupfreq(nupfreq_),
	nt_per_packet(nt_per_packet_),
	fpga_counts_per_sample(fpga_counts_per_sample_),
	pool_gb(pool_gb_)
    {
	if (nt_tot < 0)
	    throw runtime_error("rf_pipelines::chime_dummy_network_stream constructor: expected nt_tot > 0");
	if (nupfreq <= 0)
	    throw runtime_error("rf_pipelines::chime_dummy_network_stream constructor: expected nupfreq > 0");
	if (nt_per_packet <= 0)
	    throw runtime_error("rf_pipelines::chime_dummy_network_stream constructor: expected nt_per_packet > 0");
	if (!is_power_of_two(nt_per_packet))
	    throw runtime_error("rf_pipelines::chime_dummy_network_stream constructor: expected nt_per_packet to be a power of two");
	if (fpga_counts_per_sample <= 0)
	    throw runtime_error("rf_pipelines::chime_dummy_network_stream constructor: expected fpga_counts_per_sample > 0");
	if (pool_gb < 0.0)
	    throw runtime_error("rf_pipelines::chime_dummy_network_stream constructor: expected pool_gb >= 0.0");
	if (pool_gb > 100.0)
	    throw runtime_error("rf_pipelines::chime_dummy_network_stream constructor: expected pool_gb <= 100.0");

        int ndownfreq = 1024;
	double gb_per_chunk = ch_frb_io::assembled_chunk::get_memory_slab_size(nupfreq, nt_per_packet, ndownfreq);

	this->nfreq = ch_frb_io::constants::nfreq_coarse_tot * nupfreq;
	this->nt_chunk = ch_frb_io::constants::nt_per_assembled_chunk;
	this->nchunks = ssize_t(pool_gb / gb_per_chunk) + 1;
	this->chunk_pool.resize(nchunks);
    }

    virtual void _bind_stream(Json::Value &json_attrs) override
    {
	json_attrs["freq_lo_MHz"] = 400.0;
	json_attrs["freq_hi_MHz"] = 800.0;
	json_attrs["dt_sample"] = chime_file_stream_base::chime_seconds_per_fpga_count * fpga_counts_per_sample;
    }

    virtual void _start_pipeline(Json::Value &json_attrs) override
    {
	json_attrs["initial_fpga_count"] = 0;
	json_attrs["fpga_counts_per_sample"] = this->fpga_counts_per_sample;
    }

    virtual void _allocate() override
    {
	std::random_device rd;
	std::mt19937 rng(rd());

	ch_frb_io::assembled_chunk::initializer ini_params;
	ini_params.nupfreq = this->nupfreq;
	ini_params.nt_per_packet = this->nt_per_packet;
	ini_params.fpga_counts_per_sample = this->fpga_counts_per_sample;

	for (ssize_t i = 0; i < nchunks; i++) {
	    this->chunk_pool[i] = make_shared<ch_frb_io::assembled_chunk> (ini_params);
	    this->chunk_pool[i]->randomize(rng);
	}
    }

    virtual bool _fill_chunk(float *intensity, ssize_t istride, float *weights, ssize_t wstride, ssize_t pos) override
    {
	auto chunk = chunk_pool[ichunk % nchunks];
	ichunk++;

	rf_assert(this->nfreq == ch_frb_io::constants::nfreq_coarse_tot * chunk->nupfreq);
	chunk->decode(intensity, weights, istride, wstride);

	return (pos < nt_tot);
    }
    
    virtual void _deallocate() override
    {
	for (ssize_t i = 0; i < nchunks; i++)
	    this->chunk_pool[i].reset();
    }
    
    virtual Json::Value jsonize() const override
    {
	Json::Value ret;
	ret["class_name"] = "chime_dummy_network_stream";
	ret["nt_tot"] = Json::Int64(nt_tot);
	ret["nupfreq"] = nupfreq;
	ret["nt_per_packet"] = nt_per_packet;
	ret["fpga_counts_per_sample"] = fpga_counts_per_sample;
	ret["pool_gb"] = pool_gb;
	return ret;
    }

    // Note: json deserialization defined later in file.
};


shared_ptr<wi_stream> make_dummy_chime_network_stream(ssize_t nt_tot, int nupfreq, int nt_per_packet, int fpga_counts_per_sample, double pool_gb)
{
    return make_shared<chime_dummy_network_stream> (nt_tot, nupfreq, nt_per_packet, fpga_counts_per_sample, pool_gb);
}


// -------------------------------------------------------------------------------------------------


shared_ptr<wi_stream> make_chime_network_stream(const shared_ptr<ch_frb_io::intensity_network_stream> &stream, int assembler_index, float prescale)
{
    return make_shared<chime_network_stream> (stream, assembler_index, prescale);
}


shared_ptr<wi_stream> make_chime_network_stream(int udp_port, int beam_id, float prescale)
{
    ch_frb_io::intensity_network_stream::initializer ini_params;
    ini_params.nbeams = 1;

    if (udp_port > 0)
	ini_params.udp_port = udp_port;

    auto stream = ch_frb_io::intensity_network_stream::make(ini_params);
    return make_chime_network_stream(stream, 0, prescale);
}

#endif  // HAVE_CH_FRB_IO


// Json deserialization follows.
// FIXME: currently implemented for dummy_chime_network_stream, but not chime_network_stream.

static shared_ptr<wi_stream> dummy_chime_network_stream_from_json(const Json::Value &j)
{
    ssize_t nt_tot = ssize_t_from_json(j, "nt_tot");
    int nupfreq = int_from_json(j, "nupfreq");
    int nt_per_packet = int_from_json(j, "nt_per_packet");
    int fpga_counts_per_sample = int_from_json(j, "fpga_counts_per_sample");
    double pool_gb = double_from_json(j, "pool_gb");
    
    return make_dummy_chime_network_stream(nt_tot, nupfreq, nt_per_packet, fpga_counts_per_sample, pool_gb);
}

namespace {
    struct _init {
	_init() {
	    pipeline_object::register_json_deserializer("chime_dummy_network_stream", dummy_chime_network_stream_from_json);
	}
    } init;
}


}   // namespace rf_pipelines
