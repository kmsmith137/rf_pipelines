#include "rf_pipelines_internals.hpp"
#include "rf_pipelines_inventory.hpp"

#ifdef HAVE_CH_FRB_IO
#include <ch_frb_io.hpp>
#endif

using namespace std;

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif


mask_counter_transform::mask_counter_transform(int nt_chunk_, string where_) :
    wi_transform("mask_counter"),
    where(where_)
{	
    stringstream ss;
    ss << "mask_counter(nt_chunk=" << nt_chunk_ << ", where=" << where << ")";
    this->name = ss.str();
    this->nt_chunk = nt_chunk_;
    this->nds = 0;  // allows us to run in a wi_sub_pipeline

    if (nt_chunk == 0)
        throw runtime_error("rf_pipelines::mask_counter: nt_chunk must be specified");
}


void mask_counter_transform::set_callbacks(const std::shared_ptr<mask_counter_callbacks> &callbacks_)
{
    if (!callbacks_)
	throw runtime_error("mask_counter_transform::set_callbacks() called with empty pointer");
    if (this->callbacks)
	throw runtime_error("mask_counter_transform::set_callbacks() was called twice");
    if (this->state != UNBOUND)
	throw runtime_error("mask_counter_transform::set_callbacks() called after bind()");

    this->callbacks = callbacks_;
}


// virtual override
void mask_counter_transform::_bind_transform(Json::Value &json_attrs)
{
    if (callbacks)
	callbacks->_bind_transform(*this, json_attrs);
}


// virtual override
void mask_counter_transform::_start_pipeline(Json::Value &json_attrs)
{
    if (callbacks)
	callbacks->_start_pipeline(*this, json_attrs);
}


// virtual override
void mask_counter_transform::_process_chunk(float *intensity, ssize_t istride, float *weights, ssize_t wstride, ssize_t pos)
{
    rf_kernels::mask_counter_data d;
    d.nfreq = nfreq;
    d.nt_chunk = nt_chunk / nds;      // not 'nt_chunk'
    d.in = weights;                   // not 'intensity'
    d.istride = wstride;              // not 'istride'

    // Run mask-counting kernel.
    int nunmasked = callbacks ? callbacks->_run_kernel(*this,d,pos) : d.mask_count();
    int ntot = d.nfreq * d.nt_chunk;

    // cout << "mask_counter " << where << ", pos " << pos 
    // << ": N samples masked: " << (ntot - nunmasked)
    // << "/" << ntot << endl;

    this->nmasked_tot += (ntot - nunmasked);
}


// virtual override
void mask_counter_transform::_end_pipeline(Json::Value &json_attrs)
{
    string k = "nmasked_samples_" + this->where;
    json_attrs[k] = Json::Int64(this->nmasked_tot);
}


Json::Value mask_counter_transform::jsonize() const
{
    Json::Value ret;

    ret["class_name"] = "mask_counter";
    ret["nt_chunk"] = int(this->get_prebind_nt_chunk());
    ret["where"] = where;
    return ret;
}


shared_ptr<mask_counter_transform>
mask_counter_transform::from_json(const Json::Value &j)
{
    ssize_t nt_chunk = ssize_t_from_json(j, "nt_chunk");
    string where = string_from_json(j, "where");
    return make_shared<mask_counter_transform> (nt_chunk, where);
}


// -------------------------------------------------------------------------------------------------
//
// chime_mask_counter_callbacks


#ifdef HAVE_CH_FRB_IO

chime_mask_counter_callbacks::chime_mask_counter_callbacks(const shared_ptr<ch_frb_io::intensity_network_stream> &stream_, int beam_id_, int ringbuf_nhistory) :
    stream(stream_),
    ringbuf(make_shared<mask_measurements_ringbuf>(ringbuf_nhistory)),
    beam_id(beam_id_)
{ }


// virtual override
void chime_mask_counter_callbacks::_bind_transform(mask_counter_transform &self, Json::Value &json_attrs)
{
    // _bind_transform() is the earliest place we can put this assert (since self.nfreq initialized in bind())
    if (stream && (stream->ini_params.nrfifreq != self.nfreq))
	throw runtime_error("chime_mask_counter: value of 'nrfifreq' does not match the frequency resolution in the RFI transform chain");

    // Should be redundant with asserts elsewhere in rf_pipelines, but just being paranoid!
    rf_assert(self.nds == 1);
}


// virtual override
void chime_mask_counter_callbacks::_start_pipeline(mask_counter_transform &self, Json::Value &j)
{
    this->initial_fpga_count = uint64_t_from_json(j, "initial_fpga_count");
    this->fpga_counts_per_sample = int_from_json(j, "fpga_counts_per_sample");
    this->fpga_counts_initialized = true;

    if (stream && (stream->ini_params.fpga_counts_per_sample != fpga_counts_per_sample))
	throw runtime_error("chime_mask_counter: mismatched values of 'fpga_counts_per_sample' in _start_pipeline()");
}


// virtual override
int chime_mask_counter_callbacks::_run_kernel(mask_counter_transform &self, rf_kernels::mask_counter_data &d, ssize_t pos)
{
    if (!fpga_counts_initialized)
	throw runtime_error("rf_pipelines::chime_mask_counter internal error: fpga count fields were not initialized as expected");
    
    // The 'pos' argument is the current pipeline position in units of time samples (not FPGA counts)
    uint64_t fpga_counts = pos * this->fpga_counts_per_sample + this->initial_fpga_count;
    //cout << "chime_mask_counter: finding chunk for pos " << pos << " (fpga counts " << fpga_counts << ")" << endl;

    mask_measurements meas(pos, d.nfreq, d.nt_chunk);

    // Declared outside the if-statement below, so that we hold the shared_ptr<> reference while the kernel is being called.
    shared_ptr<ch_frb_io::assembled_chunk> chunk;

    if (stream) {
	// The last argument in find_assembled_chunk() is 'toplevel'.
	chunk = stream->find_assembled_chunk(beam_id, fpga_counts, true);

	// These should all be redundant with asserts in ch_frb_io, but a little paranoia never hurts.
	if (!chunk)
	    throw runtime_error("chime_mask_counter: find_assembled_chunk() returned empty pointer");
	if (!chunk->rfi_mask)
	    throw runtime_error("chime_mask_counter: find_assembled_chunk() returned chunk with no RFI mask");
	if (chunk->nrfifreq != self.nfreq)
	    throw runtime_error("chime_mask_counter: find_assembled_chunk() returned chunk with mismatched 'nrfifreq'");
	if (chunk->has_rfi_mask)
	    throw runtime_error("chime_mask_counter: find_assembled_chunk() returned chunk with has_rfi_mask=true");
	if (chunk->binning != 1)
	    throw runtime_error("chime_mask_counter: find_assembled_chunk() returned chunk with binning != 1");
    }

    d.out_bitmask = stream ? chunk->rfi_mask : nullptr;
    d.out_fcounts = meas.freqs_unmasked.get();
    d.out_bmstride = d.nt_chunk / 8;   // contiguous

    // Run mask-counting kernel.
    meas.nsamples_unmasked = d.mask_count();

    ringbuf->add(meas);

    if (stream) {
	chunk->has_rfi_mask = true;

	// Notify stream's output_devices that a chunk has had its rfi_mask filled in.
	for (auto od : stream->ini_params.output_devices)
	    od->filled_rfi_mask(chunk);
    }

    return meas.nsamples_unmasked;
}

#else  // HAVE_CH_FRB_IO

chime_mask_counter_callbacks::chime_mask_counter_callbacks(const shared_ptr<ch_frb_io::intensity_network_stream> &stream_, int beam_id_, int ringbuf_nhistory)
{
    throw runtime_error("rf_pipelines:chime_mask_counter_callbacks was constructed, but rf_pipelines was compiled without ch_frb_io");
}

#endif  // HAVE_CH_FRB_IO


// -------------------------------------------------------------------------------------------------


namespace {
    struct _init {
        _init() {
            pipeline_object::register_json_deserializer("mask_counter", mask_counter_transform::from_json);
        }
    } init;
}


// Externally callable
shared_ptr<wi_transform> make_mask_counter(int nt_chunk, string where)
{
    return make_shared<mask_counter_transform> (nt_chunk, where);
}


}  // namespace rf_pipelines
