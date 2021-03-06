#include <cassert>
#include "rf_pipelines_internals.hpp"
#include "rf_pipelines_inventory.hpp"

#ifdef HAVE_CH_FRB_IO
#include <ch_frb_io.hpp>
#else
namespace ch_frb_io { class assembled_chunk; }
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


void mask_counter_transform::set_runtime_attrs(const runtime_attrs &a)
{
    if (this->state != UNBOUND)
	throw runtime_error("mask_counter_transform::set_runtime_attrs() called after bind()");
    if (a.ringbuf_nhistory < 0)
	throw runtime_error("mask_counter_transform::set_runtime_attrs(): ringbuf_nhistory was negative");
    if (a.chime_stream && (a.chime_beam_id < 0))
	throw runtime_error("mask_counter_transform::set_runtime_attrs(): chime_stream was specified, but chime_beam_id was uninitialized or negative");

#ifndef HAVE_CH_FRB_IO
    if (a.chime_stream)
	throw runtime_error("mask_counter_transform::set_runtime_attrs(): chime_stream was specified, but rf_pipelines was compiled without ch_frb_io!");
#endif

    this->attrs = a;
    this->ringbuf.reset();

    if (a.ringbuf_nhistory > 0)
	this->ringbuf = make_shared<mask_measurements_ringbuf> (a.ringbuf_nhistory);
}


// virtual override
void mask_counter_transform::_bind_transform(Json::Value &json_attrs)
{
#ifdef HAVE_CH_FRB_IO
    if (attrs.chime_stream) {
	// This is the earliest place we can put this assert (since this->nfreq initialized in bind())
	if (attrs.chime_stream->ini_params.nrfifreq != this->nfreq)
	    throw runtime_error("mask_counter: value of 'nrfifreq' in chime_intensity_stream does not match the frequency resolution in the RFI transform chain");

	// Should be redundant with asserts elsewhere in rf_pipelines, but just being paranoid!
	rf_assert(this->nds == 1);
    }
#endif

    // Check for other mask_counters with duplicate "where" names.
    string keyname = "mask_counter_name_list";
    if (!json_attrs.isMember(keyname))
        // Add new list.
        json_attrs[keyname] = Json::Value(Json::arrayValue);
    Json::Value& names = json_attrs[keyname];
    // assert that it's an Array
    if (names.type() != Json::arrayValue)
        throw runtime_error("mask_counter: expected JSON attribute '" + keyname + "' to be an array (of 'where' entries)");
    // check for duplicates
    for (uint i=0; i<names.size(); i++) {
        if (!names[i].isString())
            throw runtime_error("mask_counter: expected JSON attribute '" + keyname + "' to contain strings");
        string othername = names[i].asString();
        if (othername == where)
            throw runtime_error("mask_counter: 'where' entry = \"" + where + "\" is duplicated in the pipeline -- it must be unique!");
    }
    names.append(Json::Value(where));

}


// virtual override
void mask_counter_transform::_start_pipeline(Json::Value &json_attrs)
{
#ifdef HAVE_CH_FRB_IO
    if (attrs.chime_stream) {
	this->chime_initial_fpga_count = uint64_t_from_json(json_attrs, "initial_fpga_count");
	this->chime_fpga_counts_per_sample = int_from_json(json_attrs, "fpga_counts_per_sample");
	this->chime_fpga_counts_initialized = true;
	
	if (attrs.chime_stream->ini_params.fpga_counts_per_sample != chime_fpga_counts_per_sample)
	    throw runtime_error("mask_counter: value of 'fpga_counts_per_sample' in chime_intensity_stream does not match the value in _start_pipeline()");
    }
#endif
}


// virtual override
void mask_counter_transform::_process_chunk(float *intensity, ssize_t istride, float *weights, ssize_t wstride, ssize_t pos)
{
    mask_measurements meas;

    rf_kernels::mask_counter_data d;
    d.nfreq = nfreq;
    d.nt_chunk = nt_chunk / nds;      // not 'nt_chunk'
    d.in = weights;                   // not 'intensity'
    d.istride = wstride;              // not 'istride'

    if (ringbuf) {
	meas = mask_measurements(pos, d.nfreq, d.nt_chunk);
	d.out_fcounts = meas.freqs_unmasked.get();
    }

#ifdef HAVE_CH_FRB_IO
    // Declared outside the if-statement below, so that we hold the shared_ptr<> reference while the kernel is being called.
    shared_ptr<ch_frb_io::assembled_chunk> chunk;

    if (attrs.chime_stream) {
	if (!chime_fpga_counts_initialized)
	    throw runtime_error("rf_pipelines::chime_mask_counter internal error: fpga count fields were not initialized as expected");
    
	// The 'pos' argument is the current pipeline position in units of time samples -- convert to FPGA counts
	uint64_t fpga_counts = pos * this->chime_fpga_counts_per_sample + this->chime_initial_fpga_count;

	// The last argument in find_assembled_chunk() is 'toplevel'.
	chunk = attrs.chime_stream->find_assembled_chunk(attrs.chime_beam_id, fpga_counts, true);

	// Reminder: find_assembled_chunk() returns an empty pointer iff stream has ended, and chunk is requested past end-of-stream.
	if (chunk) {
	    if (!chunk->rfi_mask)
		throw runtime_error("mask_counter: chime_intensity_stream::find_assembled_chunk() returned chunk with no RFI mask");
	    if (chunk->nrfifreq != this->nfreq)
		throw runtime_error("mask_counter: chime_intensity_stream::find_assembled_chunk() returned chunk with mismatched 'nrfifreq'");
	    if (chunk->has_rfi_mask)
		throw runtime_error("mask_counter: chime_intensity_stream::find_assembled_chunk() returned chunk with has_rfi_mask=true");
	    if (chunk->binning != 1)
		throw runtime_error("mask_counter: chime_intensity_stream::find_assembled_chunk() returned chunk with binning != 1");

	    d.out_bitmask = chunk->rfi_mask;
	    d.out_bmstride = d.nt_chunk / 8;   // contiguous
	}
    }
#endif

    // Run mask-counting kernel.
    int nunmasked = d.mask_count();
    this->nunmasked_tot += nunmasked;

    if (ringbuf) {
	meas.nsamples_unmasked = nunmasked;
	ringbuf->add(meas);
    }

#ifdef HAVE_CH_FRB_IO
    if (chunk) {
	chunk->has_rfi_mask = true;

	// Notify stream's output_devices that a chunk has had its rfi_mask filled in.
	for (auto od : attrs.chime_stream->ini_params.output_devices)
	    od->filled_rfi_mask(chunk);
    }
#endif
}


// virtual override
void mask_counter_transform::_end_pipeline(Json::Value &json_attrs)
{
    string k = "nunmasked_samples_" + this->where;
    json_attrs[k] = Json::Int64(this->nunmasked_tot);
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
