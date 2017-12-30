// FIXME: currently the mask_(de)serializer assumes that the pipeline defines its timestamps with
// attributes named 'fpga_counts_per_sample' and 'initial_fpga_count', as is the case for all CHIME
// streams.  Need to settle on a systematic general convention for timestamps!

#include "rf_pipelines_internals.hpp"

#ifdef HAVE_HDF5
#include <sp_hdf5.hpp>
#include <rf_kernels/quantize.hpp>
#endif

using namespace std;

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif


#ifndef HAVE_HDF5

// Externally callable factory function
shared_ptr<pipeline_object> make_mask_serializer(const string &hdf5_filename)
{
    throw runtime_error("rf_pipelines::make_mask_serializer() called, but rf_pipelines was compiled without HDF5");
}

struct mask_serializer {
    static shared_ptr<pipeline_object> from_json(const Json::Value &j)
    {
	throw runtime_error("rf_pipelines: mask_serializer appeared in json file, but rf_pipelines was compiled without HDF5");
    }
};

#else // HAVE_HDF5

struct mask_serializer : public chunked_pipeline_object
{
    // Initialized in constructor.
    const string filename;
    const rf_kernels::quantizer q;

    // Initialized in _bindc().
    // Note that _bindc() also initializes 'nt_chunk', which is defined in base class.
    shared_ptr<ring_buffer> rb_weights;
    ssize_t nfreq = 0;

    // Initialized in _allocate().
    uptr<uint8_t> tmp_buf;

    // Initialized in _start_pipeline().
    unique_ptr<sp_hdf5::hdf5_extendable_dataset<uint8_t>> output_dset;
    ssize_t initial_fpga_count = 0;
    ssize_t fpga_counts_per_sample = 0;
    ssize_t samples_written = 0;
    int verbosity = 0;
    

    mask_serializer(const string &hdf5_filename) :
	chunked_pipeline_object("mask_serializer", false),   // can_be_first = false
	filename(hdf5_filename),
	q(1)   // nbits = 1
    { }

    virtual void _bindc(ring_buffer_dict &rb_dict, Json::Value &json_attrs) override
    {
	this->rb_weights = this->get_buffer(rb_dict, "WEIGHTS");

	if (rb_weights->cdims.size() != 1)
	    _throw("expected weights array to be one-dimensional");
	if (rb_weights->nds != 1)
	    _throw("expected downsampling factor to be 1");

	this->nfreq = rb_weights->cdims[0];
	this->nt_chunk = 1024;  // FIXME hardcoded for now, should this be more flexible?
    }

    virtual void _allocate() override
    {
	this->tmp_buf = make_uptr<uint8_t> (nfreq * (nt_chunk/8));
    }

    virtual void _start_pipeline(Json::Value &json_attrs) override
    {
	if (!json_attrs.isMember("initial_fpga_count") || !json_attrs.isMember("fpga_counts_per_sample"))
	    _throw("currently, we assume that 'initial_fpga_count' and 'fpga_counts_per_sample' are defined by stream");

	this->initial_fpga_count = ssize_t_from_json(json_attrs, "initial_fpga_count");
	this->fpga_counts_per_sample = int_from_json(json_attrs, "fpga_counts_per_sample");
	this->verbosity = get_params().verbosity;
	this->samples_written = 0;

	H5::H5File f = sp_hdf5::hdf5_open_trunc(filename);
	
	sp_hdf5::hdf5_write_attribute(f, "initial_fpga_count", initial_fpga_count);
	sp_hdf5::hdf5_write_attribute(f, "fpga_counts_per_sample", fpga_counts_per_sample);
	sp_hdf5::hdf5_write_attribute(f, "nfreq", nfreq);

	// We use bitshuffle if available, and fall back to uncompressed.
	vector<string> compression = { "bitshuffle", "none" };
	vector<hsize_t> chunk_shape = { hsize_t(nfreq), hsize_t(nt_chunk/8) };
	this->output_dset = make_unique<sp_hdf5::hdf5_extendable_dataset<uint8_t>> (f, "bitmask", chunk_shape, 1, compression);
    }
    
    virtual bool _process_chunk(ssize_t pos) override
    {
	vector<hsize_t> chunk_shape = { hsize_t(nfreq), hsize_t(nt_chunk/8) };
	ring_buffer_subarray weights(rb_weights, pos, pos+nt_chunk, ring_buffer::ACCESS_READ);

	this->q.quantize(nfreq, nt_chunk, tmp_buf.get(), nt_chunk/8, weights.data, weights.stride);
	this->output_dset->write(tmp_buf.get(), chunk_shape);
	this->samples_written += nt_chunk;

	return true;
    }

    virtual void _end_pipeline(Json::Value &json_output) override
    {
	this->output_dset.reset();
	
	if (verbosity >= 2) {
	    cout << "mask_serializer: wrote " << filename 
		 << ", initial_fpga_count=" << initial_fpga_count
		 << ", fpga_counts_per_sample=" << fpga_counts_per_sample
		 << ", samples_written=" << samples_written << endl;
	}
    }

    virtual void _deallocate() override
    {
	this->tmp_buf.reset();
    }

    virtual void _unbindc() override
    {
	this->rb_weights.reset();
    }

    virtual Json::Value jsonize() const override
    {
	Json::Value ret;
	ret["class_name"] = "mask_serializer";
	ret["hdf5_filename"] = filename;
	return ret;
    }

    static shared_ptr<mask_serializer> from_json(const Json::Value &j)
    {
	const char *where = "rf_pipelines::mask_serializer::from_json()";
	
	string hdf5_filename = string_from_json(j, "hdf5_filename", where);
	return make_shared<mask_serializer> (hdf5_filename);
    }
};

// Externally callable factory function
shared_ptr<pipeline_object> make_mask_serializer(const string &hdf5_filename)
{
    return make_shared<mask_serializer> (hdf5_filename);
}

#endif  // HAVE_HDF5


namespace {
    struct _init {
	_init() {
	    pipeline_object::register_json_deserializer("mask_serializer", mask_serializer::from_json);
	}
    } init;
}

}  // namespace rf_pipelines
