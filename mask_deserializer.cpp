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
shared_ptr<pipeline_object> make_mask_deserializer(const string &hdf5_filename)
{
    throw runtime_error("rf_pipelines::make_mask_deserializer() called, but rf_pipelines was compiled without HDF5");
}

struct mask_deserializer {
    static shared_ptr<pipeline_object> from_json(const Json::Value &j)
    {
	throw runtime_error("rf_pipelines: mask_deserializer appeared in json file, but rf_pipelines was compiled without HDF5");
    }
};

#else  // HAVE_HDF5

// Helper function which returns true if an HDF5 dataset is of type 'int8' or 'uint8'.
// FIXME something equivalent should be in sp_hdf5.
static bool is_int8(const H5::DataSet &ds)
{
    H5::DataType t = ds.getDataType();

    // For a list of all H5T_class_t values, see H5Tpublic.h
    H5T_class_t c = t.getClass();

    return (c == H5T_INTEGER) && (t.getSize() == 1);
}


static bool all_zeros(const float *arr, ssize_t nfreq, ssize_t nt, ssize_t stride)
{
    for (ssize_t ifreq = 0; ifreq < nfreq; ifreq++)
	for (ssize_t it = 0; it < nt; it++)
	    if (arr[ifreq*stride + it])
		return false;

    return true;
}


struct mask_deserializer : public chunked_pipeline_object
{
    // Initialized in constructor.
    const string filename;
    const rf_kernels::dequantizer dq;

    // Initialized in _bindc().
    // Note that _bindc() also initializes 'nt_chunk', which is defined in base class.
    shared_ptr<ring_buffer> rb_weights;
    ssize_t nfreq = 0;

    // Initialized in _allocate().
    uptr<uint8_t> tmp_buf;

    // Initialized in _start_pipeline().
    H5::DataSet input_dset;
    ssize_t file_i0 = 0;
    ssize_t file_i1 = 0;
    int verbosity = 0;


    mask_deserializer(const string &hdf5_filename) :
	chunked_pipeline_object("mask_deserializer", false),   // can_be_first = false
	filename(hdf5_filename),
	dq(1)   // nbits = 1
    { 
	this->name = "mask_deserializer(" + filename + ")";
    }


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

	int pipeline_fpga_counts_per_sample = int_from_json(json_attrs, "fpga_counts_per_sample");
	ssize_t pipeline_initial_fpga_count = ssize_t_from_json(json_attrs, "initial_fpga_count");

	H5::H5File f = sp_hdf5::hdf5_open(filename);

	int file_fpga_counts_per_sample = sp_hdf5::hdf5_read_attribute<int> (f, "fpga_counts_per_sample");
	ssize_t file_initial_fpga_count = sp_hdf5::hdf5_read_attribute<ssize_t> (f, "initial_fpga_count");
	ssize_t df = file_initial_fpga_count - pipeline_initial_fpga_count;

	if (pipeline_fpga_counts_per_sample != file_fpga_counts_per_sample)
	    _throw("different values of 'fpga_counts_per_sample' in pipeline and hdf5 file");
	if (df % pipeline_fpga_counts_per_sample != 0)
	    _throw("initial_fpga_counts are not divisible by fpga_counts_per_sample");

	this->input_dset = sp_hdf5::hdf5_open_dataset(f, "bitmask");

	vector<hsize_t> shape = sp_hdf5::hdf5_get_shape(input_dset);
	
	if (!is_int8(input_dset))
	    _throw("expected 'bitmask' hdf5 dataset to be of 8-bit integer type");
	if (shape.size() != 2)
	    _throw("expected 'bitmask' hdf5 dataset to be two-dimensional");
	if (shape[0] != size_t(nfreq))
	    _throw("number of frequency channels in hdf5 file does not match pipeline");

	this->file_i0 = df / pipeline_fpga_counts_per_sample;
	this->file_i1 = file_i0 + 8 * shape[1];
	this->verbosity = get_params().verbosity;

	if (verbosity >= 2) {
	    cout << "mask_deserializer: read " << filename
		 << ", initial_fpga_count=" << file_initial_fpga_count
		 << ", fpga_counts_per_sample=" << file_fpga_counts_per_sample
		 << ", nsamples=" << (file_i1 - file_i0) << endl;
	}
    }

    
    virtual bool _process_chunk(ssize_t pos) override
    {
	ring_buffer_subarray weights(rb_weights, pos, pos+nt_chunk, ring_buffer::ACCESS_RW);

	ssize_t j0 = 0;
	ssize_t j1 = nt_chunk;

	if (file_i0 > pos) {
	    j0 = min(file_i0 - pos, nt_chunk);
	    if (!all_zeros(weights.data, nfreq, j0, weights.stride))
		_throw("hdf5 file does not cover the entire time range of the acquisition");
	}

	if (file_i1 < pos + nt_chunk) {
	    j1 = max(file_i1 - pos, ssize_t(0));
	    if (!all_zeros(weights.data + j1, nfreq, nt_chunk - j1, weights.stride))
		_throw("hdf5 file does not cover the entire time range of the acquisition");
	}

	if (j0 == j1)
	    return true;

	rf_assert((0 <= j0) && (j0 <= j1) && (j1 <= nt_chunk));
	rf_assert((j1 - j0) % dq.kernel_size == 0);

	vector<hsize_t> chunk_offset = { 0, hsize_t((pos+j0-file_i0) / 8) };
	vector<hsize_t> chunk_shape = { hsize_t(nfreq), hsize_t((j1-j0) / 8) };
	
	sp_hdf5::hdf5_read_partial_dataset(input_dset, tmp_buf.get(), chunk_offset, chunk_shape);
	dq.apply_bitmask(nfreq, j1-j0, weights.data, weights.stride, tmp_buf.get(), (j1-j0)/8);

	return true;
    }


    virtual void _end_pipeline(Json::Value &json_output) override
    {
	this->input_dset = H5::DataSet();
	this->file_i0 = 0;
	this->file_i1 = 0;
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
	ret["class_name"] = "mask_deserializer";
	ret["hdf5_filename"] = filename;
	return ret;
    }


    static shared_ptr<mask_deserializer> from_json(const Json::Value &j)
    {
	const char *where = "rf_pipelines::mask_deserializer::from_json()";

	string hdf5_filename = string_from_json(j, "hdf5_filename", where);
	return make_shared<mask_deserializer> (hdf5_filename);
    }
};

// Externally callable factory function
shared_ptr<pipeline_object> make_mask_deserializer(const string &hdf5_filename)
{
    return make_shared<mask_deserializer> (hdf5_filename);
}
 
#endif  // HAVE_HDF5


namespace {
    struct _init {
	_init() {
	    pipeline_object::register_json_deserializer("mask_deserializer", mask_deserializer::from_json);
	}
    } init;
}

}  // namespace rf_pipelines
