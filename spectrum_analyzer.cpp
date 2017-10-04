#include <rf_kernels/downsample.hpp>
#include "rf_pipelines_internals.hpp"

#ifdef HAVE_HDF5
#include <sp_hdf5.hpp>
#endif


using namespace std;

namespace rf_pipelines {
#if 0
};  // pacify emacs c-mode!
#endif


// -------------------------------------------------------------------------------------------------
//
// spectrum_analyzer: for now, just outputs an hdf5 file.
//
// Not very optimized, but shouldn't be gratuitously slow..


#ifndef HAVE_HDF5

shared_ptr<wi_transform> make_spectrum_analyzer(ssize_t Dt1, ssize_t Dt2)
{
    throw runtime_error("rf_pipelines::make_spectrum_analyzer() was called, but this rf_pipelines was compiled without HAVE_HDF5");
}

#else // HAVE_HDF5

struct spectrum_analyzer : public wi_transform
{
    const ssize_t Dt1;
    const ssize_t Dt2;

    // Downsample (nfreq, Dt1*Dt2) -> (nfreq, Dt2)
    unique_ptr<rf_kernels::wi_downsampler> ds_kernel;
    uptr<float> ds_intensity;
    uptr<float> ds_weights;

    string h5_fullname;
    H5::H5File h5_file;
    unique_ptr<sp_hdf5::hdf5_extendable_dataset<float>> h5_dset;
    uptr<float> h5_chunk;   // 1d array of length nfreq

    // Temp buffer used to calculate median.
    vector<float> median_buf;


    spectrum_analyzer(ssize_t Dt1_=16, ssize_t Dt2_=16) :
	wi_transform("spectrum_analyzer"),
	Dt1(Dt1_),
	Dt2(Dt2_)
    {
	if (Dt1 <= 0)
	    _throw("expected Dt1 > 0");
	if (Dt2 <= 0)
	    _throw("expected Dt2 > 0");

	this->nds = 1;
	this->nt_chunk = Dt1 * Dt2;  // for convenience
    }

    virtual void _allocate() override
    {
	rf_assert(nfreq > 0);

	this->ds_kernel = make_unique<rf_kernels::wi_downsampler> (1, Dt1);
	this->ds_intensity = make_uptr<float> (nfreq * Dt2);
	this->ds_weights = make_uptr<float> (nfreq * Dt2);
	this->h5_chunk = make_uptr<float> (nfreq);
    }

    virtual void _start_pipeline(Json::Value &json_attrs) override
    {
	vector<hsize_t> chunk_shape = { hsize_t(nfreq), 1 };

	// check for filename collision
	this->h5_fullname = this->out_mp->add_file("spectrum.h5");

	// open file
	this->h5_file = sp_hdf5::hdf5_open_trunc(h5_fullname);
	this->h5_dset = make_unique<sp_hdf5::hdf5_extendable_dataset<float>> (h5_file, "data", chunk_shape, 1);
    }

    virtual void _process_chunk(float *intensity, ssize_t istride, float *weights, ssize_t wstride, ssize_t pos) override
    {
	// Downsample (nfreq, Dt1*Dt2) -> (nfreq, Dt2)
	ds_kernel->downsample(nfreq, Dt2,               // nfreq_out, nt_out
			      ds_intensity.get(), Dt2,  // out_i, out_istride
			      ds_weights.get(), Dt2,    // out_w, out_wstride
			      intensity, istride,       // in_i, in_istride
			      weights, wstride);        // in_w, in_wstirde

	// Fill length-nfreq array, by computing medians in length-Dt2 blocks.
	memset(this->h5_chunk.get(), 0, nfreq * sizeof(float));

	for (int ifreq = 0; ifreq < nfreq; ifreq++) {
	    float *irow = ds_intensity.get() + ifreq * Dt2;
	    float *wrow = ds_weights.get() + ifreq * Dt2;

	    median_buf.clear();
	    
	    for (int it = 0; it < Dt2; it++)
		if (wrow[it] > 0.0)
		    median_buf.push_back(irow[it]);

	    if ((int)median_buf.size() > Dt2/4)
		h5_chunk[ifreq] = median(median_buf);
	}
	
	vector<hsize_t> chunk_shape = { hsize_t(nfreq), 1 };
	this->h5_dset->write(h5_chunk.get(), chunk_shape);
    }

    virtual void _end_pipeline(Json::Value &json_output) override
    {
	this->h5_dset.reset();
	this->h5_file.close();

	// FIXME if (verbosity >= 2) ...
	cout << "spectrum_analyzer: wrote " << h5_fullname << endl;
    }

    virtual void _deallocate() override
    {
	this->ds_kernel.reset();
	this->ds_intensity.reset();
	this->ds_weights.reset();
	this->h5_chunk.reset();
    }	

    virtual Json::Value jsonize() const override
    {
	Json::Value ret;
	ret["class_name"] = "spectrum_analyzer";
	ret["Dt1"] = Json::Int64(Dt1);
	ret["Dt2"] = Json::Int64(Dt2);
	return ret;
    }

    static shared_ptr<pipeline_object> from_json(const Json::Value &j)
    {
	ssize_t Dt1 = int_from_json(j, "Dt1");
	ssize_t Dt2 = int_from_json(j, "Dt2");

	return make_shared<spectrum_analyzer> (Dt1, Dt2);
    }
};

shared_ptr<wi_transform> make_spectrum_analyzer(ssize_t Dt1, ssize_t Dt2)
{
    return make_shared<spectrum_analyzer> (Dt1, Dt2);
}

#endif // HAVE_HDF5


// -------------------------------------------------------------------------------------------------


namespace {
    struct _init {
	_init() {
	    pipeline_object::register_json_constructor("spectrum_analyzer", spectrum_analyzer::from_json);
	}
    } init;
}


}  // namespace rf_pipelines
