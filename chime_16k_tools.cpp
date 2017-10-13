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


struct chime_16k_spike_mask : public chunked_pipeline_object
{
    static constexpr int nupfreq = 16;
    static constexpr int nfreq_c = 1024;
    static constexpr int nfreq_f = nfreq_c * nupfreq;

    std::shared_ptr<ring_buffer> rb_weights;


    chime_16k_spike_mask(ssize_t nt_chunk_=0) :
	chunked_pipeline_object("chime_16k_spike_mask", false)
    { 
	this->nt_chunk = nt_chunk_;
    }

    virtual void _bindc(ring_buffer_dict &rb_dict, Json::Value &json_attrs) override
    {
	this->rb_weights = this->get_buffer(rb_dict, "WEIGHTS");
	
	if (rb_weights->cdims.size() != 1)
	    _throw("expected weights array to be two-dimensional");
	if (rb_weights->cdims[0] != nfreq_f)
	    _throw("expected nfreq=" + to_string(nfreq_f));
    }

    virtual bool _process_chunk(ssize_t pos) override
    {
	ring_buffer_subarray weights(rb_weights, pos, pos+nt_chunk, ring_buffer::ACCESS_RW);

	for (int ifreq_c = 0; ifreq_c < nfreq_c; ifreq_c++) {
	    float *wrow = weights.data + (ifreq_c * nupfreq + 15) * weights.stride;
	    memset(wrow, 0, nt_chunk * sizeof(float));
	}

	return true;
    }	

    virtual void _unbindc() override
    {
	this->rb_weights.reset();
    }

    virtual Json::Value jsonize() const override
    {
	Json::Value ret;
	ret["class_name"] = "chime_16k_spike_mask";
	ret["nt_chunk"] = int(this->get_prebind_nt_chunk());
	return ret;
    }

    static shared_ptr<chime_16k_spike_mask> from_json(const Json::Value &j)
    {
	int nt_chunk = int_from_json(j, "nt_chunk");
	return make_shared<chime_16k_spike_mask> (nt_chunk);
    }
};


// Externally-visible factory function
shared_ptr<chunked_pipeline_object> make_chime_16k_spike_mask(ssize_t nt_chunk)
{
    return make_shared<chime_16k_spike_mask> (nt_chunk);
}


// -------------------------------------------------------------------------------------------------


struct chime_16k_derippler : public chunked_pipeline_object
{
    static constexpr int nupfreq = 16;
    static constexpr int nfreq_c = 1024;
    static constexpr int nfreq_f = nfreq_c * nupfreq;

    const double fudge_factor;
    double multiplier[nupfreq];

    std::shared_ptr<ring_buffer> rb_intensity;


    chime_16k_derippler(double fudge_factor_=1.0, ssize_t nt_chunk_=0) :
	chunked_pipeline_object("chime_16k_derippler", false),
	fudge_factor(fudge_factor_)
    { 
	constexpr float fl[nupfreq] = {
	    1.28255845837,
	    1.25644862084,
	    1.19869285334,
	    1.10585973916,
	    0.982096806832,
	    0.844561662108,
	    0.722474301172,
	    0.648123020509,
	    0.643478657879,
	    0.709964838879,
	    0.827771489632,
	    0.965153033863,
	    1.09186607691,
	    1.18902487817,
	    1.25106135586,
	    1.28086420647
	};

	this->nt_chunk = nt_chunk_;

	if ((fudge_factor < 0.0) || (fudge_factor > 2.0))
	    _throw("fudge_factor must be between 0 and 2");

	for (int i = 0; i < nupfreq; i++)
	    this->multiplier[i] = 1.0 / (1.0 + fudge_factor * (fl[i] - 1.0));
    }

    virtual void _bindc(ring_buffer_dict &rb_dict, Json::Value &json_attrs) override
    {
	this->rb_intensity = this->get_buffer(rb_dict, "INTENSITY");
	
	if (rb_intensity->cdims.size() != 1)
	    _throw("expected intensity array to be two-dimensional");
	if (rb_intensity->cdims[0] != nfreq_f)
	    _throw("expected nfreq=" + to_string(nfreq_f));
    }

    virtual bool _process_chunk(ssize_t pos) override
    {
	ring_buffer_subarray intensity(rb_intensity, pos, pos+nt_chunk, ring_buffer::ACCESS_RW);

	for (int ifreq_c = 0; ifreq_c < nfreq_c; ifreq_c++) {
	    for (int iupfreq = 0; iupfreq < nupfreq; iupfreq++) {
		float *irow = intensity.data + (ifreq_c*nupfreq + iupfreq) * intensity.stride;
		float t = multiplier[iupfreq];
		
		for (int i = 0; i < nt_chunk; i++)
		    irow[i] *= t;
	    }
	}

	return true;
    }	

    virtual void _unbindc() override
    {
	this->rb_intensity.reset();
    }

    virtual Json::Value jsonize() const override
    {
	Json::Value ret;
	ret["class_name"] = "chime_16k_derippler";
	ret["fudge_factor"] = this->fudge_factor;
	ret["nt_chunk"] = int(this->get_prebind_nt_chunk());
	return ret;
    }

    static shared_ptr<chime_16k_derippler> from_json(const Json::Value &j)
    {
	int nt_chunk = int_from_json(j, "nt_chunk");
	double fudge_factor = double_from_json(j, "fudge_factor");
	return make_shared<chime_16k_derippler> (fudge_factor, nt_chunk);
    }
};


// Externally-visible factory function
shared_ptr<chunked_pipeline_object> make_chime_16k_derippler(double fudge_factor, ssize_t nt_chunk)
{
    return make_shared<chime_16k_derippler> (fudge_factor, nt_chunk);
}


// -------------------------------------------------------------------------------------------------
//
// chime_16k_stripe_analyzer: for now, just outputs an hdf5 file.
//
// Not very optimized, but shouldn't be gratuitously slow..


#ifndef HAVE_HDF5

shared_ptr<wi_transform> make_chime_16k_stripe_analyzer(ssize_t Dt1, ssize_t Df2, ssize_t Dt2)
{
    throw runtime_error("rf_pipelines::make_chime_16k_stripe_analyzer() was called, but this rf_pipelines was compiled without HAVE_HDF5");
}

#else // HAVE_HDF5

struct chime_16k_stripe_analyzer : public wi_transform
{
    const ssize_t Dt1;
    const ssize_t Df2;
    const ssize_t Dt2;

    // Downsample (16384, Dt1*Dt2) -> (16384, Dt2)
    unique_ptr<rf_kernels::wi_downsampler> ds_kernel;
    uptr<float> ds_intensity;
    uptr<float> ds_weights;

    // Downsample (16384, Dt2) -> (1024, Dt2)
    unique_ptr<rf_kernels::wi_downsampler> ds2_kernel;
    uptr<float> ds2_intensity;
    uptr<float> ds2_weights;

    ssize_t nfreq_h5 = 0;          // = 1024/Df2
    uptr<float> h5_chunk;  // shape (16, nfreq_h5)

    string h5_fullname;
    H5::H5File h5_file;
    unique_ptr<sp_hdf5::hdf5_extendable_dataset<float>> h5_dset;

    // Temp buffer used to calculate median.
    vector<float> median_buf;


    chime_16k_stripe_analyzer(ssize_t Dt1_=16, ssize_t Df2_=16, ssize_t Dt2_=16) :
	wi_transform("chime_16k_stripe_analyzer"),
	Dt1(Dt1_),
	Df2(Df2_),
	Dt2(Dt2_)
    {
	if (Dt1 <= 0)
	    _throw("expected Dt1 > 0");
	if (Dt2 <= 0)
	    _throw("expected Dt2 > 0");
	if ((Df2 <= 0) || (1024 % Df2 != 0))
	    _throw("expected Df2 to be a positive divisor of 1024");

	this->nds = 1;
	this->nfreq = 16384;
	this->nt_chunk = Dt1 * Dt2;  // for convenience
	this->nfreq_h5 = xdiv(1024, Df2);
    }

    virtual void _allocate() override
    {
	this->ds_kernel = make_unique<rf_kernels::wi_downsampler> (1, Dt1);
	this->ds_intensity = make_uptr<float> (16384 * Dt2);
	this->ds_weights = make_uptr<float> (16384 * Dt2);

	this->ds2_kernel = make_unique<rf_kernels::wi_downsampler> (16, 1);
	this->ds2_intensity = make_uptr<float> (1024 * Dt2);
	this->ds2_weights = make_uptr<float> (1024 * Dt2);

	this->h5_chunk = make_uptr<float> (16 * nfreq_h5);
    }

    virtual void _start_pipeline(Json::Value &json_attrs) override
    {
	vector<hsize_t> chunk_shape = { 16, hsize_t(nfreq_h5), 1 };

	// check for filename collision
	this->h5_fullname = this->out_mp->add_file("stripe_analysis.h5");

	// open file
	this->h5_file = sp_hdf5::hdf5_open_trunc(h5_fullname);
	this->h5_dset = make_unique<sp_hdf5::hdf5_extendable_dataset<float>> (h5_file, "data", chunk_shape, 2);
    }

    virtual void _process_chunk(float *intensity, ssize_t istride, float *weights, ssize_t wstride, ssize_t pos) override
    {
	// Downsample (16384, Dt1*Dt2) -> (16384, Dt2)
	ds_kernel->downsample(16384, Dt2,               // nfreq_out, nt_out
			      ds_intensity.get(), Dt2,  // out_i, out_istride
			      ds_weights.get(), Dt2,    // out_w, out_wstride
			      intensity, istride,       // in_i, in_istride
			      weights, wstride);        // in_w, in_wstirde

	// Downsample (16384, Dt2) -> (1024, Dt2)
	ds2_kernel->downsample(1024, Dt2,                  // nfreq_out, nt_out
			       ds2_intensity.get(), Dt2,   // out_i, out_istride
			       ds2_weights.get(), Dt2,     // out_w, out_wstride
			       ds_intensity.get(), Dt2,    // in_i, in_istride
			       ds_weights.get(), Dt2);     // in_w, in_wstride

	// Divide 'ds' by 'ds2'.
	for (int ifreq_c = 0; ifreq_c < 1024; ifreq_c++) {
	    float *irow = ds_intensity.get() + ifreq_c * 16 * Dt2;  // shape (16,Dt2)
	    float *wrow = ds_weights.get() + ifreq_c * 16 * Dt2;    // shape (16,Dt2)
	    float *irow2 = ds2_intensity.get() + ifreq_c * Dt2;     // length Dt2

	    for (int it = 0; it < Dt2; it++) {
		if (irow2[it] > 0.0) {
		    for (int u = 0; u < 16; u++)
			irow[u*Dt2+it] /= irow2[it];
		}
		else {
		    for (int u = 0; u < 16; u++)
			wrow[u*Dt2+it] = 0.0;
		}
	    }
	}

	// Fill shape-(16,nfreq_h5) array, by computing medians in shape-(Df2,Dt2) blocks.
	memset(this->h5_chunk.get(), 0, 16 * nfreq_h5 * sizeof(float));

	for (int ifreq_h = 0; ifreq_h < nfreq_h5; ifreq_h++) {
	    for (int u = 0; u < 16; u++) {
		int ifreq_f = (ifreq_h * Df2 * 16) + u;

		// shape-(Df2,Dt2) strided arrays
		float *tmp_i = ds_intensity.get() + ifreq_f * Dt2;
		float *tmp_w = ds_weights.get() + ifreq_f * Dt2;
		int tstride = 16 * Dt2;

		median_buf.clear();

		for (int i = 0; i < Df2; i++) {
		    for (int j = 0; j < Dt2; j++) {
			if (tmp_w[i*tstride+j] > 0.0)
			    median_buf.push_back(tmp_i[i*tstride+j]);
		    }
		}

		if ((int)median_buf.size() > (Df2*Dt2)/4)
		    h5_chunk[u*nfreq_h5 + ifreq_h] = median(median_buf);
	    }
	}
	
	vector<hsize_t> chunk_shape = { 16, hsize_t(nfreq_h5), 1 };
	this->h5_dset->write(h5_chunk.get(), chunk_shape);
    }

    virtual void _end_pipeline(Json::Value &json_output) override
    {
	this->h5_dset.reset();
	this->h5_file.close();

	// FIXME if (verbosity >= 2) ...
	cout << "chime_16k_stripe_analyzer: wrote " << h5_fullname << endl;
    }

    virtual void _deallocate() override
    {
	this->ds_kernel.reset();
	this->ds_intensity.reset();
	this->ds_weights.reset();
	this->ds2_kernel.reset();
	this->ds2_intensity.reset();
	this->ds2_weights.reset();
	this->h5_chunk.reset();
    }	

    virtual Json::Value jsonize() const override
    {
	Json::Value ret;
	ret["class_name"] = "chime_16k_stripe_analyzer";
	ret["Dt1"] = Json::Int64(Dt1);
	ret["Df2"] = Json::Int64(Df2);
	ret["Dt2"] = Json::Int64(Dt2);
	return ret;
    }

    static shared_ptr<pipeline_object> from_json(const Json::Value &j)
    {
	ssize_t Dt1 = int_from_json(j, "Dt1");
	ssize_t Df2 = int_from_json(j, "Df2");
	ssize_t Dt2 = int_from_json(j, "Dt2");

	return make_shared<chime_16k_stripe_analyzer> (Dt1, Df2, Dt2);
    }
};

shared_ptr<wi_transform> make_chime_16k_stripe_analyzer(ssize_t Dt1, ssize_t Df2, ssize_t Dt2)
{
    return make_shared<chime_16k_stripe_analyzer> (Dt1, Df2, Dt2);
}

#endif // HAVE_HDF5


// -------------------------------------------------------------------------------------------------


namespace {
    struct _init {
	_init() {
	    pipeline_object::register_json_deserializer("chime_16k_derippler", chime_16k_derippler::from_json);
	    pipeline_object::register_json_deserializer("chime_16k_spike_mask", chime_16k_spike_mask::from_json);
	    pipeline_object::register_json_deserializer("chime_16k_stripe_analyzer", chime_16k_stripe_analyzer::from_json);
	}
    } init;
}


}  // namespace rf_pipelines
