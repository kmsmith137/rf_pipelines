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


struct chime_16k_destriper : public chunked_pipeline_object
{
    static constexpr int nupfreq = 16;
    static constexpr int nfreq_c = 1024;
    static constexpr int nfreq_f = nfreq_c * nupfreq;

    std::shared_ptr<ring_buffer> rb_intensity;


    chime_16k_destriper(ssize_t nt_chunk_=0) :
	chunked_pipeline_object("chime_16k_destriper", false, nt_chunk_)
    { }

    virtual void _bind_chunked(ring_buffer_dict &rb_dict, Json::Value &json_attrs) override
    {
	this->rb_intensity = this->get_buffer(rb_dict, "INTENSITY");
	
	if (rb_intensity->cdims.size() != 1)
	    _throw("expected intensity array to be two-dimensional");
	if (rb_intensity->cdims[0] != nfreq_f)
	    _throw("expected nfreq=" + to_string(nfreq_f));
    }

    virtual bool _process_chunk(ssize_t pos) override
    {
	constexpr float finv[nupfreq] = {
	    0.779691555944,
	    0.795894064757,
	    0.834242063942,
	    0.904273810309,
	    1.01822956051,
	    1.1840461684,
	    1.38413227762,
	    1.54291695921,
	    1.55405309524,
	    1.40852045797,
	    1.20806286823,
	    1.03610512003,
	    0.915863237391,
	    0.841025295903,
	    0.799321308517,
	    0.78072288612
	};

	float *intensity = rb_intensity->get(pos, pos+nt_chunk, ring_buffer::ACCESS_RW);
	ssize_t istride = rb_intensity->get_stride();

	for (int ifreq_c = 0; ifreq_c < nfreq_c; ifreq_c++) {
	    for (int iupfreq = 0; iupfreq < nupfreq; iupfreq++) {
		float t = finv[iupfreq];
		float *p = intensity + (ifreq_c*nupfreq + iupfreq) * istride;
		
		for (int i = 0; i < nt_chunk; i++)
		    p[i] *= t;
	    }
	}

	rb_intensity->put(intensity, pos, pos+nt_chunk, ring_buffer::ACCESS_RW);
	return true;
    }	

    virtual Json::Value jsonize() const override
    {
	Json::Value ret;
	ret["class_name"] = "chime_16k_destriper";
	ret["nt_chunk"] = int(this->get_orig_nt_chunk());
	return ret;
    }

    static shared_ptr<chime_16k_destriper> from_json(const Json::Value &j)
    {
	int nt_chunk = int_from_json(j, "nt_chunk");
	return make_shared<chime_16k_destriper> (nt_chunk);
    }
};


// Externally-visible factory function
shared_ptr<chunked_pipeline_object> make_chime_16k_destriper(ssize_t nt_chunk)
{
    return make_shared<chime_16k_destriper> (nt_chunk);
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

    string h5_fullpath;
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

    virtual void _bind_transform(Json::Value &json_attrs) override
    {
	// check for filename collision
	this->h5_fullpath = this->out_mp->add_file("stripe_analysis.h5");
    }

    virtual void _start_pipeline(Json::Value &json_attrs) override
    {
	vector<hsize_t> chunk_shape = { 16, hsize_t(nfreq_h5), 1 };

	// open file
	this->h5_file = sp_hdf5::hdf5_open_trunc(h5_fullpath);
	this->h5_dset = make_unique<sp_hdf5::hdf5_extendable_dataset<float>> (h5_file, "data", chunk_shape, 2);
    }

    virtual void _process_chunk(float *intensity, ssize_t istride, float *weights, ssize_t wstride, ssize_t pos) override
    {
	memset(this->h5_chunk.get(), 0, 16 * nfreq_h5 * sizeof(float));

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
	    float *irow = ds_intensity.get() + ifreq_c * 16 * Dt2;
	    float *wrow = ds_weights.get() + ifreq_c * 16 * Dt2;
	    float *irow2 = ds2_intensity.get() + ifreq_c * Dt2;

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
	
	// Compute medians in shape-(Df2,Dt2) blocks.
	for (int ifreq_h = 0; ifreq_h < nfreq_h5; ifreq_h++) {
	    for (int u = 0; u < 16; u++) {
		int ifreq_f = (ifreq_h * Df2 * 16) + u;

		// shape-(Df2,Dt2) strided arrays
		float *tmp_i = ds_intensity.get() + ifreq_f * Dt2;
		float *tmp_w = ds_weights.get() + ifreq_f * Dt2;
		int tmp_stride = 16 * Dt2;

		median_buf.clear();

		for (int i = 0; i < Df2; i++) {
		    for (int j = 0; j < Dt2; j++) {
			if (tmp_w[i*tmp_stride+j] > 0.0)
			    median_buf.push_back(tmp_i[i*tmp_stride+j]);
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
	    pipeline_object::register_json_constructor("chime_16k_destriper", chime_16k_destriper::from_json);
	    pipeline_object::register_json_constructor("chime_16k_stripe_analyzer", chime_16k_stripe_analyzer::from_json);
	}
    } init;
}


}  // namespace rf_pipelines
