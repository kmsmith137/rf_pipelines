#include "rf_pipelines_internals.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
};  // pacify emacs c-mode!
#endif


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


namespace {
    struct _init {
	_init() {
	    pipeline_object::register_json_constructor("chime_16k_destriper", chime_16k_destriper::from_json);
	}
    } init;
}


}  // namespace rf_pipelines
