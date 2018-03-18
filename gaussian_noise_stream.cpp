#include <random>
#include "rf_pipelines_internals.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
};  // pacify emacs c-mode!
#endif


class gaussian_noise_stream : public wi_stream
{
protected:
    const ssize_t nt_tot;
    const double freq_lo_MHz;
    const double freq_hi_MHz;
    const double dt_sample;
    const double sample_rms;
    const bool randomize_weights;

    std::mt19937 rng;
    std::normal_distribution<float> gdist;

public:
    gaussian_noise_stream(ssize_t nfreq_, ssize_t nt_tot_, double freq_lo_MHz_, double freq_hi_MHz_, double dt_sample_, double sample_rms_, ssize_t nt_chunk_, bool randomize_weights_) :
	wi_stream("gaussian_noise_stream"),
	nt_tot(nt_tot_),
	freq_lo_MHz(freq_lo_MHz_),
	freq_hi_MHz(freq_hi_MHz_),
	dt_sample(dt_sample_),
	sample_rms(sample_rms_),
	randomize_weights(randomize_weights_),
	rng(std::random_device{}()),
	gdist(0, sample_rms_)
    {
	this->nfreq = nfreq_;
	this->nt_chunk = nt_chunk_;

	// Sanity-checking.
	rf_assert(nfreq > 0);
	rf_assert(nt_tot > 0);
	rf_assert(freq_lo_MHz > 0.0);
	rf_assert(freq_lo_MHz < freq_hi_MHz);
	rf_assert(dt_sample > 0.0);
	rf_assert(sample_rms >= 0.0);
	rf_assert(nt_chunk >= 0);

	// Default nt_chunk
	if (nt_chunk == 0)
	    nt_chunk = 1024;

	// Improved RNG seeding.
	rng.discard(700000);
    }

    virtual ~gaussian_noise_stream() { }

    virtual void _bind_stream(Json::Value &json_attrs) override
    {
	json_attrs["freq_lo_MHz"] = this->freq_lo_MHz;
	json_attrs["freq_hi_MHz"] = this->freq_hi_MHz;
	json_attrs["dt_sample"] = this->dt_sample;
    }

    virtual bool _fill_chunk(float *intensity, ssize_t istride, float *weights, ssize_t wstride, ssize_t pos) override
    {
	// Number of Gaussian random samples to be written.
	ssize_t nt = nt_chunk;
	nt = min(nt, nt_tot - pos);
	nt = max(nt, ssize_t(0));

	// Fill the intensity array with Gaussian random numbers, and initialize the weights to 1.
	for (ssize_t ifreq = 0; ifreq < nfreq; ifreq++) {
	    for (ssize_t it = 0; it < nt; it++) {
		intensity[ifreq*istride + it] = gdist(rng);
		weights[ifreq*wstride + it] = randomize_weights ? std::uniform_real_distribution<>()(rng) : 1.0;
	    }
	}

	// The return value from _fill_stream() should be 'true' normally, or 'false' if end-of-stream has been reached.
	if (nt == nt_chunk)
	    return true;

	for (ssize_t ifreq = 0; ifreq < nfreq; ifreq++) {
	    memset(intensity + ifreq*istride + nt, 0, (nt_chunk-nt) * sizeof(float));
	    memset(weights + ifreq*wstride + nt, 0, (nt_chunk-nt) * sizeof(float));
	}

	return false;
    }

    virtual Json::Value jsonize() const override
    {
	Json::Value ret;

	ret["class_name"] = "gaussian_noise_stream";
	ret["nfreq"] = Json::Int64(nfreq);
	ret["nt_tot"] = Json::Int64(nt_tot);
	ret["freq_lo_MHz"] = freq_lo_MHz;
	ret["freq_hi_MHz"] = freq_hi_MHz;
	ret["dt_sample"] = dt_sample;
	ret["sample_rms"] = sample_rms;
	ret["nt_chunk"] = Json::Int64(this->get_prebind_nt_chunk());
	ret["randomize_weights"] = randomize_weights;

	return ret;
    }

    static shared_ptr<gaussian_noise_stream> from_json(const Json::Value &j)
    {
	ssize_t nfreq = ssize_t_from_json(j, "nfreq");
	ssize_t nt_tot = ssize_t_from_json(j, "nt_tot");
	double freq_lo_MHz = double_from_json(j, "freq_lo_MHz");
	double freq_hi_MHz = double_from_json(j, "freq_hi_MHz");
	double dt_sample = double_from_json(j, "dt_sample");
	double sample_rms = double_from_json(j, "sample_rms");
	ssize_t nt_chunk = ssize_t_from_json(j, "nt_chunk");
	bool randomize_weights = bool_from_json(j, "randomize_weights");

	return make_shared<gaussian_noise_stream> (nfreq, nt_tot, freq_lo_MHz, freq_hi_MHz, dt_sample, sample_rms, nt_chunk, randomize_weights);
    }
};


namespace {
    struct _init {
	_init() {
	    pipeline_object::register_json_deserializer("gaussian_noise_stream", gaussian_noise_stream::from_json);
	}
    } init;
}


// Factory function which returns a shared pointer to the wi_stream.  (Using a factory function
// means that details of the base class 'gaussian_noise_stream' and others can be hidden in this 
// source file to avoid overpopulating rf_pipelines.hpp with definitions of many classes, but
// this is just a preference!)

shared_ptr<wi_stream> make_gaussian_noise_stream(ssize_t nfreq, ssize_t nt_tot, double freq_lo_MHz, double freq_hi_MHz, double dt_sample, double sample_rms, ssize_t nt_chunk, bool randomize_weights)
{
    return make_shared<gaussian_noise_stream> (nfreq, nt_tot, freq_lo_MHz, freq_hi_MHz, dt_sample, sample_rms, nt_chunk, randomize_weights);
}


}   // namespace rf_pipelines
