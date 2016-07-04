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
    ssize_t nt_tot;
    double sample_rms;

public:
    gaussian_noise_stream(ssize_t nfreq_, ssize_t nt_tot_, double freq_lo_MHz_, double freq_hi_MHz_, double dt_sample_, double sample_rms_, ssize_t nt_chunk)
    {
	if (nt_chunk == 0)
	    nt_chunk = 1024;   // default

	// these arguments will be sanity-checked elsewhere (in wi_stream::run())
	this->nfreq = nfreq_;
	this->nt_maxwrite = nt_chunk;
	this->freq_lo_MHz = freq_lo_MHz_;
	this->freq_hi_MHz = freq_hi_MHz_;
	this->dt_sample = dt_sample_;

	// these arguments must be sanity-checked here
	this->sample_rms = sample_rms_;
	this->nt_tot = nt_tot_;

	rf_assert(sample_rms >= 0.0);
	rf_assert(nt_chunk > 0);
	rf_assert(nt_tot > 0);
    }

    virtual ~gaussian_noise_stream() { }

    virtual void stream_body(wi_run_state &run_state)
    {
	std::random_device rd;
	std::mt19937 rng(rd());
	std::normal_distribution<float> dist(0, sample_rms);

	run_state.start_substream(0.0);

	ssize_t it0 = 0;
	while (it0 < nt_tot) {
	    ssize_t nt = min(nt_maxwrite, nt_tot-it0);
	    float *intensity;
	    float *weights;
	    ssize_t stride;

	    bool zero_flag = false;
	    run_state.setup_write(nt, intensity, weights, stride, zero_flag);

	    for (ssize_t ifreq = 0; ifreq < nfreq; ifreq++) {
		for (ssize_t it = 0; it < nt; it++) {
		    intensity[ifreq*stride + it] = dist(rng);
		    weights[ifreq*stride + it] = 1.0;
		}
	    }

	    run_state.finalize_write(nt);
	    it0 += nt;
	}

	run_state.end_substream();
    }
};


shared_ptr<wi_stream> make_gaussian_noise_stream(ssize_t nfreq, ssize_t nt_tot, double freq_lo_MHz, double freq_hi_MHz, double dt_sample, double sample_rms, ssize_t nt_chunk)
{
    return make_shared<gaussian_noise_stream> (nfreq, nt_tot, freq_lo_MHz, freq_hi_MHz, dt_sample, sample_rms, nt_chunk);
}


}   // namespace rf_pipelines
