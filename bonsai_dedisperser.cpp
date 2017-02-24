#include "rf_pipelines_internals.hpp"

#ifdef HAVE_BONSAI
#include <bonsai.hpp>
#endif

using namespace std;

namespace rf_pipelines {
#if 0
};  // pacify emacs c-mode
#endif


#ifndef HAVE_BONSAI

shared_ptr<wi_transform> make_bonsai_dedisperser(const std::string &config_filename, bool track_global_max)
{
    throw runtime_error("make_bonsai_dedisperser() was called, but this rf_pipelines instance was compiled without bonsai");
}

#else  // HAVE_BONSAI


struct bonsai_dedisperser : public wi_transform {
    unique_ptr<bonsai::dedisperser> dedisperser;
    shared_ptr<bonsai::global_max_tracker> max_tracker;

    bonsai::config_params config;

    bonsai_dedisperser(const string &config_filename, bool track_global_max);

    virtual void set_stream(const wi_stream &stream) override;
    virtual void start_substream(int isubstream, double t0) override;
    virtual void process_chunk(double t0, double t1, float *intensity, float *weights, ssize_t stride, float *pp_intensity, float *pp_weights, ssize_t pp_stride) override;
    virtual void end_substream() override;
};



bonsai_dedisperser::bonsai_dedisperser(const string &config_filename, bool track_global_max)
{
    auto initializer = bonsai::config_initializer::make(config_filename);

    bonsai::config_params config(*initializer);
    this->dedisperser = make_unique<bonsai::dedisperser> (config);
    
    if (initializer->analytic_variance_readable)
	this->dedisperser->read_analytic_variance(*initializer);
    else
	this->dedisperser->compute_analytic_variance();

    if (track_global_max) {
	this->max_tracker = make_shared<bonsai::global_max_tracker> ();
	this->dedisperser->add_processor(max_tracker);
    }

    // initialize members of wi_transform base class
    this->name = "bonsai_dedisperser(" + config_filename + ")";
    this->nfreq = dedisperser->nfreq;
    this->nt_chunk = dedisperser->nt_chunk;
    this->nt_postpad = 0;
    this->nt_prepad = 0;

    // FIXME: write more config info?
    for (int itree = 0; itree < dedisperser->ntrees; itree++)
	this->json_persistent["max_dm"].append(dedisperser->max_dm[itree]);
}


void bonsai_dedisperser::set_stream(const wi_stream &stream)
{
    // Check that stream params match bonsai config

    if (stream.nfreq != dedisperser->nfreq)
	throw runtime_error("rf_transforms: value of 'nfreq' in stream doesn't match value in bonsai config");
    if (reldist(stream.freq_lo_MHz, dedisperser->freq_lo_MHz) > 1.0e-4)
	throw runtime_error("rf_transforms: value of 'freq_lo_MHz' in stream doesn't match value in bonsai config");
    if (reldist(stream.freq_hi_MHz, dedisperser->freq_hi_MHz) > 1.0e-4)
	throw runtime_error("rf_transforms: value of 'freq_hi_MHz' in stream doesn't match value in bonsai config");
    if (reldist(stream.dt_sample, dedisperser->dt_sample) > 1.0e-3)
	throw runtime_error("rf_transforms: value of 'dt_sample' in stream doesn't match value in bonsai config");
}


void bonsai_dedisperser::start_substream(int isubstream, double t0)
{
    // Note: it's OK to reuse the same bonsai_dedisperser object between multiple pipeline runs,
    // but we currently can't handle the case of a run which defines multiple substreams.
    if (isubstream > 0)
	throw runtime_error("bonsai_dedisperser: currently can't process a stream which defines multiple substreams");
}


void bonsai_dedisperser::process_chunk(double t0, double t1, float *intensity, float *weights, ssize_t stride, float *pp_intensity, float *pp_weights, ssize_t pp_stride)
{
    // Note: rf_pipelines and bonsai use the same frequency channel ordering (highest-to-lowest), so we can pass the arrays and stride "as is"
    dedisperser->run(intensity, weights, stride);
}


void bonsai_dedisperser::end_substream()
{
    dedisperser->end_dedispersion();

    if (max_tracker) {
	this->json_per_substream["frb_global_max_trigger"] = max_tracker->global_max_trigger;
	this->json_per_substream["frb_global_max_trigger_dm"] = max_tracker->global_max_trigger_dm;
	this->json_per_substream["frb_global_max_trigger_tfinal"] = max_tracker->global_max_trigger_arrival_time;
    }
}


shared_ptr<wi_transform> make_bonsai_dedisperser(const std::string &config_filename, bool track_global_max)
{
    return make_shared<bonsai_dedisperser> (config_filename, track_global_max);
}

#endif  // HAVE_BONSAI

}  // namespace rf_pipelines
