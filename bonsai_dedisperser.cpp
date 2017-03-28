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

shared_ptr<wi_transform> make_bonsai_dedisperser(const string &config_filename, const bonsai_initializer &ini_params)
{
    throw runtime_error("make_bonsai_dedisperser() was called, but this rf_pipelines instance was compiled without bonsai");
}

#else  // HAVE_BONSAI


struct bonsai_dedisperser : public wi_transform {
    shared_ptr<bonsai::dedisperser> dedisperser;
    shared_ptr<bonsai::global_max_tracker> max_tracker;
    bool deallocate_between_substreams;

    // Note: if 'tp' is a nonempty pointer, then caller is responsible for calling dp->add_processor(tp).
    // This is a little "fragile", but since this is an internal interface, I didn't bother improving it!
    bonsai_dedisperser(const shared_ptr<bonsai::dedisperser> &dp, const shared_ptr<bonsai::global_max_tracker> &tp, bool deallocate_between_substreams);

    virtual void set_stream(const wi_stream &stream) override;
    virtual void start_substream(int isubstream, double t0) override;
    virtual void process_chunk(double t0, double t1, float *intensity, float *weights, ssize_t stride, float *pp_intensity, float *pp_weights, ssize_t pp_stride) override;
    virtual void end_substream() override;
};


bonsai_dedisperser::bonsai_dedisperser(const shared_ptr<bonsai::dedisperser> &dp, const shared_ptr<bonsai::global_max_tracker> &tp, bool deallocate_between_substreams_) :
    dedisperser(dp), 
    max_tracker(tp),
    deallocate_between_substreams(deallocate_between_substreams_)
{ 
    // initialize members of wi_transform base class
    this->name = "bonsai_dedisperser(" + dp->config.name + ")";
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

    if (stream.nfreq != dedisperser->nfreq) {
	throw runtime_error("rf_pipelines: value of 'nfreq' in stream (=" + to_string(stream.nfreq)
			    + ") doesn't match value in bonsai config (=" + to_string(dedisperser->nfreq) + ")");
    }

    if (reldist(stream.freq_lo_MHz, dedisperser->freq_lo_MHz) > 1.0e-4) {
	throw runtime_error("rf_pipelines: value of 'freq_lo_MHz' in stream (=" + to_string(stream.freq_lo_MHz)
			    + ") doesn't match value in bonsai config (=" + to_string(dedisperser->freq_lo_MHz) + ")");
    }

    if (reldist(stream.freq_hi_MHz, dedisperser->freq_hi_MHz) > 1.0e-4) {
	throw runtime_error("rf_pipelines: value of 'freq_hi_MHz' in stream (=" + to_string(stream.freq_hi_MHz)
			    + ") doesn't match value in bonsai config (=" + to_string(dedisperser->freq_hi_MHz) + ")");
    }

    if (reldist(stream.dt_sample, dedisperser->dt_sample) > 1.0e-3) {
	throw runtime_error("rf_pipelines: value of 'dt_sample' in stream (=" + to_string(stream.dt_sample)
			    + ") doesn't match value in bonsai config (=" + to_string(dedisperser->dt_sample) + ")");
    }
}


void bonsai_dedisperser::start_substream(int isubstream, double t0)
{
    // Note: it's OK to reuse the same bonsai_dedisperser object between multiple pipeline runs,
    // but we currently can't handle the case of a run which defines multiple substreams.
    if (isubstream > 0)
	throw runtime_error("bonsai_dedisperser: currently can't process a stream which defines multiple substreams");

    if (deallocate_between_substreams)
	this->dedisperser->allocate();
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

    if (deallocate_between_substreams)
	this->dedisperser->deallocate();
}


// -------------------------------------------------------------------------------------------------


shared_ptr<wi_transform> make_bonsai_dedisperser(const std::string &config_filename, const bonsai_initializer &ini_params)
{
    bonsai::dedisperser::initializer ini2;

    ini2.file_type = ini_params.file_type;
    ini2.verbosity = ini_params.verbosity;
    ini2.allocate = !ini_params.deallocate_between_substreams;
    ini2.use_analytic_normalization = ini_params.use_analytic_normalization;

    shared_ptr<bonsai::dedisperser> dedisperser = make_shared<bonsai::dedisperser> (config_filename, ini2);
    shared_ptr<bonsai::global_max_tracker> max_tracker;

    if (ini_params.track_global_max) {
	max_tracker = make_shared<bonsai::global_max_tracker> (ini_params.dm_min, ini_params.dm_max);
	dedisperser->add_processor(max_tracker);
    }

    if (ini_params.hdf5_output_filename.size() > 0) {
	auto t = make_shared<bonsai::trigger_hdf5_file_writer> (ini_params.hdf5_output_filename, ini_params.nt_per_hdf5_file);
	dedisperser->add_processor(t);
    }

    return make_shared<bonsai_dedisperser> (dedisperser, max_tracker, ini_params.deallocate_between_substreams);
}


shared_ptr<wi_transform> make_bonsai_dedisperser(const shared_ptr<bonsai::dedisperser> &dp)
{
    if (!dp)
	throw runtime_error("rf_pipelines: empty shared_ptr<bonsai::dedisperser> was passed to make_bonsai_dedisperser()");

    shared_ptr<bonsai::global_max_tracker> tp;
    return make_shared<bonsai_dedisperser> (dp, tp, false);
}


#endif  // HAVE_BONSAI

}  // namespace rf_pipelines
