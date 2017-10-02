// FIXME will implement jsonize() after json (de)serialization is added to bonsai.

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

    // Note: if 'tp' is a nonempty pointer, then caller is responsible for calling dp->add_processor(tp).
    // This is a little "fragile", but since this is an internal interface, I didn't bother improving it!
    bonsai_dedisperser(const shared_ptr<bonsai::dedisperser> &dp, const shared_ptr<bonsai::global_max_tracker> &tp);

    virtual void _bind_transform(Json::Value &json_attrs) override;
    virtual void _allocate() override;
    virtual void _process_chunk(float *intensity, ssize_t istride, float *weights, ssize_t wstride, ssize_t pos) override;
    virtual void _end_pipeline(Json::Value &json_output) override;
    virtual void _deallocate() override;
};


bonsai_dedisperser::bonsai_dedisperser(const shared_ptr<bonsai::dedisperser> &dp, const shared_ptr<bonsai::global_max_tracker> &tp) :
    wi_transform("bonsai_dedisperser"),
    dedisperser(dp), 
    max_tracker(tp)
{ 
    // initialize members of wi_transform base class
    this->name = "bonsai_dedisperser(" + dp->config.name + ")";
    this->nfreq = dedisperser->nfreq;
    this->nt_chunk = dedisperser->nt_chunk;

    // FIXME: write more config info?
    //for (int itree = 0; itree < dedisperser->ntrees; itree++)
    //this->json_persistent["max_dm"].append(dedisperser->max_dm[itree]);
}


void bonsai_dedisperser::_allocate()
{
    if (this->dedisperser->state == bonsai::dedisperser::DEALLOCATED)
	this->dedisperser->allocate();
}


void bonsai_dedisperser::_bind_transform(Json::Value &json_attrs)
{
    if (!json_attrs.isMember("freq_lo_MHz") || !json_attrs.isMember("freq_hi_MHz"))
	throw runtime_error("bonsai_dedisperser: expected json_attrs to contain members 'freq_lo_MHz' and 'freq_hi_MHz'");

    if (!json_attrs.isMember("dt_sample"))
	throw runtime_error("bonsai_dedisperser: expected json_attrs to contain member 'dt_sample'");
    
    double freq_lo_MHz = json_attrs["freq_lo_MHz"].asDouble();
    double freq_hi_MHz = json_attrs["freq_hi_MHz"].asDouble();
    double dt_sample = json_attrs["dt_sample"].asDouble();

    // Check that pipeline params match bonsai config

    if (nfreq != dedisperser->nfreq) {
	throw runtime_error("rf_pipelines: value of 'nfreq' in stream (=" + to_string(nfreq)
			    + ") doesn't match value in bonsai config (=" + to_string(dedisperser->nfreq) + ")");
    }

    if (reldist(freq_lo_MHz, dedisperser->freq_lo_MHz) > 1.0e-4) {
	throw runtime_error("rf_pipelines: value of 'freq_lo_MHz' in stream (=" + to_string(freq_lo_MHz)
			    + ") doesn't match value in bonsai config (=" + to_string(dedisperser->freq_lo_MHz) + ")");
    }

    if (reldist(freq_hi_MHz, dedisperser->freq_hi_MHz) > 1.0e-4) {
	throw runtime_error("rf_pipelines: value of 'freq_hi_MHz' in stream (=" + to_string(freq_hi_MHz)
			    + ") doesn't match value in bonsai config (=" + to_string(dedisperser->freq_hi_MHz) + ")");
    }

    if (reldist(dt_sample, dedisperser->dt_sample) > 1.0e-3) {
	throw runtime_error("rf_pipelines: value of 'dt_sample' in stream (=" + to_string(dt_sample)
			    + ") doesn't match value in bonsai config (=" + to_string(dedisperser->dt_sample) + ")");
    }
}

void bonsai_dedisperser::_process_chunk(float *intensity, ssize_t istride, float *weights, ssize_t wstride, ssize_t pos)
{
    // Note: rf_pipelines and bonsai use the same frequency channel ordering (highest-to-lowest), so we can pass the arrays and stride "as is"
    dedisperser->run(intensity, istride, weights, wstride);
}


void bonsai_dedisperser::_end_pipeline(Json::Value &json_output)
{
    dedisperser->end_dedispersion();

    if (max_tracker) {
	json_output["frb_global_max_trigger"] = max_tracker->global_max_trigger;
	json_output["frb_global_max_trigger_dm"] = max_tracker->global_max_trigger_dm;
	json_output["frb_global_max_trigger_tfinal"] = max_tracker->global_max_trigger_arrival_time;
    }
}


void bonsai_dedisperser::_deallocate()
{
    this->dedisperser->deallocate();
}


// -------------------------------------------------------------------------------------------------


shared_ptr<wi_transform> make_bonsai_dedisperser(const std::string &config_filename, const bonsai_initializer &ini_params)
{
    bonsai::dedisperser::initializer ini2;

    ini2.fill_rfi_mask = ini_params.fill_rfi_mask;
    ini2.file_type = ini_params.file_type;
    ini2.verbosity = ini_params.verbosity;
    ini2.use_analytic_normalization = ini_params.use_analytic_normalization;
    ini2.allocate = false;  // postpone allocation until bonsai_dedisperser::_allocate()

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

    return make_shared<bonsai_dedisperser> (dedisperser, max_tracker);
}


shared_ptr<wi_transform> make_bonsai_dedisperser(const shared_ptr<bonsai::dedisperser> &dp)
{
    if (!dp)
	throw runtime_error("rf_pipelines: empty shared_ptr<bonsai::dedisperser> was passed to make_bonsai_dedisperser()");

    shared_ptr<bonsai::global_max_tracker> tp;
    return make_shared<bonsai_dedisperser> (dp, tp);
}


#endif  // HAVE_BONSAI

}  // namespace rf_pipelines

