// FIXME: currently, there are two versions of the bonsai_dedisperser, written in python and C++.
// From python, they are constructed as 'bonsai_dedisperser' and 'bonsai_dedisperser_cpp' respectively.
// In the pipeline json output, they are represented as 'bonsai_dedisperser_python' and 'bonsai_dedisperser_cpp'.
// The two versions of the bonsai_dedisperser will be combined eventually!

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
    virtual void _start_pipeline(Json::Value &json_attrs) override;
    virtual void _end_pipeline(Json::Value &json_output) override;
    virtual void _deallocate() override;
    
    virtual Json::Value jsonize() const override;
    static shared_ptr<bonsai_dedisperser> from_json(const Json::Value &j);
};


bonsai_dedisperser::bonsai_dedisperser(const shared_ptr<bonsai::dedisperser> &dp, const shared_ptr<bonsai::global_max_tracker> &tp) :
    wi_transform("bonsai_dedisperser_cpp"),
    dedisperser(dp), 
    max_tracker(tp)
{ 
    // initialize members of wi_transform base class
    this->name = "bonsai_dedisperser_cpp(" + dp->config.name + ")";
    this->nt_chunk = dedisperser->nt_chunk;
}


void bonsai_dedisperser::_allocate()
{
    // If dedisperser is already allocated, this no-ops.
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


void bonsai_dedisperser::_start_pipeline(Json::Value &json_attrs)
{
    // Pass json_attrs to all bonsai::trigger_processors, by serializing to a string
    // and using bonsai's "opaque_context" string.  If a trigger_processor wants to
    // use the attributes, it will need to deserialize the string to json.

    string s = json_stringify(json_attrs);
    dedisperser->set_opaque_context(s);
}


void bonsai_dedisperser::_end_pipeline(Json::Value &json_output)
{
    dedisperser->end_dedispersion();

    for (int itree = 0; itree < dedisperser->ntrees; itree++)
	json_output["max_dm"].append(dedisperser->max_dm[itree]);

    if (max_tracker) {
	json_output["frb_global_max_trigger"] = max_tracker->global_max_trigger;
	json_output["frb_global_max_trigger_dm"] = max_tracker->global_max_trigger_dm;
	json_output["frb_global_max_trigger_tfinal"] = max_tracker->global_max_trigger_arrival_time;
    }
}


void bonsai_dedisperser::_deallocate()
{
    // If dedisperser is already deallocated, this no-ops.
    this->dedisperser->deallocate();
}


Json::Value bonsai_dedisperser::jsonize() const
{
    Json::Value ret;

    // FIXME currently can't jsonize a bonsai_dedisperser which has not been
    // initialized directly from a config file.  We also implicitly assume
    // that the config_params hasn't been modified after it was constructed.
    // Both of these can be improved!

    if ((dedisperser->config.name == "") || (dedisperser->config.name == "config_params"))
	_throw("currently, a bonsai_transform_cpp cannot be jsonized unless it was initialized from a config file");
    if (dedisperser->get_nprocessors() > 0)
	_throw("currently, a bonsai_transform_cpp cannot be jsonized unless it has no trigger_processors (in particular, 'track_global_max' and 'hdf5_output_file' must be disabled)");

    ret["class_name"] = "bonsai_dedisperser_cpp";
    ret["config_filename"] = dedisperser->config.name;

    // FIXME the next code block will become bonsai::dedisperser::initializer::jsonize(),
    // as soon as bonsai gets jsoncpp as an (optional?) dependency.

    Json::Value jini;
    jini["fill_rfi_mask"] = dedisperser->ini_params.fill_rfi_mask;
    jini["allocate"] = dedisperser->ini_params.allocate;
    jini["verbosity"] = dedisperser->ini_params.verbosity;
    jini["file_type"] = dedisperser->ini_params.file_type;
    jini["use_analytic_normalization"] = dedisperser->ini_params.use_analytic_normalization;
    jini["estimate_cumulative_variance"] = dedisperser->ini_params.estimate_cumulative_variance;
    jini["use_unnormalized_triggers"] = dedisperser->ini_params.use_unnormalized_triggers;
    jini["force_unnormalized_triggers"] = dedisperser->ini_params.force_unnormalized_triggers;

    ret["initializer"] = jini;
    return ret;
}


shared_ptr<bonsai_dedisperser> bonsai_dedisperser::from_json(const Json::Value &j)
{
    bonsai::dedisperser::initializer ini_params;

    const Json::Value &jini = j["initializer"];
    ini_params.fill_rfi_mask = bool_from_json(jini, "fill_rfi_mask");
    ini_params.allocate = bool_from_json(jini, "allocate");
    ini_params.verbosity = int_from_json(jini, "verbosity");
    ini_params.file_type = string_from_json(jini, "file_type");
    ini_params.use_analytic_normalization = bool_from_json(jini, "use_analytic_normalization");
    ini_params.estimate_cumulative_variance = bool_from_json(jini, "estimate_cumulative_variance");
    ini_params.force_unnormalized_triggers = bool_from_json(jini, "force_unnormalized_triggers");

    string config_filename = string_from_json(j, "config_filename");
    auto dp = make_shared<bonsai::dedisperser> (config_filename, ini_params);
    
    shared_ptr<bonsai::global_max_tracker> tp;  // empty
    return make_shared<bonsai_dedisperser> (dp, tp);
}


namespace {
    struct _init {
	_init() {
	    pipeline_object::register_json_deserializer("bonsai_dedisperser_cpp", bonsai_dedisperser::from_json);
	}
    } init;
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
	auto t = bonsai::make_trigger_hdf5_writer(ini_params.hdf5_output_filename, ini_params.nt_per_hdf5_file);
	dedisperser->add_processor(t);
    }

    return make_shared<bonsai_dedisperser> (dedisperser, max_tracker);
}


shared_ptr<wi_transform> make_bonsai_dedisperser(const shared_ptr<bonsai::dedisperser> &dp)
{
    if (!dp)
	throw runtime_error("rf_pipelines: empty shared_ptr<bonsai::dedisperser> was passed to make_bonsai_dedisperser()");
    
    shared_ptr<bonsai::global_max_tracker> tp;  // empty
    return make_shared<bonsai_dedisperser> (dp, tp);
}


#endif  // HAVE_BONSAI

}  // namespace rf_pipelines
