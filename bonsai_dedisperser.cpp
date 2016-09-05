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

shared_ptr<wi_transform> make_bonsai_dedisperser(const std::string &config_hdf5_filename, const std::string &trigger_hdf5_filename, const std::string &trigger_plot_stem, int nt_per_file, int ibeam)
{
    throw runtime_error("make_bonsai_dedisperser() was called, but this rf_pipelines instance was compiled without bonsai");
}

#else  // HAVE_BONSAI


struct bonsai_dedisperser : public wi_transform {
    string config_filename;
    string trigger_hdf5_filename;
    string trigger_plot_stem;
    int nt_per_file;
    int ibeam;

    shared_ptr<bonsai::config_params> config;
    shared_ptr<bonsai::dedisperser> dedisperser;

    bonsai_dedisperser(const string &config_hdf5_filename, const string &trigger_hdf5_filename, const string &trigger_plot_stem, int nt_per_file, int ibeam);

    virtual void set_stream(const wi_stream &stream) override;
    virtual void start_substream(int isubstream, double t0) override;
    virtual void process_chunk(double t0, double t1, float *intensity, float *weights, ssize_t stride, float *pp_intensity, float *pp_weights, ssize_t pp_stride) override;
    virtual void end_substream() override;
};


// A minimal bonsai::dedisperser subclass which ensures that output files (including plots)
// get the "standard" rf_pipelines processing (e.g. filenames end up in json output).
struct my_dedisperser_subclass : public bonsai::dedisperser {
    // FIXME this bare pointer is kind of ugly, any way to improve?
    bonsai_dedisperser *transform;

    my_dedisperser_subclass(bonsai_dedisperser *transform_)
	: bonsai::dedisperser(*transform_->config, transform_->ibeam, true),  // init_weights=true
	  transform(transform_)
    { }

    virtual void _open_trigger_file(const string &basename, const string &datetime0_str, const string &datetime_str)
    {
	string filename = transform->add_file(basename);
	transform->json_misc["trigger_files"].append(filename);
	bonsai::dedisperser::_open_trigger_file(filename, datetime0_str, datetime_str);
    }

    virtual string _make_trigger_plot_filename(int itree, int ifile)
    {
	auto config = transform->config;

	int ndm = this->trigger_plot_ndm[itree];
	int nt_coarse_curr = this->trigger_plot_nt_curr[itree];
	int nt_coarse_max = this->trigger_plot_nt_max[itree];
	int nt_per_trigger = (config->nds[itree] * config->nt_per_trigger[itree]) / config->nups[itree];

	ssize_t it0 = ssize_t(ifile) * ssize_t(nt_coarse_max) * ssize_t(nt_per_trigger);
	ssize_t nt = ssize_t(nt_coarse_curr) * ssize_t(nt_per_trigger);

	string basename = dedisperser::_make_trigger_plot_filename(itree, ifile);
	string filename = transform->add_plot(basename, it0, nt, nt_coarse_curr, ndm, itree);

	return filename;
    }
};


bonsai_dedisperser::bonsai_dedisperser(const string &config_hdf5_filename, const string &trigger_hdf5_filename_, const string &trigger_plot_stem_, int nt_per_file_, int ibeam_) :
    config_filename(config_hdf5_filename),
    trigger_hdf5_filename(trigger_hdf5_filename_),
    trigger_plot_stem(trigger_plot_stem_),
    nt_per_file(nt_per_file_),
    ibeam(ibeam_)
{
    if (!endswith(config_filename,".hdf5") && !endswith(config_filename,".h5"))
	cerr << "rf_pipelines: warning: bonsai config filename doesn't end with .h5 or .hdf5, note that the bonsai_dedisperser requires an hdf5 file created with bonsai-mkweight\n";
    if (trigger_hdf5_filename.size() &&  !endswith(trigger_hdf5_filename,".hdf5") && !endswith(trigger_hdf5_filename,".h5"))
	cerr << "rf_pipelines: warning: bonsai output filename doesn't end with .h5 or .hdf5\n";

    this->config = make_shared<bonsai::config_params> (config_hdf5_filename, true);   // init_weights=true

    // initialize members of wi_transform base class
    this->name = "bonsai_dedisperser(" + config_filename + ")";
    this->nfreq = config->nchan;
    this->nt_chunk = config->nt_data;
    this->nt_postpad = 0;
    this->nt_prepad = 0;

    if (trigger_plot_stem.size() > 0) {
	if (nt_per_file <= 0)
	    throw runtime_error("rf_pipelines::bonsai_dedisperser: if the 'trigger_plot_stem' arg is specified, then 'nt_per_file' must also be specified and > 0");

	for (int itree = 0; itree < config->ntrees; itree++) {
	    int ndm_coarse = config->tree_size[itree] / config->ndm_per_trigger[itree];
	    int nt_per_trigger = (config->nds[itree] * config->nt_per_trigger[itree]) / config->nups[itree];
	    this->add_plot_group("bonsai_tree" + to_string(itree), nt_per_trigger, ndm_coarse);
	}
    }
}


void bonsai_dedisperser::set_stream(const wi_stream &stream)
{
    // Check that stream params match bonsai config

    if (stream.nfreq != config->nchan)
	throw runtime_error("rf_transforms: value of 'nfreq' in stream doesn't match value in bonsai config");
    if (reldist(stream.freq_lo_MHz, config->freq_lo_MHz) > 1.0e-4)
	throw runtime_error("rf_transforms: value of 'freq_lo_MHz' in stream doesn't match value in bonsai config");
    if (reldist(stream.freq_hi_MHz, config->freq_hi_MHz) > 1.0e-4)
	throw runtime_error("rf_transforms: value of 'freq_hi_MHz' in stream doesn't match value in bonsai config");
    if (reldist(stream.dt_sample, config->dt_sample) > 1.0e-3)
	throw runtime_error("rf_transforms: value of 'dt_sample' in stream doesn't match value in bonsai config");
}


void bonsai_dedisperser::start_substream(int isubstream, double t0)
{
    // Note: it's OK to reuse the same bonsai_dedisperser object between multiple pipeline runs,
    // but we currently can't handle the case of a run which defines multiple substreams.
    if (isubstream > 0)
	throw runtime_error("bonsai_dedisperser: currently can't process a stream which defines multiple substreams");

    // FIXME: write more config info?
    // Note that 'json_misc' is per-substream, so we can't put this initialization in set_stream() or the constructor
    for (int itree = 0; itree < config->ntrees; itree++)
	this->json_misc["max_dm"].append(config->max_dm[itree]);

    this->dedisperser = make_shared<my_dedisperser_subclass> (this);
    
    if (trigger_hdf5_filename.size())
	dedisperser->start_trigger_file(this->trigger_hdf5_filename, this->nt_per_file);
    if (trigger_plot_stem.size())
	dedisperser->start_trigger_plots(this->trigger_plot_stem, this->nt_per_file);
    
    dedisperser->global_max_trigger_active = true;
    dedisperser->global_max_trigger = 0.0;
    dedisperser->global_max_trigger_dm = 0.0;
    dedisperser->global_max_trigger_arrival_time = 0.0;
    
    dedisperser->spawn_slave_threads();
}


void bonsai_dedisperser::process_chunk(double t0, double t1, float *intensity, float *weights, ssize_t stride, float *pp_intensity, float *pp_weights, ssize_t pp_stride)
{
    // Note: rf_pipelines and bonsai use the same frequency channel ordering (highest-to-lowest), so we can pass the arrays and stride "as is"
    dedisperser->run(intensity, weights, stride);
}


void bonsai_dedisperser::end_substream()
{
    if (trigger_hdf5_filename.size())
	dedisperser->end_trigger_file();
    if (trigger_plot_stem.size())
	dedisperser->end_trigger_plots();

    this->json_misc["frb_global_max_trigger"] = dedisperser->global_max_trigger;
    this->json_misc["frb_global_max_trigger_dm"] = dedisperser->global_max_trigger_dm;
    this->json_misc["frb_global_max_trigger_tfinal"] = dedisperser->global_max_trigger_arrival_time;

    dedisperser->terminate();
}


shared_ptr<wi_transform> make_bonsai_dedisperser(const std::string &config_hdf5_filename, const std::string &trigger_hdf5_filename, const std::string &trigger_plot_stem, int nt_per_file, int ibeam)
{
    return make_shared<bonsai_dedisperser> (config_hdf5_filename, trigger_hdf5_filename, trigger_plot_stem, nt_per_file, ibeam);
}

#endif  // HAVE_BONSAI

}  // namespace rf_pipelines
