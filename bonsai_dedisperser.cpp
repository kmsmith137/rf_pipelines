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

shared_ptr<wi_transform> make_bonsai_dedisperser(const std::string &config_hdf5_filename, const std::string &output_hdf5_filename, int ibeam)
{
    throw runtime_error("make_bonsai_dedisperser() was called, but this rf_pipelines instance was compiled without bonsai");
}

#else  // HAVE_BONSAI


struct bonsai_dedisperser : public wi_transform {
    shared_ptr<bonsai::dedisperser> base;
    string trigger_filename;

    bonsai_dedisperser(const string &config_hdf5_filename, const string &output_hdf5_filename, int ibeam);

    virtual void set_stream(const wi_stream &stream);
    virtual void start_substream(int isubstream, double t0);
    virtual void process_chunk(double t0, double t1, float *intensity, float *weights, int stride, float *pp_intensity, float *pp_weights, int pp_stride);
    virtual void end_substream();
};


bonsai_dedisperser::bonsai_dedisperser(const string &config_hdf5_filename, const string &output_hdf5_filename, int ibeam)
{
    if (!endswith(config_hdf5_filename,".hdf5") && !endswith(config_hdf5_filename,".h5"))
	cerr << "rf_pipelines: warning: bonsai config filename doesn't end with .h5 or .hdf5, note that the bonsai_dedisperser requires an hdf5 file created with bonsai-mkweight\n";
    if (output_hdf5_filename.size() &&  !endswith(output_hdf5_filename,".hdf5") && !endswith(output_hdf5_filename,".h5"))
	cerr << "rf_pipelines: warning: bonsai output filename doesn't end with .h5 or .hdf5\n";

    bool init_weights = true;
    bonsai::config_params cp(config_hdf5_filename, init_weights);

    this->base = make_shared<bonsai::dedisperser> (cp, ibeam, init_weights);
    this->trigger_filename = output_hdf5_filename;

    // initialize members of wi_transform base class
    this->nfreq = base->nchan;
    this->nt_chunk = base->nt_data;
    this->nt_postpad = 0;
    this->nt_prepad = 0;
}


void bonsai_dedisperser::set_stream(const wi_stream &stream)
{
    // Check that stream params match bonsai config

    if (stream.nfreq != this->nfreq)
	throw runtime_error("rf_transforms: value of 'nfreq' in stream doesn't match value in bonsai config");
    if (reldist(stream.freq_lo_MHz, base->freq_lo_MHz) > 1.0e-4)
	throw runtime_error("rf_transforms: value of 'freq_lo_MHz' in stream doesn't match value in bonsai config");
    if (reldist(stream.freq_hi_MHz, base->freq_hi_MHz) > 1.0e-4)
	throw runtime_error("rf_transforms: value of 'freq_hi_MHz' in stream doesn't match value in bonsai config");
    if (reldist(stream.dt_sample, base->dt_sample) > 1.0e-3)
	throw runtime_error("rf_transforms: value of 'dt_sample' in stream doesn't match value in bonsai config");
}


void bonsai_dedisperser::start_substream(int isubstream, double t0)
{
    if (isubstream > 0)
	throw runtime_error("bonsai_dedisperser: currently can't process a stream which defines multiple substreams");

    bool clobber = true;
    if (trigger_filename.size())
	base->start_trigger_file(trigger_filename, clobber);

    base->spawn_slave_threads();
}


void bonsai_dedisperser::process_chunk(double t0, double t1, float *intensity, float *weights, int stride, float *pp_intensity, float *pp_weights, int pp_stride)
{
    // Note: rf_pipelines and bonsai use the same frequency channel ordering (highest-to-lowest), so we can pass the arrays and stride "as is"
    base->run(intensity, weights, stride);
}


void bonsai_dedisperser::end_substream()
{
    if (trigger_filename.size())
	base->end_trigger_file();

    base->terminate();
}


shared_ptr<wi_transform> make_bonsai_dedisperser(const std::string &config_hdf5_filename, const std::string &output_hdf5_filename, int ibeam)
{
    return make_shared<bonsai_dedisperser> (config_hdf5_filename, output_hdf5_filename, ibeam);
}

#endif  // HAVE_BONSAI

}  // namespace rf_pipelines
