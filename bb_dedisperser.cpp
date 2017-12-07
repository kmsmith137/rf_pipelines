// FIXME need to assert somehow that total stream length is a multiple of nt_chunk!
// OK for now since we only use the bb_dedisperser through the frb_olympics, which
// has its own mechanism for placing this assert.  However it exposes a general
// limitation of rf_pipelines!
//
// FIXME artificial limitation: hardcoded quantization parameters, data must have (mean,rms) ~ (0,1).

#include "rf_pipelines_internals.hpp"
#include <dedisp.h>

using namespace std;

namespace rf_pipelines {
#if 0
};  // pacify emacs c-mode
#endif

// The 'bb_dedisperser' wrapper class converts Ben Barsdell's 'dedisp' GPU code
// to the rf_pipelines API.  It does this by subclassing rf_pipelines::wi_transform,
// and defining the approporiate virtual functions.


// -------------------------------------------------------------------------------------------------
//
// Helpers


// Just a wrapper around dedisp_get_error_string) which returns std::string instead of (const char *).
static string dedisp_errmsg(dedisp_error error)
{
    return dedisp_get_error_string(error);
}


inline uint8_t quantize(float x)
{
    // FIXME hardcoded quantization scales.
    // Assumes (mean, rms) of the data is close to (0, 1).

    x = 10.0f * x + 127.0f;
    x = max(x, 0.5f);
    x = min(x, 255.5f);

    return uint8_t(x);
}


// -------------------------------------------------------------------------------------------------


struct bb_dedisperser_initializer {
    // FIXME currently, nt_in must be known in advance!
    ssize_t nt_in = 0;
    int verbosity = 1;

    // Used to determine list of trial DM's.
    // Negative initializers are to catch unintialized values.
    // FIXME is it safe to set dm_start==0 or pulse_width_ms==0?  (Just need to check dedisp source code.)
    double dm_start = -1.0;
    double dm_end = -1.0;
    double dm_tol = -1.0;
    double pulse_width_ms = -1.0;
};


struct bb_dedisperser : public wi_transform {
    // Initialized in constructor
    const bb_dedisperser_initializer ini_params;
    const ssize_t nt_in;   // same as ini_params.nt_in

    // Initialized in _bind_transform()
    double freq_lo_MHz = 0.0;
    double freq_hi_MHz = 0.0;
    double dt_sample = 0.0;

    // Initialized in _allocate()
    ssize_t ndm = 0;
    ssize_t nt_out = 0;
    ssize_t max_delay = 0;
    dedisp_plan plan = NULL;
    uptr<uint8_t> in8;    // shape (nt_in, nfreq)
    uptr<float> out32;    // shape (ndm, nt_out)

    // Set in _start_pipeline(), _process_chunk(), _end_pipeline().
    bool runflag = false;
    float max_trigger = 0.0;
    float max_trigger_dm = 0.0;
    float max_trigger_tfinal = 0.0;

    bb_dedisperser(const bb_dedisperser_initializer &ini_params);

    virtual void _bind_transform(Json::Value &json_attrs) override;
    virtual void _start_pipeline(Json::Value &json_attrs) override;
    virtual void _process_chunk(float *intensity, ssize_t istride, float *weights, ssize_t wstride, ssize_t pos) override;
    virtual void _end_pipeline(Json::Value &json_output) override;
    
    virtual void _allocate() override;
    virtual void _deallocate() override;
};


bb_dedisperser::bb_dedisperser(const bb_dedisperser_initializer &ini_params_) :
    wi_transform("bb_dedisperser"),
    ini_params(ini_params_),
    nt_in(ini_params_.nt_in)
{
    if (ini_params.dm_start <= 0.0)
	throw runtime_error("rf_pipelines::bb_dedisperser constructor: expected dm_start > 0.0");
    if (ini_params.dm_end <= ini_params.dm_start)
	throw runtime_error("rf_pipelines::bb_dedisperser constructor: expected dm_end > dm_start");
    if (ini_params.dm_tol <= 0.0)
	throw runtime_error("rf_pipelines::bb_dedisperser constructor: expected dm_tol > 0.0");
    if (ini_params.pulse_width_ms < 0.0)
	throw runtime_error("rf_pipelines::bb_dedisperser constructor: expected pulse_width_ms > 0.0");
    if (ini_params.nt_in <= 0)
	throw runtime_error("rf_pipelines::bb_dedisperser constructor: expected nt_in > 0");
}


void bb_dedisperser::_bind_transform(Json::Value &json_attrs)
{
    if (!json_attrs.isMember("freq_lo_MHz") || !json_attrs.isMember("freq_hi_MHz"))
	throw runtime_error("bonsai_dedisperser: expected json_attrs to contain members 'freq_lo_MHz' and 'freq_hi_MHz'");

    if (!json_attrs.isMember("dt_sample"))
	throw runtime_error("bonsai_dedisperser: expected json_attrs to contain member 'dt_sample'");
    
    this->freq_lo_MHz = json_attrs["freq_lo_MHz"].asDouble();
    this->freq_hi_MHz = json_attrs["freq_hi_MHz"].asDouble();
    this->dt_sample = json_attrs["dt_sample"].asDouble();

    if (freq_lo_MHz <= 0.0)
	throw runtime_error("rf_pipelines::bb_dedisperser::_bind_transform(): expected freq_lo_MHz > 0.0");
    if (freq_hi_MHz <= freq_lo_MHz)
	throw runtime_error("rf_pipelines::bb_dedisperser::_bind_transform(): expected freq_lo_MHz < freq_hi_MHz");
    if (dt_sample <= 0.0)
	throw runtime_error("rf_pipelines::bb_dedisperser::_bind_transform(): expected dt_sample > 0.0");
}


void bb_dedisperser::_allocate()
{
    rf_assert(this->plan == nullptr);

    int device_idx = 0;    // FIXME should determine this somehow
    dedisp_error error;
    
    error = dedisp_set_device(device_idx);
    if (error != DEDISP_NO_ERROR)
	throw runtime_error("ERROR: Could not set GPU device: " + dedisp_errmsg(error));

    // Channel width in MHz.
    // Should be negative (if channels ordered highest-to-lowest).
    double df = -float(this->freq_hi_MHz - this->freq_lo_MHz) / this->nfreq;
    
    error = dedisp_create_plan(&plan, this->nfreq, this->dt_sample, this->freq_hi_MHz, df);
    if (error != DEDISP_NO_ERROR)
	throw runtime_error("ERROR: Could not create dedispersion plan: " + dedisp_errmsg(error));

    error = dedisp_generate_dm_list(plan, ini_params.dm_start, ini_params.dm_end, ini_params.pulse_width_ms, ini_params.dm_tol);
    if (error != DEDISP_NO_ERROR)
	throw runtime_error("ERROR: Failed to generate dm list: " + dedisp_errmsg(error));

    if (ini_params.verbosity >= 1)
	cout << "bb_dedisperser: number of trial DM's = " << dedisp_get_dm_count(plan);

    this->ndm = dedisp_get_dm_count(plan);
    this->max_delay = dedisp_get_max_delay(plan);
    this->nt_out = nt_in - max_delay;
	
    if (nt_out <= 0)
	throw runtime_error("bb_dedisperser: length of acquisition is shorter than max dedispersion delay");

    this->in8 = make_uptr<uint8_t> (nt_in * nfreq);  // ordering (time, freq)
    this->out32 = make_uptr<float> (ndm * nt_out);   // ordering (dm, time)
}


void bb_dedisperser::_start_pipeline(Json::Value &json_attrs)
{
    this->runflag = false;
    this->max_trigger = 0.0;
    this->max_trigger_dm = 0.0;
    this->max_trigger_tfinal = 0.0;
}
	

void bb_dedisperser::_process_chunk(float *intensity, ssize_t istride, float *weights, ssize_t wstride, ssize_t pos)
{
    rf_assert(this->plan != nullptr);
    
    ssize_t nt_valid = nt_in - pos;
    nt_valid = min(nt_valid, nt_chunk);
    nt_valid = max(nt_valid, ssize_t(0));

    // Check weights.

    for (int ifreq = 0; ifreq < nfreq; ifreq++) {
	for (int it = 0; it < nt_valid; it++)
	    rf_assert(weights[ifreq*wstride + it] == 1.0);
	for (int it = nt_valid; it < nt_chunk; it++)
	    rf_assert(weights[ifreq*wstride + it] == 0.0);
    }

    // Quantize and copy.
    
    for (int ifreq = 0; ifreq < nfreq; ifreq++)
	for (int it = 0; it < nt_valid; it++)
	    in8[(pos+it)*nfreq + ifreq] = quantize(intensity[ifreq*istride+it]);

    if (runflag || (pos + nt_chunk < nt_in))
	return;
    
    // The dedisperser will be run in this iteration of _process_chunk().
    
    dedisp_error error = dedisp_execute(plan, nt_in, in8.get(), 8, (dedisp_byte *) (out32.get()), sizeof(float), DEDISP_USE_DEFAULT);

    if (error != DEDISP_NO_ERROR)
	throw runtime_error("rf_pipelines::bb_dedisperser: dedisp_execute() failed: " + dedisp_errmsg(error));
    
    // Determine rms, max, and argmax of 'out32' array
    // Mean zero is assumed!  (Currently includes a debug-print to confirm this by computing the mean and printing it.)

    float maxval = out32[0];
    ssize_t argmax = 0;
    double mean = 0.0;
    double rms = 0.0;

    for (ssize_t i = 0; i < ndm * nt_out; i++) {
	float x = out32[i];
	
	if (maxval > x) {
	    maxval = x;
	    argmax = i;
	}

	mean += x;
	rms += x*x;
    }

    // This can be removed after the mean is confirmed to be zero.
    cout << "XXX bb_dedisperser: mean: " << (mean/rms) << " sigmas" << endl;
    
    // Convert 'argmax' index to pair (idm, iout)
    ssize_t idm = argmax / nt_out;
    ssize_t it = argmax % nt_out;

    // Note: this pointer does not need to be freed (but is only temporarily valid).
    const float *dmlist = dedisp_get_dm_list(plan);

    this->max_trigger = maxval / rms;   // convert to "sigmas"
    this->max_trigger_dm = dmlist[idm];
    this->max_trigger_tfinal = it * dt_sample;    // FIXME tinitial -> tfinal
    this->runflag = true;
}


void bb_dedisperser::_end_pipeline(Json::Value &json_outputs)
{
    if (!runflag)
	throw runtime_error("bonsai::bb_dedisperser: stream ended early (_end_pipeline() called before 'nt_in' samples were processed)");

    json_outputs["frb_global_max_trigger"] = this->max_trigger;
    json_outputs["frb_global_max_trigger_dm"] = this->max_trigger_dm;
    json_outputs["frb_global_max_trigger_tfinal"] = this->max_trigger_tfinal;

    this->max_trigger = 0.0;
    this->max_trigger_dm = 0.0;
    this->max_trigger_tfinal = 0.0;
    this->runflag = false;
}


void bb_dedisperser::_deallocate()
{
    rf_assert(this->plan != nullptr);

    dedisp_destroy_plan(this->plan);
    this->plan = nullptr;    
    this->in8.reset();
    this->out32.reset();
}


// -------------------------------------------------------------------------------------------------
//
// make_bb_dedisperser()
//
// This interface (function which returns a shared pointer to the wi_transform base class) 
// allows the implementation of struct bb_dedisperser to be "hidden" from the rest of rf_pipelines.


std::shared_ptr<wi_transform> make_bb_dedisperser(const bb_dedisperser_initializer &ini_params)
{
    return std::make_shared<bb_dedisperser> (ini_params);
}


}  // namespace rf_pipelines
