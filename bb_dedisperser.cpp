// FIXME need to assert somehow that total stream length is a multiple of nt_chunk!
// OK for now since we only use the bb_dedisperser through the frb_olympics, which
// has its own mechanism for placing this assert.  However it exposes a general
// limitation of rf_pipelines!
//
// Artificial limitation: hardcoded quantization parameters, data must have (mean,rms) ~ (0,1).

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


static string dedisp_errmsg(dedisp_error error)
{
    return dedisp_get_error_string(error);    // (const char *) -> std::string
}


/* byte2quant is a utility function which coverts input to unsigned range */
inline dedisp_float tounsigned(dedisp_float f)
{
    dedisp_float v = f + 127.5f;
    dedisp_float r = v;
    if (v>255.0) {
	r = (dedisp_byte)255;
    }
    else if (v< 0.0f) {
	r = (dedisp_byte)0;
    }
    else {
	//r = (dedisp_byte)roundf(f);
    } 
    return r;
}

inline dedisp_byte bytetoquant(dedisp_float f)
{
	dedisp_float v = f + 127.5f;
	dedisp_byte r;
	if (v>255.0) 
	{
		r = (dedisp_byte)255;
	}
	else if (v<0.0f)
	{
		r = (dedisp_byte)0;
	}
	else
	{
		r = (dedisp_byte)roundf(v);
	}
	return r;
}


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

    // Initialized in _bind_transform()
    double freq_lo_MHz = 0.0;
    double freq_hi_MHz = 0.0;
    double dt_sample = 0.0;

    // Initialized in _allocate()
    ssize_t ndm = 0;
    ssize_t nt_out = 0;
    ssize_t max_delay = 0;
    dedisp_plan plan = NULL;
    uptr_t<uint8_t> in8;    // shape (nt_in, nfreq)
    uptr_t<float> out32;    // shape (ndm, nt_out)

    float max_trigger = 0.0;
    float max_trigger_dm = 0.0;
    float max_trigger_tfinal = 0.0;

    bb_dedisperser(const bb_dedisperser_initializer &ini_params);

    virtual void _bind_transform(Json::Value &json_attrs) override;
    virtual void _process_chunk(float *intensity, ssize_t istride, float *weights, ssize_t wstride, ssize_t pos) override;
    
    virtual void _allocate() override;
    virtual void _deallocate() override;
};


bb_dedisperser::bb_dedisperser(const bb_dedisperser_initializer &ini_params_) :
    wi_transform("bb_dedisperser"),
    ini_params(ini_params_)
{
    if (ini_params.dm_start <= 0.0)
	throw runtime_error("rf_pipelines::bb_dedisperser constructor: expected dm_start > 0.0");
    if (ini_params.dm_end <= ini_params.dm_start)
	throw runtime_error("rf_pipelines::bb_dedisperser constructor: must have dm_end > dm_start");
    if (ini_params.dm_tol <= 0.0)
	throw runtime_error("rf_pipelines::bb_dedisperser constructor: must have dm_tol > 0.0");
    if (ini_params.pulse_width_ms < 0.0)
	throw runtime_error("rf_pipelines::bb_dedisperser constructor: must have pulse_width_ms > 0.0");
    if (ini_params.nt_in <= 0)
	throw runtime_error("rf_pipelines::bb_dedisperser constructor: must have nt_in > 0");
    
    // Initialize this->name.
    stringstream ss;
    ss << "bb_dedisperser(dm_start=" << ini_params.dm_start
       << ",dm_end=" << ini_params.dm_end
       << ",dm_tol=" << ini_params.dm_tol
       << ",pulse_width_ms=" << ini_params.pulse_width_ms << ")";

    this->name = ss.str();
}


virtual void bb_dedispserser::_bind_transform(Json::Value &json_attrs) override
{
    if (!json_attrs.isMember("freq_lo_MHz") || !json_attrs.isMember("freq_hi_MHz"))
	throw runtime_error("bonsai_dedisperser: expected json_attrs to contain members 'freq_lo_MHz' and 'freq_hi_MHz'");

    if (!json_attrs.isMember("dt_sample"))
	throw runtime_error("bonsai_dedisperser: expected json_attrs to contain member 'dt_sample'");
    
    this->freq_lo_MHz = json_attrs["freq_lo_MHz"].asDouble();
    this->freq_hi_MHz = json_attrs["freq_hi_MHz"].asDouble();
    this->dt_sample = json_attrs["dt_sample"].asDouble();
}


void bb_dedisperser::_allocate()
{
    rf_assert(this->plan == nullptr);

    int device_idx = 0;    // FIXME should determine this somehow
    dedisp_error error;
    
    error = dedisp_set_device(device_idx);
    if (error != DEDISP_NO_ERROR)
	throw runtime_error("ERROR: Could not set GPU device: " + dedisp_errmsg(error));

    // Channel width in MHz
    xxxx;  // POSITIVE OR NEGATIVE??
    double df = float(this->freq_hi_MHz - this->freq_lo_MHz) / this->nfreq;
    
    error = dedisp_create_plan(&plan, this->nfreq, this->dt_sample, this->freq_hi_MHz, df);
    if (error != DEDISP_NO_ERROR)
	throw runtime_error("ERROR: Could not create dedispersion plan: " + dedisp_errmsg(error));

    error = dedisp_generate_dm_list(plan, ini_params.dm_start, ini_params.dm_end, ini_params.pulse_width_ms, ini_params.dm_tol);
    if (error != DEDISP_NO_ERROR)
	throw runtime_error("ERROR: Failed to generate dm list: " + dedisp_errmsg(error));

    if (verbosity >= 1)
	cout << "bb_dedisperser: number of trial DM's = " << dedisp_get_dm_count(plan);

    this->ndm = dedisp_get_dm_count(plan);
    this->max_delay = dedisp_get_max_delay(plan);
    this->nt_out = ini_params.
	
    this->in8 = make_utpr<uint8_t> (nt_in * nfreq);
    this->out32 = make_uptr<float> (ndm * nt_out);
}


void bb_dedisperser::_start_pipeline(Json::Value &json_attrs)
{
    this->runflag = false;
}
	

void bb_dedisperser::_process_chunk(float *intensity, ssize_t istride, float *weights, ssize_t wstride, ssize_t pos)
{
    rf_assert(this->plan != nullptr);
    
    int nt_valid = xx;

    // Check weights

    for (int ifreq = 0; ifreq < nfreq; ifreq++) {
	for (int it = 0; it < nt_valid; it++)
	    rf_assert(weights[ifreq*wstride + it] == 1.0);
	for (int it = nt_valid; it < nt_chunk; it++)
	    rf_assert(weights[ifreq*wstride + it] == 0.0);
    }

    // Quantize and copy
    
    for (int ifreq = 0; ifreq < nfreq; ifreq++)
	for (int it = 0; it < nt_valid; it++)
	    in8[(pos+it)*nfreq + ifreq] = bytetoquant(10. * intensity[ifreq*istride+it]);

    if (runflag || xx)
	return;
    
    error = dedisp_execute(plan, nsamps,
			   input, in_nbits,
			   (dedisp_byte *)output, out_nbits,
			   DEDISP_USE_DEFAULT);
    if( error != DEDISP_NO_ERROR ) {
	printf("\nERROR: Failed to execute dedispersion plan: %s\n",
	       dedisp_get_error_string(error));
	return -1;
    }
    
    this->runflag = true;
    
    error = dedisp_execute( this->plan, this->nsamps, input, in_nbits, (dedisp_byte *)output, out_nbits, DEDISP_USE_DEFAULT);

    if( error != DEDISP_NO_ERROR ) 
	{
		cout<< "\nERROR: Failed to execute dedispersion plan" <<endl;
    }

    //cout << "Execution done" <<endl;

    for(dedisp_size nd = 0; nd < local_dm_count; nd++) 
	{
		for (dedisp_size ns = 0; ns < local_nsamps_computed; ns++)
	   	{
			dedisp_size idx = nd*local_nsamps_computed + ns;
	       	dedisp_float val = output[idx];
		
			if(val > maximum ) 
			{
		     	this->max_trigger = val;
		     	this->max_trigger_dm = dmlist[nd];
		     	this->max_trigger_tfinal = ns * this->dt;
		     	maximum = val;
			}
		}
    }
    
    //cout << this->max_trigger << endl;
    //cout << this->max_trigger_dm << endl;
    //cout << this->max_trigger_tfinal << endl;
	int a = (dedisp_byte)240;
	//cout << a << endl;
    int b = (dedisp_byte)2400;
	//cout << b << endl;
	//cout << "Exiting process_chunk" << endl;
    
	free(output);
    free(input);
}


void bb_dedisperser::_end_substream(Json::Value &json_outputs)
{
    if (!runflag)
	xxx;
    
    // Placeholder values
    double sn = 30.0;
    //double dm = 100.0;
    //double tarr = 2.0;

    // Initialize fields in JSON output.
    //this->json_per_substream["frb_sn"] = sn;
    json_outputs["frb_global_max_trigger_dm"] = this->max_trigger_dm;
    json_outputs["frb_global_max_trigger_tfinal"] = this->max_trigger_tfinal;
    json_outputs["frb_global_max_trigger"] = this->max_trigger;
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


std::shared_ptr<wi_transform> make_bb_dedisperser(const bb_dedisperser_initializer &ini_params);
{
    return std::make_shared<bb_dedisperser> (ini_params);
}


}  // namespace rf_pipelines
