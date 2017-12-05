#include "rf_pipelines_internals.hpp"
#include "dedisp.h"

using namespace std;

namespace rf_pipelines {
#if 0
};  // pacify emacs c-mode
#endif

// The 'bb_dedisperser' wrapper class converts Ben Barsdell's 'dedisp' GPU code
// to the rf_pipelines API.  It does this by subclassing rf_pipelines::wi_transform,
// and defining the approporiate virtual functions.

/* byte2quant is a utility function which coverts input to unsigned range */
dedisp_float tounsigned(dedisp_float f)
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

dedisp_byte bytetoquant(dedisp_float f)
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

struct bb_dedisperser : public wi_transform {
    // Constructor arguments
    dedisp_float dm_start;
    dedisp_float dm_end;
    dedisp_float dm_tol;
    dedisp_float pulse_width_ms;

    dedisp_plan plan;

    float max_trigger;
    float max_trigger_dm;
    float max_trigger_tfinal;

    int device_idx;

    dedisp_size  nchans;
    dedisp_float dt; 
    dedisp_float f0;
    dedisp_float df;
    dedisp_size nsamps;
    dedisp_error error;
    // Additional class members (for example, the dedisp_plan object) can be added here.

    bb_dedisperser(double dm_start_, double dm_end_, double dm_tol_, double pulse_width_ms_);

    // The behavior of bb_dedisperser (or any subclass of wi_transform) is defined by
    // implementing these virtual functions.  For documentation on what needs to be
    // implemented, see comments in rf_pipelines.hpp.

    virtual void set_stream(const wi_stream &stream) override;
    virtual void start_substream(int isubstream, double t0) override;
    
    virtual void process_chunk(double t0, double t1, float *intensity, float *weights, ssize_t stride, float *pp_intensity, float *pp_weights, ssize_t pp_stride) override;
    virtual void end_substream() override;
};


bb_dedisperser::bb_dedisperser(double dm_start_, double dm_end_, double dm_tol_, double pulse_width_ms_) :
    dm_start(dm_start_),
    dm_end(dm_end_),
    dm_tol(dm_tol_),
    pulse_width_ms(pulse_width_ms_)
{ 
    // Initialize this->name
    //cout << "Inside bb_dedisperser constructor "<< endl;
    stringstream ss;
    ss << "bb_dedisperser(dm_start=" << dm_start << ",dm_end=" << dm_end
       << ",dm_tol=" << dm_tol << ",pulse_width_ms=" << pulse_width_ms << ")";

    this->name = ss.str();
	
    this->max_trigger = 0.0;
    this->max_trigger_dm = 0.0;
    this->max_trigger_tfinal = 0.0;
    this->device_idx = 0;

    this->nchans = 0;
    this->dt = 0;
    this->f0 = 0;
    this->df = 0;
    this->nsamps = 0;

    //cout << "Done work inside constructor. Exiting!"<< endl;

}

// Placeholders for virtual functions follow.

void bb_dedisperser::set_stream(const wi_stream &stream)
{
    //throw std::runtime_error("bb_dedisperser::set_stream() called but not implemented yet");

    //cout << "Inside set_stream" << endl;
    //cout << "\nnt_maxwrite is : " << stream.nt_maxwrite << endl;
    //cout << "\nnt_chunks is : " << stream.nt_chunk << endl;
    //cout << "\nnfreq is : " << stream.nfreq << endl;
	//cout << "\ndt_sample : " << stream.dt_sample << endl;
	//cout << "\nnsamps : " << stream.nt_tot << endl;

	//Original definition
    this->nfreq = stream.nfreq;
    this->nt_chunk = 65536;
    this->nt_prepad = 0;
    this->nt_postpad = 0;

    this->nchans = stream.nfreq;
    this->dt = stream.dt_sample; // assuming no downsampling
    this->f0 = stream.freq_hi_MHz;
    this->df = float(stream.freq_hi_MHz - stream.freq_lo_MHz)/ this->nchans;
    this->nsamps = 65536;
    
    error = dedisp_set_device(this->device_idx);
    if(error != DEDISP_NO_ERROR )
   	{
		cout << "\nERROR: Could not set GPU device : %s\n " << dedisp_get_error_string(error) << endl;	
    }


    //cout << "Device set " <<endl;


    error = dedisp_create_plan(&plan, this->nchans, this->dt, this->f0, this->df);
    if( error != DEDISP_NO_ERROR) 
	{
		cout << "\nERROR: Could not create dedispersion plan" << endl;
    }

    //cout << "Plan created " << endl;

	//cout << "DM start is "<<dm_start << endl;
    dedisp_float *dm_list = (dedisp_float*) malloc( 800 * 4);
    for (int i = 0; i < 800; i++)
    {
        dm_list[i] = i + 100;

    }
    error = dedisp_set_dm_list(plan, dm_list, 800);
    if( error != DEDISP_NO_ERROR) 
	{
	    cout << "\nERROR: Failed to generate dm list" << endl;
    }

    //cout << "DM list generated!" << endl;
    //cout << "Exiting set_stream()" << endl;
}

void bb_dedisperser::start_substream(int isubstream, double t0)
{
    
	//cout << "Inside start_substream()! " << endl;
    if( isubstream > 0) 
	{
        throw std::runtime_error("bb_dedisperser::start_substream() called with isubstream > 0");
    }

    //cout << "Exiting start_substream()! " << endl;
}

void bb_dedisperser::process_chunk(double t0, double t1, float *intensity, float *weights, ssize_t stride, float *pp_intensity, float *pp_weights, ssize_t pp_stride)
{

    //cout << "Inside process_chunk" << endl;
    
	dedisp_size local_dm_count = dedisp_get_dm_count(this->plan);
    dedisp_size local_max_delay = dedisp_get_max_delay(this->plan);
    dedisp_size local_nsamps_computed = dedisp_size(this->nsamps- local_max_delay);
    //dedisp_size local_nsamps = dedisp_size((t0 - t1)/this->dt);
    
	const dedisp_float *dmlist;
    dmlist = dedisp_get_dm_list(this->plan);
    
	//cout << this->nsamps << endl;
    //cout << local_dm_count << endl;
	//cout << local_nsamps_computed << endl;
   
    int in_nbits = 8;
	int out_nbits = 32;	
	dedisp_float *output = (dedisp_float*)malloc((this->nsamps - local_max_delay) * local_dm_count * (out_nbits/8));
    dedisp_byte *input = (dedisp_byte*)malloc(this->nsamps * this->nchans * (in_nbits/8));
    
    float maximum = 0.0;
	float minima = 0.0;
	float maxima = 0.0;
    if (input == NULL)
   	{
		cout << "\nERROR: Failed to allocate output array \n" << endl;
    }

	//cout << "nchans is : " << this->nchans << endl;
	//cout << "\n nt_chunk is : "<< this->nt_chunk << endl;
	//cout << "\n stride is : "<< stride << endl;
    
    for(dedisp_size ns = 0; ns < this->nsamps; ns++) 
	{
		for(dedisp_size nc = 0; nc < this->nchans; nc++) 
		{
			// Inverted to time major order for dedisp usage
			//cout << stride << endl;	
	       	if (minima > intensity[ns+(stride*nc)])
			{
				minima = intensity[ns+(stride*nc)];
			}
			
			if (maxima < intensity[ns+(stride*nc)])
			{
				maxima = intensity[ns+(stride*nc)];
			}

			input[(ns*this->nchans)+nc] = bytetoquant(intensity[ns+(stride*nc)]*10);
			
			//cout << input[(ns*this->nchans)+nc] << endl;
			//cout << intensity[ns+(stride*nc)] << endl;
		}
    }	
	//cout << maxima << endl;
	//cout <<  minima << endl;
    //cout << "Inversion Done" << endl;

    
    if (output == NULL) 
	{
		cout << "\nERROR: Failed to allocate output array \n" << endl;
    }

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

void bb_dedisperser::end_substream()
{
    //cout << "Inside end_substream" << endl;
    dedisp_destroy_plan(this->plan);
    

    // After the dedisperser has run to completion, you should have values of the FRB DM,
    // arrival time, and signal-to-noise.  In general, rf_pipelines transforms communicate
    // their outputs by setting fields in their JSON output (this will end up in a pipeline
    // output file named something like 'rf_pipelines_0.json').  The syntax below shows
    // how to do this!

    // Placeholder values
    double sn = 30.0;
    //double dm = 100.0;
    //double tarr = 2.0;

    // Initialize fields in JSON output.
    //this->json_per_substream["frb_sn"] = sn;
    this->json_per_substream["frb_global_max_trigger_dm"] = this->max_trigger_dm;
    this->json_per_substream["frb_global_max_trigger_tfinal"] = this->max_trigger_tfinal;
	this->json_per_substream["frb_global_max_trigger"] = this->max_trigger;
}


// -------------------------------------------------------------------------------------------------
//
// make_bb_dedisperser()
//
// This interface (function which returns a shared pointer to the wi_transform base class) 
// allows the implementation of struct bb_dedisperser to be "hidden" from the rest of rf_pipelines.


std::shared_ptr<wi_transform> make_bb_dedisperser(double dm_start, double dm_end, double dm_tol, double pulse_width_ms)
{
    return std::make_shared<bb_dedisperser> (dm_start, dm_end, dm_tol, pulse_width_ms);
}


}  // namespace rf_pipelines
