#include <random>
#include "rf_pipelines_internals.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
};  // pacify emacs c-mode!
#endif


// This is a simple stream which outputs an independent Gaussian random number
// for every (frequency, time) pair.  This code is also intended to be a reference 
// for implementing new wi_streams!


class gaussian_noise_stream : public wi_stream
{
protected:
    ssize_t nt_tot;
    double sample_rms;

public:
    //
    // The subclass constructor should initialize the following members,
    // which are inherited from the wi_stream base class:
    //
    //   nfreq            Number of frequency channels
    //   freq_lo_MHz      Lowest frequency in band (e.g. 400 for CHIME)
    //   freq_hi_MHz      Highest frequency in band (e.g. 800 for CHIME)
    //   dt_sample        Length of a time sample in seconds
    //   nt_maxwrite      Stream block size (=max number of time samples per call to setup_write())
    //
    // As an alternative to initializing these fields in the constructor, another possibility is to
    // define the function wi_stream::stream_start (overriding an empty virtual function in the base
    // class) and initialize them there.  This function gets called when the stream starts running.
    //
    // An example where stream_start() is useful is a network stream, where the value of dt_sample (say)
    // may not be known until the stream actually receives packets.  By deferring initialization of dt_sample 
    // to start_stream(), the constructor can return immediately (start_stream() would need to block until
    // the first packet is received but that's OK).
    //
    gaussian_noise_stream(ssize_t nfreq_, ssize_t nt_tot_, double freq_lo_MHz_, double freq_hi_MHz_, double dt_sample_, double sample_rms_, ssize_t nt_chunk)
    {
	if (nt_chunk == 0)
	    nt_chunk = 1024;   // default
	
	this->nfreq = nfreq_;
	this->freq_lo_MHz = freq_lo_MHz_;
	this->freq_hi_MHz = freq_hi_MHz_;
	this->dt_sample = dt_sample_;
	this->nt_maxwrite = nt_chunk;

	this->sample_rms = sample_rms_;
	this->nt_tot = nt_tot_;

	// Some sanity checking of arguments.
	// A detail: we don't sanity check the base class members { nfreq, ..., nt_maxwrite } since
	// they're automatically sanity-checked elsewhere.

	rf_assert(sample_rms >= 0.0);
	rf_assert(nt_tot > 0);
    }

    virtual ~gaussian_noise_stream() { }

    //
    // This overrides the pure virtual function wi_stream::stream_body() and defines the stream.
    // For a high-level overview, see comments in rf_pipelines.hpp (in class wi_stream).
    // The 'run_state' argument contains ring buffers which the stream will write data into.
    //
    virtual void stream_body(wi_run_state &run_state) override
    {
	std::random_device rd;
	std::mt19937 rng(rd());
	std::normal_distribution<float> dist(0, sample_rms);

	// In general a stream can be composed of multiple "substreams" (see rf_pipelines.hpp).
	// Here we put everything into a single stream, with nominal starting time t=0.
	run_state.start_substream(0.0);
	
	// Current position in stream (such that 0 <= it0 < nt_tot, where nt_tot is the stream length)
	ssize_t it0 = 0;

	while (it0 < nt_tot) {
	    // Number of samples to write in this block
	    ssize_t nt = min(nt_maxwrite, nt_tot-it0);
	    
	    // Call wi_run_state::setup_write() to reserve space in the ring buffers.
	    // For details, see comments in rf_pipelines.hpp (in class wi_run_state).
	    float *intensity;
	    float *weights;
	    ssize_t stride;
	    bool zero_flag = false;   // no need to zero buffers, since we'll overwrite them shortly
	    run_state.setup_write(nt, intensity, weights, stride, zero_flag);

	    // After setup_write() returns, the 'intensity' and 'weights' pointers point to memory
	    // regions in the ring buffers.  These are logical 2D arrays of shape (nfreq,nt), laid
	    // out in memory so that time samples are adjacent, but frequencies are separated by
	    // offset 'stride'.  Therefore, the memory location of the intensity array element with
	    // (frequency,time) indices (ifreq,it) is intensity[ifreq*stride+it], and similarly for
	    // the weights.

	    // Fill the intensity array with Gaussian random numbers, and initialize the weights to 1.
	    for (ssize_t ifreq = 0; ifreq < nfreq; ifreq++) {
		for (ssize_t it = 0; it < nt; it++) {
		    intensity[ifreq*stride + it] = dist(rng);
		    weights[ifreq*stride + it] = 1.0;
		}
	    }

	    // Call wi_run_state::finalize_write() after filling the arrays, to advance the ring buffers.
	    // After finalize_write() returns, the intensity and weights pointers are no longer valid.
	    run_state.finalize_write(nt);

	    it0 += nt;
	}

	run_state.end_substream();
    }
};


//
// Factory function which returns a shared pointer to the wi_stream.  (Using a factory function
// means that details of the base class 'gaussian_noise_stream' and others can be hidden in this 
// source file to avoid overpopulating rf_pipelines.hpp with definitions of many classes, but
// this is just a preference!)
//
shared_ptr<wi_stream> make_gaussian_noise_stream(ssize_t nfreq, ssize_t nt_tot, double freq_lo_MHz, double freq_hi_MHz, double dt_sample, double sample_rms, ssize_t nt_chunk)
{
    return make_shared<gaussian_noise_stream> (nfreq, nt_tot, freq_lo_MHz, freq_hi_MHz, dt_sample, sample_rms, nt_chunk);
}


}   // namespace rf_pipelines
