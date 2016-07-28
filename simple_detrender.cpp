#include "rf_pipelines_internals.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif


//
// simple_detrender: this the simplest possible detrending algorithm.  We really
// need something better here!  It just divides the data into chunk, and subtracts
// the time-average of the data for every (chunk, frequency_channel) pair.
//

struct simple_detrender : public wi_transform {
    // The 'nt_chunk' constructor argument is the detrending chunk size (in number of samples).
    simple_detrender(ssize_t nt_detrend);

    // Note that we define start_substream() and end_substream() to be empty functions.
    // Since the simple_detrender doesn't maintain state between chunks, we don't need to
    // know when substreams begin and end.

    virtual void set_stream(const wi_stream &stream);
    virtual void start_substream(int isubstream, double t0) { }
    virtual void end_substream() { }

    virtual void process_chunk(double t0, double t1, float *intensity, float *weights, ssize_t stride, float *pp_intensity, float *pp_weights, ssize_t pp_stride);
};


simple_detrender::simple_detrender(ssize_t nt_detrend)
{
    if (nt_detrend <= 1)
	throw runtime_error("simple_detrender: nt_detrend must be > 1");

    // By initializing 'nt_chunk' to nt_detrend, rf_pipelines will deliver the data
    // in chunks of size nt_detrend.  This is convenient since the simple_detrender
    // can just process each chunk independently.
    this->nt_chunk = nt_detrend;

    // Note that 'nt_prepad' and nt_postpad are automatically initialized to zero, which
    // is what we want.  We still need to initialize 'nfreq'.  This will be done in set_stream().
}


void simple_detrender::set_stream(const wi_stream &stream)
{
    // this->nfreq is "the number of frequency channels which the transform expects".  We just
    // initialize this to the actual number of channels in the stream.  Most transforms will do
    // the same thing, but in a few cases (e.g. bonsai) nfreq is predetermined, so instead we do:
    //   assert(stream.nfreq == this->nfreq)
    // to get a runtime check.

    this->nfreq = stream.nfreq;
}


void simple_detrender::process_chunk(double t0, double t1, float *intensity, float *weights, ssize_t stride, float *pp_intensity, float *pp_weights, ssize_t pp_stride)
{
    //
    // In process_chunk(), all we do is compute the time-average of the intensities independently
    // for each frequency channel and subtract it.  A few minor things to note:
    // 
    //  - we use a weighted average, which in particular means that intensities which are masked
    //    (i.e. weight zero) aren't used
    //  - make sure to handle the corner case where all weights are zero
    //  - note the use of 'stride' when indexing the intensity and weights arrays
    //

    for (ssize_t ifreq = 0; ifreq < this->nfreq; ifreq++) {
	float num = 0.0;  // sum of weighted intensities
	float den = 0.0;  // sum of weights

	for (ssize_t it = 0; it < this->nt_chunk; it++) {
	    num += weights[ifreq*stride + it] * intensity[ifreq*stride + it];
	    den += weights[ifreq*stride + it];
	}
	
	// Avoid division by zero (which can happen if a chunk is entirely masked, i.e. weight zero)
	if (den <= 0.0)
	    continue;

	// Note that (num/den) is the weighted mean of all the intensity samples.
	for (ssize_t it = 0; it < this->nt_chunk; it++)
	    intensity[ifreq*stride + it] -= (num/den);
    }
}

// Externally visible factory function declared in rf_transforms.hpp
shared_ptr<wi_transform> make_simple_detrender(ssize_t nt_chunk)
{
    return make_shared<simple_detrender> (nt_chunk);
}


}  // namespace rf_pipelines
