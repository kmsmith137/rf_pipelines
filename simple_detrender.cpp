#include "rf_pipelines_internals.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif


// simple_detrender: just divides the data into chunks and subtracts the mean in each chunk


struct simple_detrender : public wi_transform {
    simple_detrender(int nt_chunk_);

    virtual void set_stream(const wi_stream &stream);
    virtual void start_substream(double t0) { }
    virtual void end_substream() { cerr << "simple_detrender::end_substream() called!\n"; }

    virtual void process_chunk(double t0, double t1, float *intensity, float *weight, int stride, float *pp_intensity, float *pp_weight, int pp_stride);
};


simple_detrender::simple_detrender(int nt_chunk_)
{
    // Note that 'nt_prepad' and nt_postpad are automatically initialized to zero.
    // We still need to initialize 'nfreq'.  This will be done in set_stream().
    this->nt_chunk = nt_chunk_;

    if (nt_chunk <= 1)
	throw runtime_error("simple_detrender: nt_chunk must be > 1");
}


void simple_detrender::set_stream(const wi_stream &stream)
{
    this->nfreq = stream.nfreq;
}


void simple_detrender::process_chunk(double t0, double t1, float *intensity, float *weight, int stride, float *pp_intensity, float *pp_weight, int pp_stride)
{
    for (int ifreq = 0; ifreq < this->nfreq; ifreq++) {
	float num = 0.0;  // sum of weighted intensities
	float den = 0.0;  // sum of weights

	for (int it = 0; it < this->nt_chunk; it++) {
	    num += weight[ifreq*stride + it] * intensity[ifreq*stride + it];
	    den += weight[ifreq*stride + it];
	}
	
	// Avoid division by zero (which can happen if a chunk is entirely masked, i.e. weight zero)
	if (den <= 0.0)
	    continue;

	// Note that (num/den) is the weighted mean of all the intensity samples.
	for (int it = 0; it < this->nt_chunk; it++)
	    intensity[ifreq*stride + it] -= (num/den);
    }
}

// Externally visible factory function
shared_ptr<wi_transform> make_simple_detrender(int nt_chunk)
{
    return make_shared<simple_detrender> (nt_chunk);
}


}  // namespace rf_pipelines
