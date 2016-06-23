// Note: this code isn't cleaned up and might be hard to understand!

#include "rf_pipelines_internals.hpp"

#ifdef HAVE_PSRFITS
extern "C" {
#include <psrfits.h>
}
#endif // HAVE_PSRFITS


using namespace std;

namespace rf_pipelines {
#if 0
};  // pacify emacs c-mode!
#endif


#ifndef HAVE_PSRFITS

shared_ptr<wi_stream> make_psrfits_stream(const string &filename)
{
    throw runtime_error("make_psrfits_stream() was called, but this rf_pipelines instance was compiled without psrfits");
}

#else  // HAVE_PSRFITS

//
// The psrfits reader is factored as two classes:
//   psrfits_wrapper: lightweight standalone C++ wrapper for 'struct psrfits'.
//   psrfits_stream: wraps around the psrfits_wrapper, and implements the rf_pipeline stream API.
//
// I factored things this way because I thought the psrfits_wrapper might some day be useful
// outside rf_pipelines.
//
// Note: I haven't learned how the psrfits library works yet, and the code below is
// brain-dead cut-and-paste from a toy program Scott Ransom sent me!
//


// psrfits_wrapper: a lightweight standalone C++ wrapper for 'struct psrfits' from Scott Ransom C library
struct psrfits_wrapper {
    struct psrfits pf;
    string filename;

    //
    // Note: when the psrfits_wrapper object is created, it reads the first row of the FITS file into this->data.
    // Each call to read_next_row() will read the next row into this->data, setting this->eof if necessary.
    //
    // The data is represented as a shape-(nt_per_row, nfreq) unsigned 8-bit integer array.
    // Note that time is the slowest varying index.  This is the opposite of the rf_pipelines 
    // convention, so there is a transpose in psrfits_stream::get_next_chunk() below.
    //
    unsigned char *data;    // bare pointer to a buffer which lives in the 'struct psrfits'
    bool eof;

    // The psrfits file format supports a per-frequency floating-point weight.
    // Note that there is no way to get a per-(frequency,time_sample) weight.
    float *freq_weights;     // bare pointer to a buffer which lives in the 'struct psrfits'

    int nfreq;
    int npol;
    int nt_per_row;

    double freq_lo_MHz;
    double freq_hi_MHz;
    double dt_sample;     // in seconds

    // make noncopyable
    psrfits_wrapper(const psrfits_wrapper &) = delete;
    psrfits_wrapper &operator=(const psrfits_wrapper &) = delete;

    explicit psrfits_wrapper(const string &filename_)
	: filename(filename_)
    {
	memset(&pf, 0, sizeof(pf));   // This just seems like a good idea
	this->eof = false;

	// FIXME not sure if strdup() is needed here; does psrfits hang on to the pointer?
	// (Also, there is currently a small memory leak here.)
	char *f = strdup(filename.c_str());

	char *filenames[2] = { f, nullptr };
	psrfits_set_files(&pf, 1, filenames);

	pf.tot_rows = 0;
	pf.N = 0;
	pf.T = 0;
	pf.status = 0;

	int rv = psrfits_open(&pf);
	if (rv) { 
	    // FIXME print error to string and include text in exception, rather than printing to stderr
	    fits_report_error(stderr, rv);
	    throw runtime_error(filename + ": file not found or wrong format");
	}

	// I think this case should be OK, but this warning seems like a good idea...
	if (pf.hdr.nbits != 8)
	    cerr << (filename + ": warning: bit depth is not equal to 8, this is a case that has never been tested so there may be bugs!\n");

	this->nfreq = pf.hdr.nchan;
	this->npol = pf.hdr.npol;

	rf_assert(nfreq >= 2);
	rf_assert(npol == 1);  // FIXME not sure how to handle the case where npol != 1

	// Allocate buffers
	pf.sub.dat_freqs = aligned_alloc<float> (nfreq);
	pf.sub.dat_weights = aligned_alloc<float> (nfreq);
	pf.sub.dat_offsets = aligned_alloc<float> (nfreq * npol);
	pf.sub.dat_scales = aligned_alloc<float> (nfreq * npol);
	
	int bits_per_subint = pf.sub.bytes_per_subint * 8;
	if (bits_per_subint % pf.hdr.nbits)
	    throw runtime_error(filename + ": bits_per_subint is not divisible by hdr.nbits, not quite sure what do in this case");

	this->data = aligned_alloc<unsigned char> (bits_per_subint / pf.hdr.nbits);
	pf.sub.data = data;
	pf.sub.rawdata = data;

	this->freq_weights = pf.sub.dat_weights;

	rv = psrfits_read_subint(&pf);
	if (rv)
	    throw runtime_error(filename + ": couldn't read first row of FITS file");

	// Not sure whether it's necessary to call psrfits_read_subint() before making the assignments below

	this->nt_per_row = pf.hdr.nsblk;
	this->dt_sample = pf.hdr.dt;

	float freq0 = pf.sub.dat_freqs[nfreq-1];
	float freq1 = pf.sub.dat_freqs[0];

	// If this ever fails, it should be trivial to fix, but I don't want to add code to handle this case
	// until I have an actual file to work with for debugging.
	if (freq0 > freq1)
	    throw runtime_error(filename + ": psrfits wrapper code currently assumes frequencies are ordered from highest to lowest");

	this->freq_lo_MHz = freq0 - (0.5/(nfreq-1)) * (freq1-freq0);
	this->freq_hi_MHz = freq1 + (0.5/(nfreq-1)) * (freq1-freq0);
    }

    // Reads next FITS row into this->data array (as 8-bit unsigned chars).  Sets this->eof if row couldn't be read.
    void read_next_row()
    {
	if (this->eof)
	    return;

	cerr << ".";
	int rv = psrfits_read_subint(&pf);
	this->eof = (rv != 0);  // FIXME how to differentiate errors from end-of-file?
    }

    ~psrfits_wrapper()
    {
	free(pf.sub.dat_freqs);
	free(pf.sub.dat_weights);
	free(pf.sub.dat_offsets);
	free(pf.sub.dat_scales);
	free(data);
	
	pf.sub.dat_freqs = pf.sub.dat_weights = pf.sub.dat_offsets = pf.sub.dat_scales = nullptr;
	data = pf.sub.data = pf.sub.rawdata = nullptr;
    }
};


// psrfits_stream: wraps around the psrfits_wrapper class, and implements the rf_pipeline stream API.
struct psrfits_stream : public wi_stream
{
    shared_ptr<psrfits_wrapper> p;

    psrfits_stream(const shared_ptr<psrfits_wrapper> &p);

    virtual ~psrfits_stream() { }
    virtual void stream_body(wi_run_state &run_state);
};


psrfits_stream::psrfits_stream(const shared_ptr<psrfits_wrapper> &p_)
{ 
    this->p = p_;
    this->nfreq = p->nfreq;
    this->freq_lo_MHz = p->freq_lo_MHz;
    this->freq_hi_MHz = p->freq_hi_MHz;
    this->dt_sample = p->dt_sample;
    this->nt_maxwrite = p->nt_per_row;
}


void psrfits_stream::stream_body(wi_run_state &run_state)
{
    // FIXME psrfits_stream currently sets the initial time of the stream to zero.
    // What's a sensible way to determine an initial time from a 'struct psrfits'?
    double t0 = 0.0;
    run_state.start_substream(t0);

    while (!p->eof) {
	float *intensity;
	float *weights;
	int stride;
	bool zero_flag = false;
	run_state.setup_write(this->nt_maxwrite, intensity, weights, stride, zero_flag);

	// Transpose and convert uint8 -> float
	for (int it = 0; it < nt_maxwrite; it++)
	    for (int ifreq = 0; ifreq < nfreq; ifreq++)
		intensity[ifreq*stride + it] = (float)p->data[it*nfreq + ifreq];

	// psrfits weights are per-(frequency,chunk), not per-(frequency,sample)
	for (int ifreq = 0; ifreq < nfreq; ifreq++) {
	    float w = p->freq_weights[ifreq];
	    if (w < 0.0)
		throw runtime_error(p->filename + ": negative weight in file, this is currently treated as an error");
	
	    for (int it = 0; it < nt_maxwrite; it++)
		weights[ifreq*stride + it] = w;
	}

	p->read_next_row();
    }

    run_state.end_substream();
}


// Factory function returning new stream
shared_ptr<wi_stream> make_psrfits_stream(const string &filename)
{
    shared_ptr<psrfits_wrapper> p = make_shared<psrfits_wrapper> (filename);
    return make_shared<psrfits_stream> (p);
}


#endif  // HAVE_PSRFITS


}   // namespace rf_pipelines
