#include <algorithm>
#include "rf_pipelines_internals.hpp"

#ifdef HAVE_CH_FRB_IO
#include <ch_frb_io.hpp>
#endif

using namespace std;

namespace rf_pipelines {
#if 0
};  // pacify emacs c-mode!
#endif


#ifndef HAVE_CH_FRB_IO

shared_ptr<wi_transform> make_chime_file_writer(const string &filename, bool clobber, int bitshuffle, ssize_t nt_chunk)
{
    throw runtime_error("rf_pipelines::make_chime_file_writer() was called, but rf_pipelines was compiled without ch_frb_io");
}

#else  // HAVE_CH_FRB_IO


struct chime_file_writer : public wi_transform {
    // Constructor args
    string filename;
    bool clobber;
    int bitshuffle;

    // Stream params (not available until set_stream() gets called)
    double freq_lo_MHz;
    double freq_hi_MHz;
    double dt_sample;

    std::unique_ptr<ch_frb_io::intensity_hdf5_ofile> ofile;
    int ichunk = 0;
    
    // Used to repack 'intensity' and 'weights' into contiguous arrays
    // FIXME these can go away when ch_frb_io::intensity_hdf5_ofile::append_chunk() supports a 'stride' argument
    std::vector<float> intensity_contig_buf;
    std::vector<float> weights_contig_buf;


    chime_file_writer(const string &filename_, bool clobber_, int bitshuffle_, ssize_t nt_chunk_)
	: filename(filename_), clobber(clobber_), bitshuffle(bitshuffle_)
    {
	if (nt_chunk_ == 0)
	    nt_chunk_ = 128;  // default value
	
	// Initialization of nfreq postponed to set_stream()
	this->nt_chunk = nt_chunk_;
	this->nt_prepad = 0;
	this->nt_postpad = 0;
    }
    
    virtual ~chime_file_writer() { }

    virtual void set_stream(const wi_stream &stream) override
    {
	this->name = "chime_file_writer(" + filename + ")";
	this->nfreq = stream.nfreq;
	this->freq_lo_MHz = stream.freq_lo_MHz;
	this->freq_hi_MHz = stream.freq_hi_MHz;
	this->dt_sample = stream.dt_sample;

	this->intensity_contig_buf.resize(nfreq * nt_chunk);
	this->weights_contig_buf.resize(nfreq * nt_chunk);
    }

    virtual void start_substream(int isubstream, double t0) override
    {
	if (isubstream > 0)
	    throw runtime_error("rf_pipelines: multiple substreams are not currently supported in class chime_file_writer");

	if (!clobber && file_exists(filename))
	    throw runtime_error(filename + ": file already exists and clobber=false was specified in the the chime_file_writer constructor");

	// Not really correct but that's OK
	vector<string> pol = { "XX" };

	// Note swapped ordering of freq_hi_MHz and freq_lo_MHz.  This is because rf_pipelines always orders frequency channels from highest to lowest.
	this->ofile = make_unique<ch_frb_io::intensity_hdf5_ofile> (filename, nfreq, pol, freq_hi_MHz, freq_lo_MHz, dt_sample, 0, t0, bitshuffle, nt_chunk);
	this->ichunk = 0;
    }

    virtual void process_chunk(double t0, double t1, float *intensity, float *weights, ssize_t stride, float *pp_intensity, float *pp_weights, ssize_t pp_stride) override
    {
	// Repack to contiguous arrays
	for (int ifreq = 0; ifreq < nfreq; ifreq++) {
	    memcpy(&intensity_contig_buf[0] + ifreq*nt_chunk, intensity + ifreq*stride, nt_chunk * sizeof(float));
	    memcpy(&weights_contig_buf[0] + ifreq*nt_chunk, weights + ifreq*stride, nt_chunk * sizeof(float));
	}

	this->ofile->append_chunk(nt_chunk, &intensity_contig_buf[0], &weights_contig_buf[0], ichunk * nt_chunk, t0);
	this->ichunk++;
    }

    virtual void end_substream() override
    {
	// Resetting this pointer will close file
	this->ofile.reset();
    }
};


// See rf_pipelines.hpp for an explanation of the arguments
shared_ptr<wi_transform> make_chime_file_writer(const string &filename, bool clobber, int bitshuffle, ssize_t nt_chunk)
{
    return make_shared<chime_file_writer> (filename, clobber, bitshuffle, nt_chunk);
}


}  // namespace rf_pipelines

#endif  // HAVE_CH_FRB_IO
