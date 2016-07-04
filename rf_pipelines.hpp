#ifndef _RF_PIPELINES_HPP
#define _RF_PIPELINES_HPP

#if (__cplusplus < 201103) && !defined(__GXX_EXPERIMENTAL_CXX0X__)
#error "This source file needs to be compiled with C++0x support (g++ -std=c++0x)"
#endif

#include <iostream>
#include <vector>
#include <memory>

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif

struct wi_stream;
struct wi_transform;
class wi_run_state;


// -------------------------------------------------------------------------------------------------
//
// First, the library of built-in streams and transforms.
//
// The C++ syntax for applying a sequence of transforms to a stream is something like this:
//
//   shared_ptr<wi_stream> stream = make_chime_stream(...);
//
//   vector<shared_ptr<wi_transform> > transform_list;
//   transform_list.push_back(make_simple_detrender(...));
//   transform_list.push_back(make_bonsai_dedispserser(...));  // in rf_pipelines_bonsai.hpp
//
//   stream->run(transform_list);
//


// PSRFITS file stream (e.g. gbncc)
extern std::shared_ptr<wi_stream> make_psrfits_stream(const std::string &filename);


//
// CHIME file streams, either from single file, explciit filename, or acquisition directory.
// In the 'acqusition directory' case, the directory is scanned for filenames of the form NNNNNNNN.h5, where N=[0,9].
//    
// The 'nt_chunk' arg is the chunk size used internally when moving data from hdf5 file
// into the rf_pipelines buffer.  If unspecified or zero, it will default to a reasonable value.
//
// Note: a quick way to inspect a CHIME hdf5 file is using the 'ch-show-intensity-file' program,
// in the ch_frb_io github repo.
//
extern std::shared_ptr<wi_stream> make_chime_stream_from_acqdir(const std::string &filename, ssize_t nt_chunk=0);
extern std::shared_ptr<wi_stream> make_chime_stream_from_filename(const std::string &filename, ssize_t nt_chunk=0);
extern std::shared_ptr<wi_stream> make_chime_stream_from_filename_list(const std::vector<std::string> &filename_list, ssize_t nt_chunk=0);


// Simple stream which simulates Gaussian random noise
// If 'nt_chunk' is unspecified or zero, it will default to a reasonable value.
extern std::shared_ptr<wi_stream> make_gaussian_noise_stream(ssize_t nfreq, ssize_t nt_tot, double freq_lo_MHz, double freq_hi_MHz, double dt_sample, double sample_rms=1.0, ssize_t nt_chunk=0);


// Simplest possible detrender: just divides the data into chunks and subtracts the mean in each chunk
extern std::shared_ptr<wi_transform> make_simple_detrender(ssize_t nt_chunk);

//
// For now, the bonsai_dedisperser must be initialized from an hdf5 config file (created using bonsai-mkweight)
//
// Coarse-grained triggers are "processed" by writing them to an hdf5 output file for later analysis 
// (to disable the output file, use an empty string as filename)
//
// FIXME eventually we'll want to generalize this interface!
//
extern std::shared_ptr<wi_transform> make_bonsai_dedisperser(const std::string &config_hdf5_filename, const std::string &output_hdf5_filename, int ibeam=0);


// -------------------------------------------------------------------------------------------------
//
// Second, the 'wi_stream' and 'wi_transform' virtual base classes.
//
// These define the API that you'll need to implement, in order to make new steams and transforms.


struct wi_stream {
    //
    // The subclass is responsible for initializing the fields { nfreq, ..., nt_maxwrite }
    // in its constructor.  (As a detail, we don't require this initialization to be done 
    // via a base class constructor, since the subclass constructor sometimes needs to do 
    // a little processing first, e.g. opening the first file in a sequence.)
    //
    // Note: don't set nt_maxwrite to an excessively large value, since there is an internal
    // buffer of approximate size (24 bytes) * nfreq * nt_maxwrite.
    //
    ssize_t nfreq;
    double freq_lo_MHz;
    double freq_hi_MHz;
    double dt_sample;     // in seconds
    ssize_t nt_maxwrite;      // max number of time samples per call to setup_write()

    wi_stream() :
	nfreq(0), freq_lo_MHz(0.), freq_hi_MHz(0.), dt_sample(0.), nt_maxwrite(0)
    { }

    virtual ~wi_stream() { }

    // run() is intended to be the main interface for running a sequence of transforms on a stream
    void run(const std::vector<std::shared_ptr<wi_transform> > &transforms, bool noisy=true);
    
    //
    // This is the function which must be implemented to define a stream. 
    // Schematically it should look something like this:
    //
    //   for (...) {   // loop over substreams
    //      run_state.start_substream();
    //      while (cond) {
    //          run_state.setup_write();
    //          run_state.finalize_write();
    //      }
    //      run_state.end_substream()
    //   }
    //
    // See below for the definition of 'class wi_run_state'.
    //
    virtual void stream_body(wi_run_state &run_state) = 0;
};


struct wi_transform {
    //
    // The subclass is responsible for initializing the fields { nfreq, ..., nt_postpad },
    // but the initialization can either be done in the subclass constructor, or in the member 
    // function set_stream() below.  The latter option may be more convenient since the value
    // of wi_stream::nfreq can be used to initialize wi_transform::nfreq.
    //
    ssize_t nfreq;
    ssize_t nt_chunk;
    ssize_t nt_prepad;
    ssize_t nt_postpad;
    
    wi_transform() : 
	nfreq(0), nt_chunk(0), nt_prepad(0), nt_postpad(0)
    { }

    virtual ~wi_transform() { }

    // This is the API which must be implemented to define a transform.
    virtual void set_stream(const wi_stream &stream) = 0;
    virtual void start_substream(int isubstream, double t0) = 0;
    virtual void process_chunk(double t0, double t1, float *intensity, float *weights, ssize_t stride, float *pp_intensity, float *pp_weights, ssize_t pp_stride) = 0;
    virtual void end_substream() = 0;
};


// -------------------------------------------------------------------------------------------------
//
// Third: low-level classes which are probably not needed from the outside world.
//
// Exception: if you're implementing a new wi_stream, then you'll probably want to look at the
// definition of wi_run_state.


struct wraparound_buf {
    // specified at construction
    ssize_t nfreq;
    ssize_t nt_contig;
    ssize_t nt_ring;

    // 2d arrays of shape (nfreq, nt_tot)
    std::vector<float> intensity;
    std::vector<float> weights;
    ssize_t nt_tot;

    ssize_t ipos;

    // Main constructor syntax
    wraparound_buf(ssize_t nfreq, ssize_t nt_contig, ssize_t nt_ring);

    // Alternate syntax: use default constuctor, then call construct()
    wraparound_buf();

    void construct(ssize_t nfreq, ssize_t nt_contig, ssize_t nt_ring);
    void reset();

    void setup_write(ssize_t it0, ssize_t nt, float* &intensityp, float* &weightp, ssize_t &stride);
    void setup_append(ssize_t nt, float* &intensityp, float* &weightp, ssize_t &stride, bool zero_flag);
    void append_zeros(ssize_t nt);

    void finalize_write(ssize_t it0, ssize_t nt);
    void finalize_append(ssize_t nt);

    void _copy(ssize_t it_dst, ssize_t it_src, ssize_t nt);
    void _check_integrity();

    static void run_unit_tests();
};


class wi_run_state {
protected:
    friend void wi_stream::run(const std::vector<std::shared_ptr<wi_transform> > &transforms, bool noisy);

    // make noncopyable
    wi_run_state(const wi_run_state &) = delete;
    wi_run_state& operator=(const wi_run_state &) = delete;

    // stream params
    const ssize_t nfreq;
    const ssize_t nt_stream_maxwrite;

    // transform list
    const int ntransforms;
    const std::vector<std::shared_ptr<wi_transform> > transforms;

    // timeline (times are in seconds, relative to arbitrary stream-defined origin)
    double dt_sample;                  // initialized in constructor
    double substream_start_time;       // initialized in start_substream()
    double stream_curr_time;           // set in every call to setup_write()

    // sample counts
    std::vector<ssize_t> transform_ipos;   // satisfies transform_ipos[0] >= transform_ipos[1] >= ...
    ssize_t stream_ipos;
    
    // state=0: initialized
    // state=1: start_substream() called, but first call to setup_write() hasn't happened yet
    // state=2: setup_write(), matching call to finalize_write() hasn't happened yet
    // state=3: finalize_write() called
    // state=4: end_substream() called
    int state;
    int isubstream;
    ssize_t nt_pending;  // only valid in state 2
    bool noisy;

    // buffers
    wraparound_buf main_buffer;
    std::vector<wraparound_buf> prepad_buffers;

public:
    wi_run_state(const wi_stream &stream, const std::vector<std::shared_ptr<wi_transform> > &transforms, bool noisy);

    // Called by wi_stream::run()
    void start_substream(double t0);
    void setup_write(ssize_t nt, float* &intensityp, float* &weightp, ssize_t &stride, bool zero_flag);
    void setup_write(ssize_t nt, float* &intensityp, float* &weightp, ssize_t &stride, bool zero_flag, double t0);
    void finalize_write(ssize_t nt);
    void end_substream();
};


}  // namespace rf_pipelines

#endif // _RF_PIPELINES_HPP
