//
// rf_pipelines: plugin-based radio astronomy pipelines.
//
// Note: This code is best "documented by example", so if you're seeing it for the first
// time, I recommend starting with the example programs in the examples/ directory!
//
// Warning: currently, the C++ and python API's are pretty well documented, so if you're
// sticking to one or the other then the code should look reasonable.  However, if
// you want to python-wrap a class written in C++, then you'll find that the code in 
// rf_pipelines/rf_pipelines_c.cpp is kind of a mess!  I hope to improve this soon.
// In the meantime, if you want to python-wrap a C++ class, just email me and I'll help
// navigate the mess!
//
// There are two fundamental classes, wi_streams ("weighted intensity stream")
// and wi_transforms ("weighted intensity transform").  These are both virtual
// base classes, with various subclasses defined.  For example, there are
// subclasses of wi_stream to simulate data, read captured data from files,
// or capture real-time data from the network.
//
// A wi_stream incrementally generates chunks of intensity data which is channelized,
// meaning that the data lives in a (freq channel, time sample) 2D array, and weighted,
// meaning that there is a parallel array of floating-point weights.
//
// A wi_transform operates on weighted intensity data, modifying it in place.  For
// example, a detrender or an RFI removal stage could be a wi_transform.  There are
// also pseudo-transforms which process the data without actually modifying it, for
// example plotters or dedispersers.  Transforms run in sequence: each transform's
// output is the input to the next transform.
//
// A generally important thing to be aware of!  Throughout rf_pipelines, we always
// use a frequency channel ordering which goes from highest frequency to lowest.
// This is the same ordering used in the CHIME data and in the bonsai code.
// If you're writing a stream or transform class which "naturally" uses the 
// opposite ordering, then you'll need to do some indexology to translate when
// reading and writing from the rf_pipelines arrays.
//
// Here is a list of the streams and transforms which are currently available in C++.
// More transforms are available in python (see rf_pipelines module docstring for details).
// The python transforms are intended for rapid prototyping and should eventually be 
// translated to C++ for speed.
//
// Factory functions which return streams (std::shared_ptr<wi_stream>):
//
//   make_psrfits_stream(f)        -> reads data from psrfits file (e.g. GBNCC data)
//   make_chime_stream()           -> reads data from file in CHIME hdf5 format
//   make_gaussian_noise_stream()  -> outputs gaussian random data
//
// Factory functions which return transforms (std::shared_ptr<wi_transform>):
//
//   make_simple_detreneder()
//      -> a really boneheaded detrending transform which just subtracts the mean in chunks
//   make_chime_file_writer()
//      -> writes stream to a single file in CHIME hdf5 format
//   make_bonsai_dedisperser() 
//      -> runs data through the bonsai dedisperser (dedispersion output is written to an hdf5 file)
//

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


namespace constants {
    //
    // Throughout the CHIMEFRB backend, we represent times in seconds, but the raw packets use timestamps
    // constructed from FPGA counts.  We convert by assuming that each FPGA count is exactly 2.56e-6 seconds.
    // The precise conversion matters (to machine precision!) when predicting the location of the noise
    // source edges from the timestamps.  Therefore, the 2.56e-6 "magic number" must be used consistently
    // throughout rf_pipelines and ch_vdif_assembler.
    //
    static constexpr double chime_seconds_per_fpga_count = 2.56e-6;
};


// -------------------------------------------------------------------------------------------------
//
// Factory functions returning wi_streams


// PSRFITS file stream (e.g. gbncc)
extern std::shared_ptr<wi_stream> make_psrfits_stream(const std::string &filename);


//
// CHIME file streams, either from single file, explciit filename, or acquisition directory.
// In the 'acqusition directory' case, the directory is scanned for filenames of the form NNNNNNNN.h5, where N=[0,9].
//    
// The 'nt_chunk' arg is the chunk size used internally when moving data from hdf5 file
// into the rf_pipelines buffer.  If unspecified or zero, it will default to a reasonable value.
//
// If 'noise_source_align' is nonzero, then it should be equal to the DETRENDER chunk size (not the chime_file_stream nt_chunk).
// In this case, the stream will align the noise source edges with the detrender chunks, by discarding initial data if necessary.
//
// Note: a quick way to inspect a CHIME hdf5 file is using the 'ch-show-intensity-file' and 'ch-plot-intensity-file'
// programs, in the ch_frb_io github repo.
//
extern std::shared_ptr<wi_stream> make_chime_stream_from_acqdir(const std::string &filename, ssize_t nt_chunk=0, ssize_t noise_source_align=0);
extern std::shared_ptr<wi_stream> make_chime_stream_from_filename(const std::string &filename, ssize_t nt_chunk=0, ssize_t noise_source_align=0);
extern std::shared_ptr<wi_stream> make_chime_stream_from_filename_list(const std::vector<std::string> &filename_list, ssize_t nt_chunk=0, ssize_t noise_source_align=0);


// Simple stream which simulates Gaussian random noise
//
//   nfreq            Number of frequency channels
//   nt_tot           Total number of time samples written before stream ends.
//   freq_lo_MHz      Lowest frequency in band (e.g. 400 for CHIME)
//   freq_hi_MHz      Highest frequency in band (e.g. 800 for CHIME)
//   dt_sample        Length of a time sample in seconds
//   nt_chunk         Stream block size (if zero, will default to a reasonable value)
//
extern std::shared_ptr<wi_stream> make_gaussian_noise_stream(ssize_t nfreq, ssize_t nt_tot, double freq_lo_MHz, double freq_hi_MHz, double dt_sample, double sample_rms=1.0, ssize_t nt_chunk=0);


// -------------------------------------------------------------------------------------------------
//
// Factory functions returning wi_transforms


//
// simple_detrender: this the simplest possible detrending algorithm.  We really
// need something better here!  It just divides the data into chunk, and subtracts
// the time-average of the data for every (chunk, frequency_channel) pair.
//
// The 'nt_detrend' constructor argument is the detrending chunk size (in number of samples).
//
extern std::shared_ptr<wi_transform> make_simple_detrender(ssize_t nt_detrend);


//
// This is a pseudo-transform which doesn't actually modify the data, it just writes it to a file in
// CHIME hdf5 format.  (For now, the entire stream is written to a single file, I'll generalize later
// to break the stream into multiple files.)
//
// If 'clobber' is false, and the target file already exists, an exception will be thrown rather than clobbering the old file.
// If 'nt_chunk' is set to zero, a default chunk size will be chosen.
//
// The meaning of the 'bitshuffle' arg is:
//   0 = no compression
//   1 = try to compress, but if plugin fails then just write uncompressed data instead
//   2 = try to compress, but if plugin fails then print a warning and write uncompressed data instead
//   3 = compression mandatory
//
// Note: a quick way to inspect a CHIME hdf5 file is using the 'ch-show-intensity-file' and 'ch-plot-intensity-file'
// programs, in the ch_frb_io github repo.
//
std::shared_ptr<wi_transform> make_chime_file_writer(const std::string &filename, bool clobber=false, int bitshuffle=2, ssize_t nt_chunk=0);


//
// Returns a "transform" which doesn't actually modify the data, it just runs the bonsai dedisperser.  
// The output is a stream of coarse-grained triggers which are written to an output hdf5 file.  
// The dedisperser must be initialized from a config hdf5 file produced with the program 
// 'bonsai-mkweight' in the bonsai github repo.
//
// Note that the program 'bonsai-plot-triggers.py' in the bonsai github repo may be useful
// for quick visual inspection of the bonsai output.
//
// If 'nt_per_file' is zero, then all triggers will be written to a single "monster file".  
// Otherwise multiple files will be written.  Note that nt_per_file is the number of input 
// time samples (the number of coarse-grained triggers is usually much smaller).
//
// The 'ibeam' argument determines the assignment of threads to cores and can probably
// be zero except in special situations.
//
// FIXME 1: Currently the only trigger "processing" which can be done is writing the triggers
// to an hdf5 file for later analysis.
//
// FIXME 2: Currently the dedisperser must be initialized from a config hdf5 file (rather than
// the simpler config text file) since we use analytic weights to normalize the triggers.
// Since the analytic weights are only correct for unit-variance noise, the trigger normalization
// will be wrong for a real experiment, and the triggers won't be meaningfully normalized to
// "sigmas".  All of this is just a placeholder until Monte Carlo trigger variance estimation
// is implemented in bonsai.
//
extern std::shared_ptr<wi_transform> make_bonsai_dedisperser(const std::string &config_hdf5_filename, const std::string &output_hdf5_filename, int nt_per_file=0, int ibeam=0);


// -------------------------------------------------------------------------------------------------
//
// The 'wi_stream' and 'wi_transform' virtual base classes.
//
// These define the API that you'll need to implement, in order to make new steams and transforms.


// Note: for a reference example showing how to implement a wi_stream, check out gaussian_noise_stream.cpp
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
    ssize_t nfreq;            // number of frequency channels
    double freq_lo_MHz;       // lowest frequency in band (e.g. 400 for CHIME)
    double freq_hi_MHz;       // highest frequency in band (e.g. 400 for CHIME)
    double dt_sample;         // length of a sample in seconds
    ssize_t nt_maxwrite;      // block size of stream (defined as max number of time samples per call to setup_write())

    wi_stream() :
	nfreq(0), freq_lo_MHz(0.), freq_hi_MHz(0.), dt_sample(0.), nt_maxwrite(0)
    { }

    virtual ~wi_stream() { }
    
    //
    // This is the virtual function which must be implemented to define a stream.
    //
    // The 'run_state' argument is an object containing ring buffers which the stream should
    // write its data to.  For the definition of 'class wi_run_state' see below.  The stream_body()
    // function will probably consist of a loop which moves blocks of data from some source (a file
    // or network stream) into the ring buffers.
    //
    // "Moving" a block of data is done in two steps.  First, call wi_run_state::setup_write() to request 
    // space in the ring buffers.  This will return bare pointers to chunks of memory inside the ring buffers.
    // (There are two pointers, one for the 'intensity' array and one for the 'weights'.)  After filling
    // these memory areas with data, call wi_run_state::finalize_write() to advance the ring buffers.
    //
    // The stream can also define multiple "substreams" by calling wi_run_state::start_substream() and
    // wi_run_state::end_substream().  The downstream transforms should reset state between substreams.
    // At the moment this feature isn't very well-supported, so it's probably best for all streams to
    // represent their data as a single substream.
    //
    // Summarizing, wi_stream::stream_body() should look something like this.
    // See 'class wi_run_state' below for more details on setup_write(), finalize_write(), etc!
    //
    //   for (...) {                          // outer loop over substreams
    //      run_state.start_substream();
    //      for (...) {                       // inner loop over data blocks
    //          run_state.setup_write();
    //          run_state.finalize_write();
    //      }
    //      run_state.end_substream()
    //   }
    //
    virtual void stream_body(wi_run_state &run_state) = 0;

    // This non-virtual function isn't defined by the wi_stream subclass, it's called 
    // "from the outside" to run the rf_pipeline.
    void run(const std::vector<std::shared_ptr<wi_transform> > &transforms, bool noisy=true);
};


//
// Note: for a reference example showing how to implement a wi_stream, check out simple_detrender.cpp.
// (This may be too simple to be an ideal example, I might suggest something different later!)
//
struct wi_transform {

    // The following members must be initialized by the subclass.  The initialization may be
    // done either in the subclass constructor, or in the member function wi_stream::set_stream()
    // below.  The latter option may be more convenient if it's useful to inspect the stream
    // parameters.  For example, the transform may wish to set 'nfreq' to the value of stream.nfreq.

    ssize_t nfreq = 0;        // number of frequency channels which the transform expects to receive from the stream
    ssize_t nt_chunk = 0;     // chunk size for process_chunk(), see below
    ssize_t nt_prepad = 0;    // prepad size for process_chunk(), see below
    ssize_t nt_postpad = 0;   // postpad size for process_chunk(), see below
    
    wi_transform() { }

    virtual ~wi_transform() { }

    // The subclass must define the four virtual functions which follow.

    //
    // set_stream(): this is called once, at the beginning of a pipeline run.
    // Note that the following members of 'wi_stream' may be useful:
    //
    //     stream.nfreq         number of frequency channels in stream
    //     stream.freq_lo_MHz   lower limit of frequency band
    //     stream.freq_hi_MHz   upper limit of frequency band
    //     stream.dt_sample     sample length in seconds
    //     stream.nt_maxwrite   internal chunk size of stream, this probably won't be useful.
    //
    virtual void set_stream(const wi_stream &stream) = 0;

    //
    // start_substream(): Each stream can divide its output into one or more "substreams" which 
    // are bracketed with start_substream() and end_substream() calls.  Currently I don't use the 
    // "substreaming"  feature much: all streams just define a single substream, and not all 
    // transforms support multiple substreams.
    //
    // However, I anticipate it being a useful feature when we implement real-time network streams,
    // since we'll want a way to finalize state when the correlator goes down (or repoints) and
    // restart when it comes back.
    //
    // The 'isubstream' arg is 0 for the first substream, 1 for the second substream, etc.
    // The 't0' arg is the initial time of the substream in seconds, relative to an arbitrary 
    // stream-defined origin.
    //
    virtual void start_substream(int isubstream, double t0) = 0;

    //
    // process_chunk(): This routine is called to deliver chunks of data to the transform.  Each 
    // transform defines three buffer sizes (see above): nt_chunk, nt_prepad, and nt_postpad.
    //
    // Each call to process_chunk() is responsible for processing one "chunk" of data with 2D shape
    // (nfreq, nt_chunk).  The 'intensity' and 'weights' arguments are floating-point arrays with
    // this logical shape.
    //
    // The 'stride' argument is the memory stride between logically adjacent frequencies in the
    // 'intensity' and 'weights' 2D arrays.  In other words, the memory location of the array
    // element with indices (ifreq,it) is intensity[ifreq*stride + it] and similarly for the weights.
    //  
    // Important: the frequency axis of the arrays is ordered from highest frequency to lowest!
    // This is the same ordering used by 'bonsai', and in both GBNCC and CHIME data, so this ordering
    // seemed most convenient.
    //
    // The 't0' and 't1' args are timestamps in seconds at the beginning and end of the chunk.
    // (To be precise, t0 is the start time of the first sample in the chunk, and t1 is the end
    // time of the last sample in the chunk.)  This information is mostly redundant since the
    // initial time of the substream is passed in start_substream(), and the time sample length
    // is available in set_stream().  However, I'm anticipating that for long-running streams it 
    // may be useful to allow for a small amount of timestamp "drift" via the t0/t1 args.
    //
    // Some transforms will need to do read-only inspection of data outside the chunk.  For example,
    // a detrending transform may need to look at values of the data a little bit before and after
    // the chunk, in order to detrend data in the chunk.  
    //
    // A transform which needs to "see into the future" can set nt_postpad to a nonzero value.  
    // In this case, it will be safe to index the 'intensity' and 'weights' with time indices in
    // the range 0 <= it < (nt_chunk + nt_postpad).  Important: the transform is only allowed to 
    // modify the first 'nt_chunk' samples!
    //
    // A transform which needs to "see into the past" can set self.nt_prepad to a nonzero value.
    // In this case, the prepadding data is passed as separate arrays, via the 'pp_intensity' and
    // 'pp_weights' args.  These are 2D arrays with logical shapes (nfreq, nt_prepad).  Note that
    // these arrays have their own memory stride 'pp_stride' which will not be equal to 'stride'.
    // Note also that prepadding works differently from the postpadding case, where the extra data 
    // is passed by extending the chunk array.  If nt_prepad is zero, then the 'pp_intensity' and 
    // 'pp_weights' arguments are NULL.
    //
    // Each transform can initialize its nt_chunk, nt_prepad, and nt_postpad independently of
    // the other transforms, and the rf_pipelines library will handle the necessary buffering 
    // and rechunking.
    //
    // Transforms should make sure to handle the case where there are many zeroes in the 'weights'
    // array, including the extreme case where the weights are all zeros.  One situation where
    // many zeros arise is at the end of a stream, where the total stream length may not be
    // a multiple of nt_chunk, and so the rf_pipelines library will append zero-weight data.
    //
    virtual void process_chunk(double t0, double t1, float *intensity, float *weights, ssize_t stride, float *pp_intensity, float *pp_weights, ssize_t pp_stride) = 0;

    // end_substream(): counterpart to start_substream() above
    virtual void end_substream() = 0;
};


// -------------------------------------------------------------------------------------------------
//
// Low-level classes.


// Helper class for wi_run_state (probably not useful from the outside world)
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


//
// This class contains ring buffers which hold the intensity data and weights as they move
// through the transform chain.  The details are hidden from the wi_transforms, but if you're
// implementing a new wi_stream, you'll need to call the public member functions below.
//
// For more details, see comments in class wi_stream, or see gaussian_noise_stream.cpp
// for a reference example.
// 
class wi_run_state {
public:
    wi_run_state(const wi_stream &stream, const std::vector<std::shared_ptr<wi_transform> > &transforms, bool noisy);

    // stream params
    const ssize_t nfreq;
    const ssize_t nt_stream_maxwrite;

    // The 't0' arg is the substream start time in seconds, relative to an arbitrary stream-defined origin.
    void start_substream(double t0);
    
    //
    // This is called to reserve space in the ring buffers for a block of new data from the stream.
    //
    // The 'nt' arg is the number of time samples to be written (i.e. the intensity and weights arrays
    // will be 2d arrays with shape (nfreq,nt))
    //
    // The 'intensityp' and 'weightp' args will be initialized to bare pointers inside the ring buffers.
    // These are logical arrays of shape (nfreq, nt) but the memory layout is non-contiguous: consecutive
    // time indices are adjacent in memory, but consecutive frequencies are separated by a large offset.
    // The 'stride' argument will be initialized to the value of this offset.  In other words, the 
    // memory location of the intensity array element with (frequency,time) indices (ifreq,it) is
    //
    //    intensityp + ifreq*stride + it
    //
    // and likewise for the weights.
    //
    // If the 'zero_flag' argument is set, the intensity and weights arrays will be zeroed before returning
    // pointers to the caller.
    //
    // The 't0' argument is the start time of the block in seconds.  This is optional (see below) since it
    // can be inferred from the substream start time, the value of wi_stream::dt_sample, and the number of
    // samples written so far.  However, it may be useful to specify t0 occasionally in order to keep track
    // of slow timestamp drifts over time.  For example in the chimefrb pipeline, the intensity samples 
    // always correspond to a fixed number of FPGA counts, and the fpga clock drifts on long timescales.
    // 
    void setup_write(ssize_t nt, float* &intensityp, float* &weightp, ssize_t &stride, bool zero_flag, double t0);

    // A version of setup_write() which infers the value of t0 from the substream start time, the value
    // of wi_stream::dt_sample, and the number of samples written so far, assuming ideal regularly spaced
    // samples with no drift.
    void setup_write(ssize_t nt, float* &intensityp, float* &weightp, ssize_t &stride, bool zero_flag);
    
    // The 'nt' arg is the number of samples written (currently this must be the same as the corresponding
    // argument to setup_write(), but this will be generalized so that finalize_write() can write fewer
    // samples than originally requested).
    void finalize_write(ssize_t nt);

    void end_substream();

protected:
    friend void wi_stream::run(const std::vector<std::shared_ptr<wi_transform> > &transforms, bool noisy);

    // make noncopyable
    wi_run_state(const wi_run_state &) = delete;
    wi_run_state& operator=(const wi_run_state &) = delete;

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
};


}  // namespace rf_pipelines

#endif // _RF_PIPELINES_HPP
