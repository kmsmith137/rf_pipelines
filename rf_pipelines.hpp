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
//   make_chime_stream_from_acqdir()
//   make_chime_stream_from_filename()
//   make_chime_stream_from_filename_list()
//   make_gaussian_noise_stream()
//   make_psrfits_stream()
//
// Factory functions which return transforms (std::shared_ptr<wi_transform>):
//
//   make_bonsai_dedisperser()     Runs data through bonsai dedisperser.
//   make_chime_file_writer()      Writes stream to a single file in CHIME hdf5 format.
//   make_chime_packetizer()       Converts a stream to UDP packets, and sends them over the network.
//   make_intensity_clipper()      "Clips" an array by masking outlier intensities.
//   make_polynomial_detrender()   Detrends along time or frequenciy axis, by subtracting a best-fit polynomial.
//   make_std_dev_clipper()        "Clips" an array by masking rows/cols whose standard deviation is an outlier.
//
// See below for more info on all these functions!


#ifndef _RF_PIPELINES_HPP
#define _RF_PIPELINES_HPP

#if (__cplusplus < 201103) && !defined(__GXX_EXPERIMENTAL_CXX0X__)
#error "This source file needs to be compiled with C++0x support (g++ -std=c++0x)"
#endif

#include <set>
#include <vector>
#include <memory>
#include <iostream>
#include <json/json.h>

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif

struct wi_stream;
struct wi_transform;
class wi_run_state;
struct outdir_manager;   // declared in rf_pipelines_internals.hpp
struct plot_group;       // declared in rf_pipelines_internals.hpp


namespace constants {

    // Throughout the CHIMEFRB backend, we represent times in seconds, but the raw packets use timestamps
    // constructed from FPGA counts.  We convert by assuming that each FPGA count is exactly 2.56e-6 seconds.
    // The precise conversion matters (to machine precision!) when predicting the location of the noise
    // source edges from the timestamps.  Therefore, the 2.56e-6 "magic number" must be used consistently
    // throughout rf_pipelines and ch_vdif_assembler.

    static constexpr double chime_seconds_per_fpga_count = 2.56e-6;

    // Some of our assembly language transforms have parameters, such as max allowed downsampling,
    // which must be known at compile time.  If you need to change them, note the following:
    //
    //   - make sure that the downsampling-related parameters are still powers of two
    //   - the only downside is increased compile time and librf_pipelines.so file size
    //   - you will need to recompile rf_pipelines after the change

    static constexpr int polynomial_detrender_max_degree = 20;
    static constexpr int max_frequency_downsampling = 32;
    static constexpr int max_time_downsampling = 32;

    // Number of single-precision floats which fit into a SIMD word on this machine.
    // For now, we just hardcode 8, assuming a CPU with the AVX instruction set but not AVX-512,
    // but we define a compile-time constant here as a hook for future improvement.

    static constexpr int single_precision_simd_length = 8;
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


//
// CHIME network stream.  Receives UDP packets in "CHIME L0-L1 format".
//
// This interface is less general than the low-level interface in ch_frb_io: 
// only one beam can be received, and not all boolean options are supported.
//
// If the 'udp_port' argument is zero, then the default chimefrb port will be used.
//
extern std::shared_ptr<wi_stream> make_chime_network_stream(int udp_port=0, int beam_id=0);


// Simple stream which simulates Gaussian random noise
//
//   nfreq               Number of frequency channels
//   nt_tot              Total number of time samples written before stream ends.
//   freq_lo_MHz         Lowest frequency in band (e.g. 400 for CHIME)
//   freq_hi_MHz         Highest frequency in band (e.g. 800 for CHIME)
//   dt_sample           Length of a time sample in seconds
//   sample_rms          RMS of intensity samples (Gaussian distributed)
//   nt_chunk            Stream block size (if zero, will default to a reasonable value)
//   randomize_weights   If true, weights will be uniform random numbers (if false, all weights will be 1.0)
//
extern std::shared_ptr<wi_stream> make_gaussian_noise_stream(ssize_t nfreq, ssize_t nt_tot, double freq_lo_MHz, double freq_hi_MHz, 
							     double dt_sample, double sample_rms=1.0, ssize_t nt_chunk=0, bool randomize_weights=false);


// -------------------------------------------------------------------------------------------------
//
// Factory functions returning wi_transforms


enum axis_type {
    AXIS_FREQ = 0,
    AXIS_TIME = 1,
    AXIS_NONE = 2
};

// In misc.cpp
extern std::ostream &operator<<(std::ostream &os, axis_type axis);

//
// polynomial_detrender: detrends along either the time or frequency axis,
// by subtracting a best-fit polynomial.  The detrending is independent in
// every "row" (where "row" means "frequency channel" in the case of time-axis
// detrending, or "time sample" in the case of frequency-axis detrending).
//
// If the fit is poorly conditioned then the entire row will be masked (by
// setting its weights to zero).  The threshold is controlled by the parameter
// 'epsilon'.  I think that 1.0e-2 is a reasonable default here, but haven't
// experimented systematically.
//
extern std::shared_ptr<wi_transform> make_polynomial_detrender(int nt_chunk, axis_type axis, int polydeg, double epsilon=1.0e-2);


// A "simple detrender" is a time-axis polynomial fitter with degree zero.
// FIXME: this will be removed soon, in favor of calling make_polynomial_detrender() directly.
inline std::shared_ptr<wi_transform> make_simple_detrender(ssize_t nt_detrend)
{
    std::cerr << "make_simple_detrender(): this function is now deprecated in favor of make_polynomial_detrender()\n";
    return make_polynomial_detrender(nt_detrend, AXIS_TIME, 0);   // polydeg=0
}


//
// intensity_clipper: this "clips" an array by masking outlier intensities.
// The masking is performed by setting elements of the weights array to zero.
//
// The 'sigma' argument is the threshold (in sigmas from the mean) for clipping.  Note
// that the weights are used when calculating both the mean and rms intensity.
//
// The (Df,Dt) args are downsampling factors on the frequency/time axes.
// If no downsampling is desired, set Df=Dt=1.
//
// The 'axis' argument has the following meaning:
//   axis=AXIS_FREQ   clip along frequency axis, with an outer loop over time samples
//   axis=AXIS_TIME   clip along time axis, with an outer loop over frequency samples
//   axis=AXIS_NONE   2-d clipper
//
// If niter > 1, then the mean/rms intensity will be computed using iterated clipping,
// with threshold 'iter_sigma'.  If the 'iter_sigma' argument is zero, then it defaults
// to 'sigma', but the two thresholds need not be the same.
//
// If the 'two_pass' flag is set, a more numerically stable but slightly slower algorithm will be used.
//

extern std::shared_ptr<wi_transform> make_intensity_clipper(int nt_chunk, axis_type axis, double sigma, int niter=1, 
							    double iter_sigma=0.0, int Df=1, int Dt=1, bool two_pass=false);


//
// std_dev_clipper: this "clips" an array by masking rows/columns whose standard deviation is an outlier.
//
// The 'axis' argument has the following meaning:
//   axis=AXIS_FREQ   clip time samples whose variance in frequency is high
//   axis=AXIS_TIME   clip frequency channels whose variance in time is high
//
// The (Df,Dt) args are downsampling factors on the frequency/time axes.
// If no downsampling is desired, set Df=Dt=1.
//
// The 'sigma' argument is the threshold (in sigmas from the mean) for clipping.
//
// If the 'two_pass' flag is set, a more numerically stable but slightly slower algorithm will be used.
//
std::shared_ptr<wi_transform> make_std_dev_clipper(int nt_chunk, axis_type axis, double sigma, int Df=1, int Dt=1, bool two_pass=false);


// Standalone functions with the equivalent functionality to the polynomial_detrender,
// intensity_clipper, and std_dev_clipper transforms.  (See comments above for documentation.)
//
// Note that we generally assume the intensity/weights arrays have the same stride, but
// this isn't really necessary in the kernels, and this assumption could be relaxed if
// necessary.


extern void apply_polynomial_detrender(float *intensity, float *weights, int nfreq, int nt, 
				       int stride, axis_type axis, int polydeg, double epsilon);

extern void apply_intensity_clipper(const float *intensity, float *weights, int nfreq, int nt, 
				    int stride, axis_type axis, double sigma, int niter=1, 
				    double iter_sigma=0.0, int Df=1, int Dt=1, bool two_pass=false);

extern void apply_std_dev_clipper(const float *intensity, float *weights, int nfreq, int nt, int stride,
				  axis_type axis, double sigma, int Df=1, int Dt=1, bool two_pass=false);


// Helper routines for the RFI transforms above, factored out as standalone functions.
//
// wi_downsample(): downsamples an (intensity, weights) pair.  The downsampling factors (Df,Dt)
// must be powers of two.  Note that the normalization of the downsampled weights array differs
// (by a factor of Df*Dt) from the python version of wi_downsample().
//
// weighted_mean_and_rms(): computes weighted mean/rms of a 2D intensity array.
// If the 'niter' argument is >1, then the calculation will be iterated, clipping
// outlier samples which differ from the mean by the specified number of "sigmas".
// If the 'two_pass' flag is set, a more numerically stable but slightly slower algorithm will be used.


extern void wi_downsample(float *out_intensity, float *out_weights, int out_stride, const float *in_intensity, 
			  const float *in_weights, int in_nfreq, int in_nt, int in_stride, int Df, int Dt);


extern void weighted_mean_and_rms(float &mean, float &rms, const float *intensity, const float *weights, 
				  int nfreq, int nt, int stride, int niter=1, double sigma=3.0, bool two_pass=false);


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
// Converts a stream to UDP packets in "CHIME L0_L1" format, and sends them over the network.
// This interface is less general than the low-level interface in ch_frb_io: only one beam can
// be sent, and not all boolean options are supported.
//
// Some artificial restrictions: the stream 'nfreq' value must be a multiple of 1024, and
// the stream 'dt_sample' value must be an integer multiple of 2.56e-6 seconds.  This is because
// the packet protocol doesn't include a count of total frequency channels, or the fpga clock
// rate, so these parameters are frozen to the CHIME instrumental values.
//
// The 'dstname' argument is a string of the form HOSTNAME:PORT.  For example 'localhost:13178' or
// 'chimer.physics.ubc.ca:13178'.  If the port is omitted then the default chimefrb port is used.
// (Be careful sending packets over the internet since the bandwidth can be very high!)
//
// The 'wt_cutoff' argument is used to convert the rf_pipelines 'weights' array to a boolean mask.
// This conversion is necessary because the CHIME L0_L1 packet format doesn't support a floating-point
// weight array.  Samples with weight below the cutoff will be masked.
//
// If the 'target_gbps' argument is nonzero, then output will be "throttled" to the target bandwidth, specified
// in Gbps.  If target_gbps=0, then packets will be sent as quickly as possible.
//
// The nfreq_coarse_per_packet, nt_per_packet arguments define the amount of data sent per packet.
// The nt_per_chunk arg just determines an internal chunk size and isn't very important (must be
// a multiple of nt_per_packet; suggest a value like 512).
//
extern std::shared_ptr<wi_transform> make_chime_packetizer(const std::string &dstname, int nfreq_coarse_per_packet, int nt_per_chunk, 
							   int nt_per_packet, float wt_cutoff, double target_gbps);


//
// Returns a "transform" which doesn't actually modify the data, it just runs the bonsai dedisperser.  
// The dedisperser must be initialized from a config hdf5 file produced with the program 
// 'bonsai-mkweight' in the bonsai github repo.
//
// If 'trigger_hdf5_filename' is a nonempty string, then triggers will be written to one
// or more HDF5 output files.  If 'nt_per_file' is zero, then all triggers will be written
// to a single "monster file".  Otherwise multiple files will be written.  Note that nt_per_file
// is the number of input time samples (the number of coarse-grained triggers is usually
// much smaller).
//
// If 'trigger_plot_stem' is a nonempty string, then realtime trigger plots will be written.  
// In this case, the nt_per_file arg must be positive.  Filenames are of the form
//   ${trigger_plot_stem}_${plot_number}_tree${tree_index}.png.
//
// FIXME: Currently the dedisperser must be initialized from a config hdf5 file (rather than
// the simpler config text file) since we use analytic weights to normalize the triggers.
// Since the analytic weights are only correct for unit-variance noise, the trigger normalization
// will be wrong for a real experiment, and the triggers won't be meaningfully normalized to
// "sigmas".  All of this is just a placeholder until Monte Carlo trigger variance estimation
// is implemented in bonsai.
//
extern std::shared_ptr<wi_transform> make_bonsai_dedisperser(const std::string &config_hdf5_filename, const std::string &trigger_hdf5_filename, 
							     const std::string &trigger_plot_stem, int nt_per_file=0);


// Some day, this factory function will return a C++ implementation of the 'badchannel_mask' class.
// Right now, it is a placeholder which throws an exception if called.

extern std::shared_ptr<wi_transform> make_badchannel_mask(const std::string &maskpath, int nt_chunk=1024);


// -------------------------------------------------------------------------------------------------
//
// The 'wi_stream' and 'wi_transform' virtual base classes.
//
// These define the API that you'll need to implement, in order to make new steams and transforms.


// Note: for a reference example showing how to implement a wi_stream, check out gaussian_noise_stream.cpp
struct wi_stream {
    //
    // The fields { nfreq, ... , nt_maxwrite } must be initialized by the base class, either
    // in its constructor or in stream_start().
    //
    // As a detail, it's sometimes convenient to defer initialization of some fields until stream_start(),
    // since they may not be known until the stream starts running.  An example is a network stream, which
    // may not get the value of dt_sample (say) until the first packet is received.
    //
    // Note: don't set nt_maxwrite to an excessively large value, since there is an internal
    // buffer of approximate size (24 bytes) * nfreq * nt_maxwrite.
    //
    ssize_t nfreq = 0;              // number of frequency channels
    double freq_lo_MHz = 0.0;       // lowest frequency in band (e.g. 400 for CHIME)
    double freq_hi_MHz = 0.0;       // highest frequency in band (e.g. 800 for CHIME)
    double dt_sample = 0.0;         // length of a sample in seconds
    ssize_t nt_maxwrite = 0;        // block size of stream (defined as max number of time samples per call to setup_write())

    wi_stream() { }

    virtual ~wi_stream() { }
    
    //
    // There are two virtual functions which can be implemented to define a stream.
    //
    // First, a function stream_start() which is called to start the stream.  More precisely, it is called
    // after wi_stream::run(), but before any of the calls to wi_transform::set_stream().  This function
    // serves as a "last chance" to set the stream parameters { nfreq, ..., nt_maxwrite } as discussed above.
    //

    virtual void stream_start() { }   // note non pure virtual, defaults to empty function

    //
    // Second, a pure virtual function stream_body() which generates incremental chunks of data.
    //
    // The 'run_state' argument is an object containing ring buffers which the stream should
    // write its data to.  For the definition of 'class wi_run_state' see below.  The stream_body()
    // function will consist of a loop which moves blocks of data from some source (a file
    // or network connection) into the rf_pipelines ring buffer.
    //
    // "Moving" a block of data is done in two steps.  First, call wi_run_state::setup_write() to request 
    // space in the ring buffers.  This will return bare pointers to chunks of memory inside the ring buffers.
    // (There are two pointers, one for the 'intensity' array and one for the 'weights'.)  After filling
    // these memory areas with data, call wi_run_state::finalize_write() to advance the ring buffers.
    //
    // The stream can also define multiple "substreams" by calling wi_run_state::start_substream() and
    // wi_run_state::end_substream().  The downstream transforms should reset state between substreams.
    // At the moment the "multiple-substream" feature isn't very well-supported, so it's probably best 
    // for all streams to represent their data as a single substream.
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

    //
    // This non-virtual function runs the rf_pipeline.
    //
    // 'outdir' is the rf_pipelines output directory, where the rf_pipelines json file will
    // be written, in addition to other transform-specific output files such as plots. 
    //
    // If 'outdir' is an empty string, then the json file will not be written, and 
    // any transform which tries to write an output file (such as a plotter_transform)
    // will throw an exception.
    //
    // If the 'json_output' pointer is non-null, then the pipeline's json output is
    // written there (i.e. same data which is written to rf_pipelines.json)
    //
    // If 'clobber' is false, then an exception will be thrown if the pipeline tries to
    // overwrite an old rf_pipelines.json file.
    //
    // The meaning of the 'verbosity' argument is:
    //   0 = no output
    //   1 = high-level summary output (names of transforms, number of samples processed etc.)
    //   2 = show all output files
    //   3 = debug trace through pipeline
    //
    void run(const std::vector<std::shared_ptr<wi_transform> > &transforms, 
	     const std::string &outdir = ".", 
	     Json::Value *json_output = nullptr,
	     int verbosity=2, bool clobber=true);
};


struct wi_transform {
    // Subclass should initialize the transform name in its constructor.
    // The transform name will appear in the json output, and in python __str__().
    std::string name;

    // The following members must be initialized by the subclass.  The initialization may be
    // done either in the subclass constructor, or in the member function wi_stream::set_stream()
    // below.  The latter option may be more convenient if it's useful to inspect the stream
    // parameters.  For example, the transform may wish to set 'nfreq' to the value of stream.nfreq.

    ssize_t nfreq = 0;        // number of frequency channels which the transform expects to receive from the stream
    ssize_t nt_chunk = 0;     // chunk size for process_chunk(), see below
    ssize_t nt_prepad = 0;    // prepad size for process_chunk(), see below
    ssize_t nt_postpad = 0;   // postpad size for process_chunk(), see below

    //
    // Each transform can define key/value pairs which get written to the pipeline json output file.
    // This data is always written on a per-substream basis, but it's convenient not to reinitialize it
    // for every substream.  Therefore we define three json objects which differ in when they get cleared.
    //
    // Example syntax for modifying these (from a member function of wi_transform):
    //    this->json_persistent["key"] = value;
    //    this->json_per_stream["key"]["key2"].append(value);
    //
    // FIXME should these be Json::Object instead of Json::Value?
    //
    // FIXME there is currently no way to modify these from python.
    // One solution would be to python-wrap the Json objects but that sounds a little painful, and I'm not
    // sure it's really necessary.  Another possibility is just to have parallel json_* members in the
    // Python class py_wi_substream, and add a (non-pure) virtual functions to wi_transform to retrieve
    // and clear json.  In the python subclass these could make appropriate calls into the python interpreter.
    // I think this would work!
    //
    // The transforms's per-substream json output is obtained by combining
    //   - the json_persistent, json_per_stream, and json_per_substream objects
    //   - data from the plot_groups (see below), which goes into "plots"
    //
    Json::Value json_persistent;       // never cleared
    Json::Value json_per_stream;       // cleared just before start_stream()
    Json::Value json_per_substream;    // cleared just before start_substream()

    //
    // These helper functions are used by wi_transforms which write output files (e.g. hdf5, png).
    //
    // make_plot_group(): Each transform's output plots are divided into one or more "plot groups".
    //   For example, the bonsai dedisperser can write one plot group per internally defined tree.
    //   The 'nt_per_pix' arg is the number of pipeline time samples per x-pixel in the plot.
    //   The 'ny' arg is the number of y-pixels (assumed to be the same for all plots in the group).
    //   The return value is the group_id arg needed in add_plot(), and group_ids always go 0,1,...
    //
    // add_plot(): Call just before writing a plot.
    //   The range of time samples in the plot is [it0:it0+nt).
    //   The pixel dimensions of the plot are (nx,ny).  These are redundant since they can be deduced
    //     from (it0,nt) but we use them for error checking.
    //   The return value is the full pathname ('basename' with the stream output_dir prepended)
    //
    // add_file(): Call just before writing a non-plot file, to check for filename collisions between transforms.
    //   The return value is the full pathname ('basename' with stream output_dir prepended)
    //
    // The "plots" part of the transform's (per-substream) json output consists of
    //
    //   - a list with one element per plot_group, consisting of
    //
    //   - an object with scalar keys { 'name', 'nt_per_pix', 'ny', 'it0', 'it1', 'files' }.
    //
    //       name = name of the plot_group (e.g. 'waterfall')
    //       nt_per_pix = number of time samples per x-pixel
    //       ny = number of y-pixels in each image
    //       [it0:it1] = sample index range covered by plot_group
    //       files = a list with one element per image file, consisting of
    //
    //   - an object with keys { "filename", "it0", "nx" }
    //
    //       it0 = first sample index in range covered by plot
    //       nx = number of x-pixels in image
    //
    int add_plot_group(const std::string &name, int nt_per_pix, int ny);   // returns group id
    std::string add_plot(const std::string &basename, int64_t it0, int nt, int nx, int ny, int group_id=0);
    std::string add_file(const std::string &basename);


    // Data used internally by rf_pipelines -- probably a bad idea to use these fields directly!
    // Note: the outdir_manager is a nonempty pointer if and only if the transform is currently running.
    std::vector<std::shared_ptr<plot_group> > plot_groups;
    std::shared_ptr<rf_pipelines::outdir_manager> outdir_manager;
    double time_spent_in_transform = 0.0;


    wi_transform() { }

    virtual ~wi_transform() { }


    // --------------- The subclass must define the five pure virtual functions which follow ---------------

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
    wi_run_state(const wi_stream &stream, 
		 const std::vector<std::shared_ptr<wi_transform> > &transforms, 
		 const std::shared_ptr<outdir_manager> &manager, 
		 Json::Value *json_output, int verbosity);

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
    // always correspond to a fixed number of FPGA counts, and the FPGA clock drifts on long timescales.
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
    friend void wi_stream::run(const std::vector<std::shared_ptr<wi_transform> > &transforms, 
			       const std::string &outdir, Json::Value *json_outputs, 
			       int verbosity, bool clobber);

    // make noncopyable
    wi_run_state(const wi_run_state &) = delete;
    wi_run_state& operator=(const wi_run_state &) = delete;

    // outputs
    const std::shared_ptr<outdir_manager> manager;
    Json::Value *json_output;

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
    int verbosity;

    // buffers
    wraparound_buf main_buffer;
    std::vector<wraparound_buf> prepad_buffers;
    
    void output_substream_json();
    void clear_per_substream_data();
};


}  // namespace rf_pipelines

#endif // _RF_PIPELINES_HPP
