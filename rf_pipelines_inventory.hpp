#ifndef _RF_PIPELINES_INVENTORY_HPP
#define _RF_PIPELINES_INVENTORY_HPP

#if (__cplusplus < 201103) && !defined(__GXX_EXPERIMENTAL_CXX0X__)
#error "This source file needs to be compiled with C++11 support (g++ -std=c++11)"
#endif

// Abstract base classes
// ---------------------
//   pipeline_object
//   chunked_pipeline_object
//   wi_stream
//   wi_transform
//
// Container classes
// -----------------
//   pipeline
//   wi_sub_pipeline
//
// Streams
// -------
//   chime_stream_from_acqdir
//   chime_stream_from_filename
//   chime_stream_from_filename_list
//   chime_frb_stream_from_filename
//   chime_frb_stream_from_filename_list
//   chime_frb_stream_from_glob
//   chime_network_stream
//   gaussian_noise_stream
//
// Detrenders
// ----------
//   spline_detrender
//   polynomial_detrender
//
// Clippers
// --------
//   intensity_clipper
//   std_dev_clipper
//   mask_expander
//
// CHIME-specific
// --------------
//   chime_file_writer
//   chime_packetizer
//   chime_16k_spike_mask
//   chime_16k_derippler
//   chime_16k_stripe_analyzer
//
// Miscellaneous transforms
// ------------------------
//   adversarial_masker (*)
//   badchannel_mask (*)
//   bonsai_dedisperser (**) 
//   frb_injector_transform (*)
//   mask_filler (*) 
//   noise_filler (*)
//   plotter_transform (*)
//   variance_estimator (*)
//
// (*) = python-only
// (**) = the bonsai_dedisperser currently has both python and C++ versions

#include <mutex>

// ring_buffer, pipeline_object, etc.
#include "rf_pipelines_base_classes.hpp"

// enum axis_type
#include <rf_kernels/core.hpp>

// A little hack so that all definitions still compile if optional dependencies are absent.
namespace bonsai { class dedisperser; }

namespace ch_frb_io {
    class intensity_network_stream;
    class output_device_pool;
    class memory_slab_pool;
}

namespace rf_pipelines {
#if 0
}  // emacs pacifier
#endif


// -------------------------------------------------------------------------------------------------
//
// pipeline: This ubiquitous container class is used to chain pipeline_objects together.


class pipeline : public pipeline_object {
public:
    std::vector<std::shared_ptr<pipeline_object>> elements;
    
    explicit pipeline(const std::string &name="");
    explicit pipeline(const std::vector<std::shared_ptr<pipeline_object>> &elements, const std::string &name="");
    explicit pipeline(const std::string &class_name, const std::string &name);   // for subclasses (e.g. wi_sub_pipeline)

    void add(const std::shared_ptr<pipeline_object> &p);
    inline int size() const { return elements.size(); }

    virtual Json::Value jsonize() const override;
    static std::shared_ptr<pipeline> from_json(const Json::Value &x);
    
protected:
    virtual void _bind(ring_buffer_dict &rb_dict, Json::Value &json_attrs) override;
    virtual ssize_t _advance() override;
        
    virtual void _allocate() override;
    virtual void _deallocate() override;
    virtual void _start_pipeline(Json::Value &j) override;
    virtual void _end_pipeline(Json::Value &j) override;
    virtual void _reset() override;
    virtual void _unbind() override;
    virtual void _get_info(Json::Value &j) override;
    virtual void _visit_pipeline(std::function<void(const std::shared_ptr<pipeline_object>&,int)> f, const std::shared_ptr<pipeline_object> &self, int depth) override;
    
    virtual ssize_t get_preferred_chunk_size() override;
};


// -------------------------------------------------------------------------------------------------
//
// wi_sub_pipeline: this more specialized container class is used to run a "sub-pipeline"
// at lower (freqency, time) resolution, then upsample and apply the resulting mask.


class wi_sub_pipeline : public pipeline {
public:
    // The initializer allows a flexible syntax where some fields can be specified (i.e. nonzero)
    // and others unspecified (i.e. zero).  For example:
    //
    //    - If 'nfreq_out' is specified and 'Df' is not, then Df will be set to (nfreq_in / nfreq_out).
    //    - If 'nfreq_out' is unspecified and 'Dt' is specified, then nfreq_out will be set to (nfreq_in / Df).
    //    - If 'nfreq_out' and 'Df' are both specified, then an exception will be raised unless nfreq_in = (nfreq_out * Df)
    //    - If neither 'nfreq_out' nor 'Df' are specified, then an exception will be raised.
    //
    // The parameter pair (nds_out, Dt) behaves similarly.

    struct initializer {
	double w_cutoff = 0.0;
	ssize_t nt_chunk = 0;
	ssize_t nfreq_out = 0;  // number of frequency channels after downsampling to sub-pipeline
	ssize_t nds_out = 0;    // cumulative time downsampling (relative to input data) after downsampling to sub-pipeline
	ssize_t Df = 0;         // frequency downsampling factor (between input pipeline and sub-pipeline)
	ssize_t Dt = 0;         // time downsampling factor (between input pipeline and sub-pipeline)
    };

    wi_sub_pipeline(const std::shared_ptr<pipeline_object> &sub_pipeline, const initializer &ini_params);

    virtual Json::Value jsonize() const override;
    static std::shared_ptr<wi_sub_pipeline> from_json(const Json::Value &x);

    // Constructor arguments (saved for use in jsonize())
    const initializer ini_params;
    const std::shared_ptr<pipeline_object> sub_pipeline;

protected:
    virtual void _bind(ring_buffer_dict &rb_dict, Json::Value &json_attrs) override;
    virtual void _visit_pipeline(std::function<void(const std::shared_ptr<pipeline_object>&,int)> f, const std::shared_ptr<pipeline_object> &self, int depth) override;
    
    virtual ssize_t get_preferred_chunk_size() override;
};


// -------------------------------------------------------------------------------------------------
//
// "Utility" classes: mask_expander, pipeline_fork


// mask_expander
//
// Experimental: expands the RFI mask, in a way which is intended to "fill gaps".
//
// It is assuemd that the caller has saved the weights at a previous point in the pipeline
// (using pipeline_fork, see below).  We use the term "prev_mask" to mean the RFI mask at
// this previous point in the pipeline, and "delta_mask" to mean the set of pixels which
// are currently masked in the pipeline, but were not masked in the prev_mask.
//
// By default, the mask_expander actually expands the delta-mask, but this behavior
// can be modified (see the 'alpha' parameter below).
//
// The expansion is done by computing exponential moving averages of the delta-mask in
// both directions, and masking pixels when both averages are above a threshold.  This
// will be written up in more detail later!
//
// Constructor arguments
// ---------------------
//
// 'axis': currently, only AXIS_FREQ is implemented.  AXIS_TIME is coming soon!
//
// 'prev_wname': pipeline bufname of the saved weights (a string).  Note that in order
//   to save the weights at a previous point in the pipeline, you can use a pipeline_fork
//   whose input_bufname parameter is "WEIGHTS" and whose output_bufname is a string which
//   uniquely identifies the saved weights (e.g. "WEIGHTS_SAVE1").  The 'prev_wname' argument
//   to the mask_expander should be the same as the 'output_bufname' argument of the
//   pipeline_fork.
//
// 'width': the decay width of the exponential moving average.  In the AXIS_FREQ case,
//   this is expressed as a fraction of the frequency band, i.e. width=0.1 means that
//   the characteristic width of the mask_expander is 10% of the full frequency band.
//
// 'threshold': value between 0 and 1 which determines how aggressive the mask_expander is.
//   Low values correspond to more masking.  The numerical value can be roughly interpreted
//   as the fraction of data which must be delta-masked before mask expansion will
//   occur.  For example, if threshold=0.1, then mask expansion will occur in regions of
//   the data where ~10% or more of the pixels are delta-masked.
//
// 'alpha': to explain this parameter, we first note that delta-masked pixels are
//   "sources" for the mask_expander, and unmasked pixels are "sinks".  That is,
//   mask expansion occurs in regions where the number of delta-masked pixels
//   relative to the number of unmasked pixels is above a threshold.
//
//   The alpha paramaeter determines how the prev_mask is handled by the mask_expander.
//   By default (alpha=0), prev-masked pixels are "neutral", i.e. they are neither
//   sources nor sinks for the mask_expander.  
//
//   If 0 < alpha < 1, then prev-masked pixels are sinks for the mask_expander, i.e.
//   they reduce the amount of mask expansion, and the amount of reduction is proportional
//   to alpha.  If alpha=1, then prev_masked pixels are equivalent to unmasked pixels.
//
//   If -1 < alpha < 0, then prev-masked pixels are sources for the mask_expander, i.e.
//   they increase the amount of mask expansion, and the amount of reduction is proportional
//   to (-alpha).  If alpha=-1, then prev_masked pixels are equivalent to delta_masked pixels.

extern std::shared_ptr<pipeline_object> make_mask_expander(rf_kernels::axis_type axis, const std::string &prev_wname, double width, double threshold, double alpha=0.0, ssize_t nt_chunk=0);


// pipeline_fork
//
// Creates one or more new pipeline ring_buffers, by copying existing ring_buffers.
// The 'bufnames' argument should be a list of (input_bufname, output_bufname) pairs.
// Frequently, the input_bufname will be one of the built-in names "INTENSITY" or "WEIGHTS".


extern std::shared_ptr<pipeline_object> make_pipeline_fork(const std::vector<std::pair<std::string,std::string>> &bufnames);


// -------------------------------------------------------------------------------------------------
//
// Detrenders.


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
// Note: the 'axis' argument should be one of
//   rf_kernels::AXIS_FREQ
//   rf_kernels::AXIS_TIME

extern std::shared_ptr<wi_transform>
make_polynomial_detrender(int nt_chunk, rf_kernels::axis_type axis, int polydeg, double epsilon=1.0e-2);


// Experimental: spline_detrender.
// I suspect this will work better than the polynomial_detrender, and it will definitely be faster!
// Currently, the only allowed axis type is rf_kernels::AXIS_FREQ.

extern std::shared_ptr<wi_transform>
make_spline_detrender(int nt_chunk, rf_kernels::axis_type axis, int nbins, double epsilon=3.0e-4);


// -------------------------------------------------------------------------------------------------
//
// Clippers.


// The badchannel_mask transform sets bad freq channels of a weights array to 0.
//
// 'mask_path' is the full path to a mask file that contains affected freq
// intervals, written in rows with the following format: e.g., 420.02,423.03.
// If 'maskpath' is an empty string, then no mask file will be read.
//
// 'mask_ranges' is a list of (freq_lo, freq_hi) pairs, which define additional
// frequency ranges to be masked.


// List of (freq_lo, freq_hi) pairs, for badchannel_mask.
using bc_mask_t = std::vector<std::pair<double,double>>;

extern std::shared_ptr<wi_transform> 
make_badchannel_mask(const std::string &mask_path, 
		     const bc_mask_t &mask_ranges = bc_mask_t());


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
//   axis=rf_kernels::AXIS_FREQ   clip along frequency axis, with an outer loop over time samples
//   axis=rf_kernels::AXIS_TIME   clip along time axis, with an outer loop over frequency samples
//   axis=rf_kernels::AXIS_NONE   2-d clipper
//
// If niter > 1, then the mean/rms intensity will be computed using iterated clipping,
// with threshold 'iter_sigma'.  If the 'iter_sigma' argument is zero, then it defaults
// to 'sigma', but the two thresholds need not be the same.
//
// If the 'two_pass' flag is set, a more numerically stable but slightly slower algorithm will be used.

extern std::shared_ptr<wi_transform>
make_intensity_clipper(int nt_chunk, rf_kernels::axis_type axis, double sigma, int niter=1, 
		       double iter_sigma=0.0, int Df=1, int Dt=1, bool two_pass=false);


// std_dev_clipper: this "clips" an array by masking rows/columns whose standard deviation is an outlier.
//
// The 'axis' argument has the following meaning:
//   axis=rf_kernels::AXIS_FREQ   clip time samples whose variance in frequency is high
//   axis=rf_kernels::AXIS_TIME   clip frequency channels whose variance in time is high
//
// The (Df,Dt) args are downsampling factors on the frequency/time axes.
// If no downsampling is desired, set Df=Dt=1.
//
// The 'sigma' argument is the threshold (in sigmas from the mean) for clipping.
//
// If the 'two_pass' flag is set, a more numerically stable but slightly slower algorithm will be used.

extern std::shared_ptr<wi_transform>
make_std_dev_clipper(int nt_chunk, rf_kernels::axis_type axis, double sigma, int Df=1, int Dt=1, bool two_pass=false);


// -------------------------------------------------------------------------------------------------
//
// bonsai_dedisperser: a "transform" which doesn't actually modify the data, it just runs the bonsai dedisperser.  
//
// FIXME: currently, there are two versions of the bonsai_dedisperser, written in python and C++.
// From python, they are constructed as 'bonsai_dedisperser' and 'bonsai_dedisperser_cpp' respectively.
// In the pipeline json output, they are represented as 'bonsai_dedisperser_python' and 'bonsai_dedisperser_cpp'.
// The two versions of the bonsai_dedisperser will be combined eventually!
//
// If the 'track_global_max' flag is set to true, then the following json output will be written:
//   frb_global_max_trigger
//   frb_global_max_trigger_dm
//   frb_global_max_trigger_tfinal
//
// We don't currently define any mechanism for the C++ bonsai_transform to write plots,
// but this should be easy to change if needed.  The python bonsai_transform does contain plotter logic.


struct bonsai_initializer {
    bool fill_rfi_mask = false;                  // If true, then online_mask_filler will be run (this makes a big difference!)
    bool use_analytic_normalization = false;     // If true, then unit-variance toy model is assumed (not suitable for real data!)
    bool track_global_max = false;               // If true, then global max trigger info will be written to pipeline json file
    int verbosity = 1;                           // Print some informational output (in constructor only)
    int dm_min = 0.0;                            // Only meaningful if track_global_max = True
    int dm_max = 0.0;                            // Only meaningful if track_global_max = True.  Zero means "no max DM".
    std::string file_type;                       // Filetype of config file.  Allowed values are { txt, hdf5 }.  Empty string means "infer from filename".
    std::string hdf5_output_filename;            // If this string is nonempty, then trigger HDF5 files will be written.
    int nt_per_hdf5_file = 0;                    // Only meaningful if hdf5_output_filename is nonempty.  Zero means "one big file".

    bonsai_initializer() { }
};

// This interface is similar to the python make_bonsai_dedisperser().
extern std::shared_ptr<wi_transform> make_bonsai_dedisperser(const std::string &config_filename, const bonsai_initializer &ini_params = bonsai_initializer());

// This interface may be more suitable for low-level use.
extern std::shared_ptr<wi_transform> make_bonsai_dedisperser(const std::shared_ptr<bonsai::dedisperser> &d);


// -------------------------------------------------------------------------------------------------
//
// gaussian_noise_stream: simple stream which simulates Gaussian random noise.
//
//   nfreq               Number of frequency channels
//   nt_tot              Total number of time samples written before stream ends.
//   freq_lo_MHz         Lowest frequency in band (e.g. 400 for CHIME)
//   freq_hi_MHz         Highest frequency in band (e.g. 800 for CHIME)
//   dt_sample           Length of a time sample in seconds
//   sample_rms          RMS of intensity samples (Gaussian distributed)
//   nt_chunk            Stream block size (if zero, will default to a reasonable value)
//   randomize_weights   If true, weights will be uniform random numbers (if false, all weights will be 1.0)


extern std::shared_ptr<wi_stream> make_gaussian_noise_stream(ssize_t nfreq, ssize_t nt_tot, double freq_lo_MHz, double freq_hi_MHz, 
							     double dt_sample, double sample_rms=1.0, ssize_t nt_chunk=0, bool randomize_weights=false);


// -------------------------------------------------------------------------------------------------
//
// CHIME streams.


// File streams, either from single file, explicit filename, or acquisition directory.
// In the 'acqusition directory' case, the directory is scanned for filenames of the form NNNNNNNN.h5, where N=[0,9].
// The 'nfiles' optional argument can be used to limit the acquisition to the first N files.
//    
// The 'nt_chunk' arg is the chunk size used internally when moving data from hdf5 file
// into the rf_pipelines buffer.  If unspecified or zero, it will default to a reasonable value.
//
// If 'noise_source_align' is nonzero, then it should be equal to the DETRENDER chunk size (not the chime_file_stream nt_chunk).
// In this case, the stream will align the noise source edges with the detrender chunks, by discarding initial data if necessary.
//
// Note: functions beginning 'make_chime_stream..' are HDF5 streams, whereas functions beginning 'make_chime_frb_stream...' are msgpack.\n"
// For example, make_chime_stream_from_filename() and make_chime_frb_stream_from_filename() create streams from a single HDF5 or msgpack\n"
// file, respectively.
//
// Note: a quick way to inspect a CHIME hdf5 file is using the 'ch-show-intensity-file' and 'ch-plot-intensity-file'
// programs, in the ch_frb_io github repo.

extern std::shared_ptr<wi_stream> make_chime_stream_from_filename(const std::string &filename, ssize_t nt_chunk=0, ssize_t noise_source_align=0);
extern std::shared_ptr<wi_stream> make_chime_stream_from_acqdir(const std::string &dirname, ssize_t nt_chunk=0, ssize_t noise_source_align=0, ssize_t nfiles=0);
extern std::shared_ptr<wi_stream> make_chime_stream_from_filename_list(const std::vector<std::string> &filename_list, ssize_t nt_chunk=0, ssize_t noise_source_align=0);

// CHIME assembled_chunk file stream, in msgpack format.
extern std::shared_ptr<wi_stream> make_chime_frb_stream_from_glob(const std::string &glob_pattern, ssize_t nt_chunk=0, ssize_t noise_source_align=0);
extern std::shared_ptr<wi_stream> make_chime_frb_stream_from_filename(const std::string &filename, ssize_t nt_chunk=0, ssize_t noise_source_align=0);
extern std::shared_ptr<wi_stream> make_chime_frb_stream_from_filename_list(const std::vector<std::string> &filename_list, ssize_t nt_chunk=0, ssize_t noise_source_align=0);

// CHIME network stream.  Receives UDP packets in "CHIME L0-L1 format".
// Assumes the ch_frb_io::intensity_network_stream object is already constructed (but not started).
// If the 'prescale' argument is specified, all intensity values will be multiplied by its value.
// This is a temporary workaround for 16-bit overflow issues in bonsai.
extern std::shared_ptr<wi_stream> make_chime_network_stream(const std::shared_ptr<ch_frb_io::intensity_network_stream> &stream, int beam_id, float prescale=1.0);

// A higher-level interface which constructs a ch_frb_io::intensity_network_stream expecting a single beam_id.  
// If the 'udp_port' argument is zero, then the default chimefrb port will be used.
extern std::shared_ptr<wi_stream> make_chime_network_stream(int udp_port=0, int beam_id=0, float prescale=1.0);

// "Dummy" CHIME network stream, intended for timing.
// Returns a stream which decodes a preallocated ch_frb_io::assembled_chunk pool as the pipeline progresses,
// but does not actually receive packets over the network.  This allows the CPU cost of assembled_chunk
// decoding to be included in pipeline timings.  Default parameter values are appropriate for full CHIME.
extern std::shared_ptr<wi_stream> make_dummy_chime_network_stream(ssize_t nt_tot, int nupfreq=16, int nt_per_packet=16, int fpga_counts_per_sample=384, double pool_gb=1.0);

// Experimental: masks "spikes" in 16K data.
extern std::shared_ptr<chunked_pipeline_object> make_chime_16k_spike_mask(ssize_t nt_chunk=0);

// Experimental: removes "ripples" from 16K data.
extern std::shared_ptr<chunked_pipeline_object> make_chime_16k_derippler(double fudge_factor=1.0, ssize_t nt_chunk=0);

// Experimental: analyzes 16k-ripples and writes result to HDF5 file for follow-up analysis.
extern std::shared_ptr<wi_transform> make_chime_16k_stripe_analyzer(ssize_t Dt1=16, ssize_t Df2=16, ssize_t Dt2=16);

// Experimental: writes HDF5 file containing the intensity spectrum.
extern std::shared_ptr<wi_transform> make_spectrum_analyzer(ssize_t Dt1=16, ssize_t Dt2=16);


// -------------------------------------------------------------------------------------------------
//
// CHIME output streams (these are "transforms" as far as rf_pipelines is concerned.)


// chime_file_writer.
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

std::shared_ptr<wi_transform> make_chime_file_writer(const std::string &filename, bool clobber=false, int bitshuffle=2, ssize_t nt_chunk=0);


// chime_packetizer.
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

extern std::shared_ptr<wi_transform> make_chime_packetizer(const std::string &dstname, int nfreq_coarse_per_packet, int nt_per_chunk, 
							   int nt_per_packet, float wt_cutoff, double target_gbps, int beam_id=0);


// -------------------------------------------------------------------------------------------------
//
// Mask counter transform -- counts masked data samples


struct mask_measurements {
    ssize_t pos = 0;   // index of first time sample (relative to start of pipeline, without any time downsampling factor applied)
    int nf = 0;        // number of frequencies (may be downsampled relative to "toplevel" frequency resolution in pipeline)
    int nt = 0;        // number of time samples (may be downsampled relative to "toplevel" frequency resolution in pipeline)
    int nsamples = 0;  // always equal to (nf * nt)
    int nsamples_unmasked = 0;  // number of unmasked samples (satisfies 0 <= nunmasked <= nsamples)

    // For each frequency, how many of the "nt" samples are masked?
    std::shared_ptr<int> freqs_unmasked;
};


class mask_measurements_ringbuf {
public:
    mask_measurements_ringbuf(int nhistory=300);

    std::unordered_map<std::string, float> get_stats(float period);
    std::vector<rf_pipelines::mask_measurements> get_all_measurements();

    void add(rf_pipelines::mask_measurements& meas);
    
protected:
    std::vector<rf_pipelines::mask_measurements> ringbuf;
    std::mutex mutex;
    int current;
    int maxsize;
};


class mask_counter_transform : public wi_transform {
public:
    std::string where;
    std::shared_ptr<mask_measurements_ringbuf> ringbuf;

    mask_counter_transform(int nt_chunk_, std::string where_, std::string class_name_="mask_counter");
    virtual ~mask_counter_transform() { }
    virtual void _process_chunk(float *intensity, ssize_t istride, float *weights, ssize_t wstride, ssize_t pos) override;
    virtual Json::Value jsonize() const override;
    static std::shared_ptr<mask_counter_transform> from_json(const Json::Value &j);

    virtual void process_measurement(mask_measurements& meas);
    void init_measurements(mask_measurements& meas);
    
    std::shared_ptr<mask_measurements_ringbuf> get_ringbuf();
};

class chime_mask_counter : public mask_counter_transform {
public:
    chime_mask_counter(std::string where_);
    virtual ~chime_mask_counter() { }

    virtual void _bind_transform(Json::Value &json_attrs) override;
    virtual void _start_pipeline(Json::Value &j) override;
    virtual void _process_chunk(float *intensity, ssize_t istride, float *weights, ssize_t wstride, ssize_t pos) override;
    
    virtual Json::Value jsonize() const override;
    static std::shared_ptr<chime_mask_counter> from_json(const Json::Value &j);

    // Called "externally" by L1 server
    void set_stream(std::shared_ptr<ch_frb_io::intensity_network_stream> stream, int beam);
                    
protected:
    std::shared_ptr<ch_frb_io::intensity_network_stream> stream;
    int beam = -1;

    // The FPGA count related fields are initialized in _start_pipeline().
    bool fpga_counts_initialized = false;
    uint64_t initial_fpga_count = 0;
    int fpga_counts_per_sample = 0;
};

// Externally callable
std::shared_ptr<wi_transform> make_mask_counter(int nt_chunk, std::string where);
std::shared_ptr<wi_transform> make_chime_mask_counter(std::string where);


// -------------------------------------------------------------------------------------------------
//
// chime_slow_pulsar_writer


struct chime_slow_pulsar_writer : public wi_transform
{
    struct real_time_state {
	// Throughout the CHIMEFRB pipeline, beams are identified by a beam_id between 0 and 1024.
	int beam_id = -1;

	// Use this to request memory from inside the CHIMEFRB L1 server (see ch_frb_io/ch_frb_io.hpp)
	std::shared_ptr<ch_frb_io::memory_slab_pool> memory_pool;
	
	// Use this to queue a write_request, for writing to disk by the L1 server I/O threads.
	std::shared_ptr<ch_frb_io::output_device_pool> output_devices;
    };

    struct output_file_params {
	int nfreq_out = 0;   // number of frequency channels in output file
	int nds_out = 0;     // time downsampling factor in output file
	int nbits_out = 0;   // bit depth in output file
    };
    
    real_time_state rt_state;
    output_file_params of_params;

    chime_slow_pulsar_writer(ssize_t nt_chunk);

    // Called by RPC thread, once during initialization.
    void init_real_time_state(const real_time_state &rt_state);

    // Called by RPC thread, intermittently while pipeline is running.
    void set_output_file_params(const output_file_params &of_params);

    // Called by rf_pipelines thread.
    virtual void _process_chunk(float *intensity, ssize_t istride, float *weights, ssize_t wstride, ssize_t pos) override;
    virtual void _end_pipeline(Json::Value &json_output) override;

    virtual Json::Value jsonize() const override;
    static std::shared_ptr<chime_slow_pulsar_writer> from_json(const Json::Value &j);
};


}  // namespace rf_pipelines

#endif // _RF_PIPELINES_INVENTORY_HPP
