#ifndef _RF_PIPELINES_BASE_CLASSES_HPP
#define _RF_PIPELINES_BASE_CLASSES_HPP

#if (__cplusplus < 201103) && !defined(__GXX_EXPERIMENTAL_CXX0X__)
#error "This source file needs to be compiled with C++11 support (g++ -std=c++11)"
#endif

#include <string>
#include <vector>
#include <memory>
#include <climits>
#include <iostream>
#include <unordered_map>
#include <json/json.h>


namespace rf_pipelines {
#if 0
}  // emacs pacifier
#endif


// -------------------------------------------------------------------------------------------------
//
// Here is a good place to explain the rf_pipelines state model.
//
// The first thing to say is that high-level users of rf_pipelines probably don't need to know about 
// the state model; they should just be able to call pipeline::run(), which can advance the state model 
// from its initial state (UNBOUND) to its final state (DONE).
//
// When a new pipeline_object is constructed, its state is UNBOUND, meaning that pipeline-wide parameters
// have not been determined yet (number of frequency channels, level of downsampling in each ring buffer, 
// etc.)  The pipeline_object only "knows" about parameters which were specified in its constructor.
//
// Next, the pipline is "bound", by calling pipeline_object::bind() for each pipeline_object in the pipeline.
// This advances the state from UNBOUND to BOUND, and does the following:
//
//   - Propagates pipeline attributes through the pipeline.  Any pipeline_object may define
//     attributes, and subsequent pipeline_objects will have access to them.  For example, a
//     stream object may define "dt_sample", the timestream sample length in seconds, and a
//     subsequent transform object can access the value.  The pipeline attributes are represented
//     as a JSON object.
//
//   - Determines which ring buffers exist in the pipeline, which pipeline_objects access them,
//     chunk sizes for each pipeline_object, and latency for each pipeline_object.  Thus, after
//     the pipeline is bound, its total memory usage and latency are known, although the memory
//     has not actually been allocated yet.
//
// Next, the pipeline is "allocated", which just means that its state advances from BOUND
// to ALLOCATED, and all ring buffers are allocated.  (Exception: in pipeline_objects which
// emit zoomable plots for the web viewer, some allocation happens later, either in start_pipeline(),
// or on-the-fly while the pipeline is running.  This is nontrivial to fix!  In realistic
// cases, I think the total memory footprint of the plotting buffers should be small.)
//
// Next, the pipeline is "started", which advances its state from ALLOCATED to RUNNING,
// and does the following:
//
//   - A second round of pipeline attributes are propagated through the pipeline (the first
//     round of attributes was propagated in bind()).  This is to allow for pipeline attributes
//     which are not known until the pipeline is actually started.  For example, in CHIME,
//     the initial timestamp (FPGA count) isn't known until the pipeline is started, since
//     we don't actually start listening for FRB network packets until then.
//
// Next, the pipeline runs, and its state advances from RUNNING to DONE.
// All plots and json output are written to disk here.
// 
// After the pipeline finishes, if you want to run it again (or re-use its constituent
// pipeline_objects in a different pipeline), you'll need to call one of the following:
//
//     pipeline_object::reset()        reverts state to ALLOCATED
//     pipeline_object::deallocate()   reverts state to BOUND
//     pipeline_object::unbind()       reverts state to UNBOUND
//
// Note that resetting the pipeline (without fully unbinding) is useful in a real-time pipeline,
// but probably doesn't make sense in an offline analysus.


// Defined in rf_pipelines_internals.hpp
struct outdir_manager;
struct plot_group;
struct zoomable_tileset_state;


// -------------------------------------------------------------------------------------------------
//
// run_params: parameters which are specified before running a pipeline.
//
// When a pipeline is run, these parameters get incorporated into its 'json attributes',
// which are available to all pipeline_objects in the _bind() virtual.  This is explained
// in more detail below.


struct run_params {
    
    // 'outdir' is the rf_pipelines output directory, where the rf_pipelines json file will
    // be written, in addition to other transform-specific output files such as plots. 
    //
    // If 'outdir' is an empty string, then the json file will not be written, and 
    // any transform which tries to write an output file (such as a plotter_transform)
    // will throw an exception.
    //
    // If 'clobber' is false, then an exception will be thrown if the pipeline tries to
    // overwrite an old rf_pipelines.json file.
    //
    // Plot-related params:
    //   img_nzoom = number of zoom levels (FIXME currently hardcoded, should switch to adaptive)
    //   img_nds = time downsampling factor of plots at lowest zoom level
    //   img_nx = number of x-pixels (i.e. time axis) in each plot tile
    //
    // The meaning of the 'verbosity' argument is:
    //   0 = no output
    //   1 = high-level summary output (names of transforms, number of samples processed etc.)
    //   2 = log all output files
    //   3 = debug trace through pipeline
    //
    //
    // FIXME: verbosity not actually implemented yet.
    
    std::string outdir = ".";
    bool clobber = true;
    ssize_t img_nzoom = 4;
    ssize_t img_nds = 16;
    ssize_t img_nx = 256;
    int verbosity = 2;
    bool debug = false;

    void check() const;
};


// -------------------------------------------------------------------------------------------------
//
// ring_buffer, ring_buffer_dict, ring_buffer_subarray
//
// An important note to keep in mind when using the ring_buffer!  The ring buffer can be
// downsampled relative to the "native" time resolution of the pipeline.  The amount of
// downsampling is given by the parameter 'nds'.
//
// In general, time indices are always sample counts at native (non-downsampled) resolution.
// For example, the 'nt_contig' and 'nt_maxlag' arguments to ring_buffer::update_params(),
// and the 'pos0' and 'pos1' arguments to ring_buffer::get().
//
// However, array dimensions always have the downsampling factor applied.  For example,
// the pointer returned by ring_buffer::get() points to an array whose last dimension
// is (pos1-pos0)/nds, not (pos1-pos0).


class ring_buffer {
public:
    // "Complementary" dimensions (all dimensions except time axis)
    const std::vector<ssize_t> cdims;
    const ssize_t csize;  // product of all cdims
    
    // Downsampling factor
    const ssize_t nds;                 

    ring_buffer(const std::vector<ssize_t> &cdims, ssize_t nds);
    ~ring_buffer();

    // The 'nt_contig' and 'nt_maxlag' arguments do not have the downsampling factor 'nds' applied.
    void update_params(ssize_t nt_contig, ssize_t nt_maxlag);

    void allocate();
    void deallocate();
    void reset();
    
    // The access_mode is optional, but enables some debug checks, and can
    // also help performance.  The numerical values are chosen for convenient
    // bitmasking.  The ACCESS_NONE value is a placeholder which throws an exception.

    static constexpr int ACCESS_NONE = 0;
    static constexpr int ACCESS_READ = 0x01;
    static constexpr int ACCESS_WRITE = 0x02;
    static constexpr int ACCESS_RW = 0x03;
    static constexpr int ACCESS_APPEND = 0x06;

    // Each call to get() must be followed by a call to put(), which "returns" the memory reference.
    // However, rather than using ring_buffer::get(), ring_buffer::put() directly, it's usually
    // preferable to the RAII-wrapper class ring_buffer_subarray below.
    //
    // The 'pos0' and 'pos1' arguments do not have the downsampling factor 'nds' applied.
    // However, the pointer returned by get() points to an array whose dimensions do have the
    // downsampling factor applied, i.e. the last dimension of the array is (pos1-pos0)/nds,
    // not (pos1-pos0).

    float *get(ssize_t pos0, ssize_t pos1, int access_mode);
    void put(float *p, ssize_t pos0, ssize_t pos1, int access_mode);

    ssize_t get_stride() const;

    // The ring_buffer is noncopyable, since it contains a bare pointer.
    // Move constructor/assignment would be trivial to implement, but I haven't needed it yet.
    // (I always use ring_buffers via a shared_ptr<ring_buffer>.)

    ring_buffer(ring_buffer &&) = delete;
    ring_buffer(const ring_buffer &) = delete;
    ring_buffer& operator=(ring_buffer &&) = delete;
    ring_buffer& operator=(const ring_buffer &) = delete;

    static std::string access_mode_to_string(int access_mode);
    static void check_cdims(const std::vector<ssize_t> &cdims);

protected:
    // These parameters are initialized by repeated calls to update_params(),
    // before the ring buffer is allocated.
    ssize_t nt_contig = 0;   // no downsampling factor applied
    ssize_t nt_maxlag = 0;   // no downsampling factor applied

    // These parameters are initialized when the buffer is allocated.
    ssize_t period = 0;     // downsampling factor applied
    ssize_t stride = 0;     // downsampling factor applied
    float *buf = nullptr;

    // Runtime state
    ssize_t curr_pos = 0;            // downsampling factor applied
    ssize_t first_valid_sample = 0;  // downsampling factor applied
    ssize_t last_valid_sample = 0;   // downsampling factor applied

    // Is there an active pointer?
    float *ap = nullptr;
    float ap_pos0 = 0;  // no downsampling factor applied
    float ap_pos1 = 0;  // no downsampling factor applied
    float ap_mode = ACCESS_NONE;

    // Helper functions for internal use.
    // All time indices have the downsampling factor applied.
    void _mirror_initial(ssize_t it0);
    void _mirror_final(ssize_t it1);
    void _copy(ssize_t it_dst, ssize_t it_src, ssize_t n);

    friend struct ring_buffer_subarray;
};


// RAII wrapper for ring_buffer::get() and ring_buffer::put().
struct ring_buffer_subarray {
    std::shared_ptr<ring_buffer> buf;
    float *data = nullptr;
    ssize_t stride = 0;

    ssize_t access_mode = ring_buffer::ACCESS_NONE;
    ssize_t pos0 = 0;
    ssize_t pos1 = 0;

    ring_buffer_subarray() { }

    ring_buffer_subarray(const std::shared_ptr<ring_buffer> &buf_, ssize_t pos0_, ssize_t pos1_, int access_mode_) :
	buf(buf_),
	data(buf_->get(pos0_,pos1_,access_mode_)),
	stride(buf_->stride),
	access_mode(access_mode_),
	pos0(pos0_),
	pos1(pos1_)
    { }

    ~ring_buffer_subarray()
    {
	// FIXME how to handle case where ring_buffer::put() throws an exception?
	if (buf)
	    buf->put(data, pos0, pos1, access_mode);
    }

    // Make noncopyable (This is a lightweight class designed to live briefly on the call stack.)
    // Move constructor/assignment would be trivial to implement, but I haven't needed it yet.
    ring_buffer_subarray(ring_buffer_subarray &&) = delete;
    ring_buffer_subarray(const ring_buffer_subarray &) = delete;
    ring_buffer_subarray &operator=(ring_buffer_subarray &&) = delete;
    ring_buffer_subarray &operator=(const ring_buffer_subarray &) = delete;

    inline void get(const std::shared_ptr<ring_buffer> &buf_, ssize_t pos0_, ssize_t pos1_, int access_mode_)
    {
	this->reset();
	this->data = buf_->get(pos0_, pos1_, access_mode_);
	this->stride = buf_->stride;
	this->buf = buf_;
	this->pos0 = pos0_;
	this->pos1 = pos1_;
	this->access_mode = access_mode_;
    }

    inline void reset()
    {
	if (buf) {
	    buf->put(data, pos0, pos1, access_mode);
	    buf.reset();
	}

	data = nullptr;
    }
};


using ring_buffer_dict = std::unordered_map<std::string, std::shared_ptr<ring_buffer>>;


// -------------------------------------------------------------------------------------------------
//
// zoomable_tileset
//
// Virtual base class which represents one set of tiled plots for the web_viewer.
//
// The high-level logic is as follows:
// 
//    - Each zoomable_tileset maintains a vector of one or more ring buffers (rbvec).
//      The length of the vector (i.e. number of ring buffers) is subclass-dependent.
//      For example, a bonsai plotter might maintain one ring buffer per dedispersion
//      tree, and an intensity plotter might maintain two ring buffers (intensity, weights).
//
//    - At each zoom level, there is an additional rbvec which is downsampled by a
//      factor two, relative to the previous rbvec.  The zoomable_tileset defines a
//      virtual downsample_rbvec() which fills each rbvec_t by downsampling the previous one.
//
//    - The bottom-level (i.e. lowest zoom level) rbvec is filled incrementally as
//      the pipeline runs.  Rather than defining a new virtual for this, we just do
//      it in the usual pipeline-processing virtual (i.e. pipeline_object::_advance()
//      or chunked_pipeline_object::_process_chunk()).
//
//    - At each zoom level, when enough data has accumulated to emit a plot, the
//      zoomable_tileset defines a virtual which fills an RGB array from the rbvec.


struct zoomable_tileset {
    
    // At each zoom level, the zoomable_tileset maintains a vector of ring buffers (rbvec_t).  
    // All ring buffers in the same rbvec_t have the same downsmapling factor 'nds', but the 
    // value of nds differs by a factor two in consectutive zoom levels.

    using rbvec_t = std::vector<std::shared_ptr<ring_buffer>>;

    // Recall that the 'cdims' of a ring buffer are a vector<ssize_t> corresponding to all
    // array dimensions other than time.  A zoomable_tileset has one set of cdims for each
    // of its ring buffers, so zoomable_tileset::cdims is a vector<vector<ssize_t>>.
    //
    // Note that the number of ring buffers in the zoomable_tileset is cdims.size(), and
    // the array dimension of ring buffer 'i' is (cdims[i].size()+1).

    const std::vector<std::vector<ssize_t>> cdims;
    
    // Here is the complete list of plotting-related parameters:
    //
    //   img_prefix: determines tile filenames as ${outdir}/${img_prefix}_${izoom}_${ifile}.png
    //   img_nds: time downsampling of tiles at lowest zoom level (i.e. number of time samples per x-pixel)
    //   img_nzoom: number of zoom levels plotted
    //   img_nx: number of x-pixels per tile
    //   img_ny: number of y-pixels per tile
    //   nds_arr: time downsampling of ring buffer and RGB arrays at lowest zoom level
    //   ny_arr: number of y-pixels in RGB arrays
    // 
    // Note 1: the parameters (img_nzoom, img_nds, img_nx) are the same for every tileset
    //   emitted by the pipeline.  These parameters are in 'struct run_params' (which is
    //   pipeline-global), and are not available in 'struct zoomable_tileset'.
    //
    // Note 2: the global parameter 'img_nds' need not be the same as 'nds_arr', although
    //   there is a constraint that both must be powers of two.  From the persepective of
    //   the zoomable_tileset, nds_arr is the parameter that "matters".  (More precisely,
    //   when the lowest-level ring buffers are filled during pipeline processing, the
    //   downsampling level of the buffers is nds_arr.)
    //
    //   In the case where img_nds and nds_arr are not equal, rf_pipelines will appropriately 
    //   downsample or upsample before emitting plots.  This gives each pipeline_object the
    //   flexibility to choose its internal downsampling level depending on what is most
    //   convenient.
    //
    // Note 3: similarly, the parameters 'img_ny' and 'ny_arr' need not be equal, and
    //   from the perspective of the zoomable_tileset, ny_arr is the one that "matters".
    //   (More precisely, when the plot_rbvec() virtual is called, its output RGB array
    //   has y-dimension ny_arr, not img_ny.)
    //
    //   Here there is a constraint that ny_arr is equal to, or a divisor of, img_ny.
    //   This gives the pipeline_object the flexibility to choose ny_arr < img_ny, and
    //   rf_pipelines will automatically upsample the images.  However, in the case where 
    //   the "natural" y-dimension is larger than img_ny, there is no way to get automatic
    //   downsampling - the pipeline_object must supply its own downsampling logic in this case.

    const std::string img_prefix;
    const ssize_t img_ny;
    const ssize_t ny_arr;
    const ssize_t nds_arr;

    // Base class constructor.
    zoomable_tileset(const std::vector<std::vector<ssize_t>> &cdims, const std::string &img_prefix,
		     ssize_t img_ny, ssize_t ny_arr, ssize_t nds_arr);
	
    // Virtuals.
    //
    // Following a general convention in rf_pipelines, the 'pos' and 'nt' arguments
    // to these routines do not have any downsampling factors applied!
    //
    // plot_rbvec(rgb_out, rb_in, pos, nt)
    //   Creates an RGB image from ring buffer contents.
    //   Can be called at any zoom level, so the downsampling factor of 'rb_in' can be arbitrary.
    //   The output 'rgb_out' has shape (ny_arr, nt/nds, 3), where nds is the downsampling factor.
    //
    // downsample_rbvec(rb_out, rb_in, pos, nt)
    //   Note that 'rb_out' has twice the downsampling factor of 'rb_in'.
    //
    // extend_rbvec(rb_out, pos, nt)
    //   By default, this zeroes the buffers, but this can be overridden.

    virtual void plot_rbvec(uint8_t *rgb_out, rbvec_t &rb_in, ssize_t pos, ssize_t nt) = 0;
    virtual void downsample_rbvec(rbvec_t &rb_out, rbvec_t &rb_in, ssize_t pos, ssize_t nt);
    virtual void extend_rbvec(rbvec_t &rb_out, ssize_t pos, ssize_t nt);
};


// -------------------------------------------------------------------------------------------------
//
// Some abstract base classes: pipeline_object, chunked_pipeline_object, wi_stream, wi_transform.
//
// Note that all members/methods are public, due to complications in python-wrapping protected methods
// that I hope to fix eventually.


struct pipeline_object {
public:
    // Constructor for this abstract base class.
    // The name must be initialized at construction (possibly to something simple like the class name),
    // but can be changed later to something more verbose.
    explicit pipeline_object(const std::string &name);

    // High-level API: to run a pipeline, just call run().	
    //
    // If 'callback' is specified, it will be called periodically as callback(pos_lo, pos_hi),
    // where pos_lo is the number of samples which have been finalized by the pipeline, and
    // pos_hi is the number of partially processed samples.  (Note pos_lo <= pos_hi.)

    using callback_t = std::function<void(ssize_t,ssize_t)>;
    Json::Value run(const run_params &params = run_params(), const callback_t &callback = callback_t());

    // A more fine-grained high-level API.
    // See comment at top of this source file for more explanation!

    void bind(const run_params &params);  // does global pipeline initializations (advances state UNBOUND -> BOUND)
    void allocate();                      // allocates all buffers (advances state BOUND -> ALLOCATED)
    void reset();                         // reverts state to ALLOCATED
    void deallocate();                    // reverts state to BOUND
    void unbind();                        // reverts state to UNBOUND

    // Everything which follows is low-level stuff, which should not be needed by high-level users
    // of rf_pipelines, but may be needed if you're writing your own pipeline_object.

    enum {
	UNBOUND = 0,
	BOUND = 1,
	ALLOCATED = 2,
	RUNNING = 3,
	DONE = 4
    } state;

    std::string name;

    // General note: all time indices (for example pipeline_object::nt_chunk_*, or pipeline_object::pos_*)
    // use the "native" time resolution of the pipeline, i.e. they are not divided by any downsampling factor 'nds'.
    //
    // Runtime state:
    //   pos_lo = number of time samples which have been processed by this pipeline_object
    //   pos_hi = number of time samples which are ready for processing by this pipeline_object
    //   pos_max = number of time samples which have entered the pipeline
    //
    // Note that pos_lo <= pos_hi <= pos_max.

    ssize_t pos_lo = 0;   // always a multiple of nt_chunk_out
    ssize_t pos_hi = 0;   // always a multiple of nt_chunk_in
    ssize_t pos_max = 0;

    // These parameters control the flow of data into the pipeline_object.
    // They are set "externally", just before the virtual function _bind() is called.
    // Don't initialize or change them in the pipeline_object subclass, or strange things will happen!

    ssize_t nt_chunk_in = 0;   // "input" chunk size (more precisely, step size in pos_hi between calls to _advance())
    ssize_t nt_maxlag = 0;     // maximum lag when _advance() is called (more precisely, max difference between pos_hi and pos_max)

    // These parameters control the flow of data out of the pipeline object.
    // The subclass must initialize them, in the virtual function _bind().
    // Don't change them after _bind() returns, or strange things will happen!
    
    ssize_t nt_chunk_out = 0;  // "output" chunk size (more precisely, step size in pos_lo)
    ssize_t nt_maxgap = -1;    // maximum gap after _advance() returns (more precisely, max allowed value of (pos_hi-pos_lo))
    ssize_t nt_contig = 0;     // max contiguous chunk size this pipeline_object will request from ring buffers

    // Used internally, to keep track of which ring_buffers are bound to this pipeline_object.
    std::vector<std::shared_ptr<ring_buffer>> new_ring_buffers;    // ring buffers created by this pipeline object
    std::vector<std::shared_ptr<ring_buffer>> all_ring_buffers;    // all ring buffers used by this pipeline object (including new_ring_buffers)
    std::vector<std::shared_ptr<zoomable_tileset_state>> zoomable_tilesets;

    // Initialized in bind().
    // Note: if the pipeline is run without an output directory, then 'out_mp' will be 
    // a nonempty pointer, but out_mp->outdir will be an empty string.
    std::shared_ptr<outdir_manager> out_mp;
    double time_spent_in_transform = 0.0;

    // New plot_groups should be created in _start_pipeline(), by calling pipeline_object::add_plot_group().
    std::vector<plot_group> plot_groups;

    // Used internally, to save pipeline attributes from bind() and start_pipeline() respectively.
    // Note: only in the "top-level" pipeline!  (where run() is called)
    Json::Value json_attrs1;
    Json::Value json_attrs2;

    // Here are the virtuals which can/must be implemented, in order to define a pipeline_object subclass.
    // The meaning of each of the virtuals needs some explanation, see the long comment below for details!

    // Mandatory (pure virtual).
    virtual void _bind(ring_buffer_dict &rb_dict, Json::Value &json_attrs) = 0;    
    virtual ssize_t _advance() = 0;

    // Optional, but defaults throw exceptions.
    virtual ssize_t get_preferred_chunk_size();  // must be defined for stream-type objects which are first in pipeline
    virtual Json::Value jsonize() const;         // must be defined in order to serialize to json

    // Optional, and defaults do nothing.
    virtual void _allocate();
    virtual void _deallocate();
    virtual void _start_pipeline(Json::Value &json_attrs);
    virtual void _end_pipeline(Json::Value &json_output);
    virtual void _reset();
    virtual void _unbind();

    // Each of the following methods is a wrapper around the corresponding virtual function.
    // For example, bind() contains "generic" logic, and wraps _bind() which contains 
    // subclass-dependent logic.
    //
    // Note that "container" objects (for example, class pipeline) need to 
    // recurse into their children.  This should be done using the wrappers.
    // For example, pipeline::_bind() calls bind() in each child (not _bind()).

    void bind(ring_buffer_dict &rb_dict, ssize_t nt_chunk_in, ssize_t nt_maxlag, Json::Value &json_attrs, const std::shared_ptr<outdir_manager> &out_mp);
    void start_pipeline(Json::Value &json_attrs);
    void end_pipeline(Json::Value &json_output);
    ssize_t advance(ssize_t pos_hi, ssize_t pos_max);    


    // Here is a long comment explaning each of the virtuals above!
    //
    //
    // _bind(rb_dict, json_attrs): does global pipeline initializations, such as determining ring buffer sizes and latency.
    //
    //     Before _bind() is called, the following members are initialized
    //        ssize_t nt_chunk_in = 0;      "input" chunk size (more precisely, step size in pos_hi between calls to _advance())
    //        ssize_t nt_maxlag = 0;        maximum lag when _advance() is called (more precisely, max difference between pos_hi and pos_max)
    //
    //     The _bind() virtual is responsible for initializing the following members:
    //        ssize_t nt_chunk_out = 0;     "output" chunk size (more precisely, step size in pos_lo)
    //        ssize_t nt_maxgap = -1;       maximum gap after _advance() returns (more precisely, max allowed value of (pos_hi-pos_lo))
    //        ssize_t nt_contig = 0;        max contiguous chunk size this pipeline_object will request from ring buffers
    //
    //     It is also responsible for calling get_buffer() or create_buffer(), for each ring buffer which will be
    //     accessed by the pipeline_object.  (See below for more info on these member functions.)
    //
    //     Similarly, bind must call add_zoomable_tileset() for each zoomable_tileset which will be emitted
    //     by the pipeline_object.  (Note that there are currently two redundant API's in rf_pipelines for
    //     managing plots, the "zooomable_tileset" API which is called in bind(), and the "plot_group"
    //     API which is called in start_pipeline().  I'm trying to phase out the plot_group API in favor
    //     of the zoomable_tileset API.)
    //
    //     Finally, _bind() may optionally add new pipeline attributes to 'json_attrs', or read existing attributes.
    //     Note that each field of 'struct run_params' is incorported into the json_attrs, e.g. json_attrs['verbosity']
    //     is defined.
    //
    //
    // _allocate(): does subclass-specific buffer allocation.
    //
    //     Note that pipeline ring buffers are allocated separately, in the "generic" part of the pipeline,
    //     and do not need to be allocated in _allocate().  Most pipeline_object subclasses will not need to 
    //     define _allocate().
    //
    //
    // _start_pipeline(json_attrs): called just before pipeline starts running.
    //
    //     For pipeline_objects which use the deprecated plot_group API, the _start_pipeline() virtual
    //     is responsible for calling add_plot_group().  (Note that there are currently two redundant 
    //     API's in rf_pipelines for managing plots, the "zooomable_tileset" API which is called in bind(), 
    //     and the "plot_group" API which is called in start_pipeline().  I'm trying to phase out the 
    //     plot_group API in favor of the zoomable_tileset API.)
    //
    //     Optionally, _start_pipeline() can add new pipeline attributes to 'json_attrs', or read values of existing attributes.
    //     Note that there are two "rounds" of pipeline attribute propagation, one in _bind() and one in _start_pipeline().
    //     The purpose of the second round is to allow for pipeline attributes which aren't known until the pipeline is started.
    //
    //
    // _advance(): this pure virtual is the "core" routine which processes data.
    //
    //     Before _advance() is called, the following members are initalized:
    //        ssize_t pos_lo;     number of time samples which have been processed by this pipeline_object
    //        ssize_t pos_hi;     number of time samples which are ready for processing by this pipeline_object
    //        ssize_t pos_max;    number of time samples which have entered the pipeline
    //
    //     Note that pos_lo <= pos_hi <= pos_max.
    //
    //     The _advance() virtual is then responsible for processing data, and advancing pos_lo
    //     appropriately.  (Subject to the constraints that pos_lo is advanced by a multiple of
    //     nt_chunk_out, and that the final value of pos_lo is >= pos_hi - nt_maxlag.)
    //
    //     The return value from _advance() is a time sample count 'nt_end'.  The pipeline ends when the number
    //     of processed samples exceeds the min(nt_end), where the minimum is taken over all _advance() calls
    //     to all pipeline_objects in the pipeline.  For example:
    //
    //       - a "transform-type" object which does not define its own end-of-stream
    //         always returns SSIZE_MAX
    //
    //       - a "stream-type" object will return SSIZE_MAX until end-of-stream is 
    //         encountered, then returns the number of samples in the stream
    //
    //       - if a pipeline_object wants to "hard-terminate" the pipeline (but does not
    //         want to throw an exception), it can return zero.
    //
    //
    // _end_pipeline(json_output)
    //
    //     The _end_pipeline() virtual is responsible for adding fields to the pipeline_object's json output.
    //     By default, the pipeline defines a few generic fields ('name', 'cpu_time', etc.), so defining 
    //     _end_pipeline() is optional.
    //
    //
    // _reset(), _deallocate(), _unbind()
    //
    //     These are all similar and are called by the wrappers reset(), deallocate(), unbind().
    //
    //     Note that deallocation of the pipeline ring buffers is done separately, and this does
    //     not need to be done in _deallocate().
    //
    // get_preferred_chunk_size(): defines chunk size for stream-type object
    //
    //    This is only called on the first pipeline_object in the pipeline (subsequent chunk sizes are determined
    //    automatically).  The default virtual returns 0, which results in an exception ""...: this pipeline_object 
    //    cannot appear first in pipeline".  Stream-type pipeline_objects which can appear first in a pipeline 
    //    should override get_preferred_chunk_size() to return a nonzero value.
    //
    //
    // jsonize(): defines json seralization for pipeline_object.
    //
    //     A pipeline_object subclass which implements jsonize() will also want to implement 
    //     and "register" a deserializer.  See 
    //     from_json() static member function, and "register" it with register_json_deserializer().


    // Helper functions called by _bind(), to manage ring buffers which are used/created by the pipeline_object.
    //   get_buffer(): called for each pre-existing pipeline ring buffer which the pipeline_object uses.
    //   create_buffer(): called for each new pipeline ring buffer which the pipeline_object creates.
    //   add_zoomable_tileset(): called to create a zoomable_tileset for the web viewer.

    std::shared_ptr<ring_buffer> get_buffer(ring_buffer_dict &rb_dict, const std::string &bufname);
    std::shared_ptr<ring_buffer> create_buffer(ring_buffer_dict &rb_dict, const std::string &bufname, const std::vector<ssize_t> &cdims, ssize_t nds);

    // Called from bind().  
    // The 'zt' arg should be a new zoomable_tileset object constructed by the caller.
    // The 'json_attrs' arg should be the same as the 'json_attrs' argument to bind().
    // The return value is a vector containing the ring buffers at the lowest zoom level (see comments in 'struct zoomable_tileset').
    std::vector<std::shared_ptr<ring_buffer>> add_zoomable_tileset(const std::shared_ptr<zoomable_tileset> &zt, const Json::Value &json_attrs);

    
    // This is the deprecated "plot_group" API for managing plots, which is called from start_pipeline().
    //
    // FIXME: phase this out in favor of the "zoomable_tileset" API!
    //
    // add_plot_group(): 
    //
    //   Each pipeline_objects's output plots are divided into one or more "plot groups".
    //   For example, the bonsai_dedisperser can write one plot group per internally defined tree.
    //
    //   The 'nt_per_pix' arg is the number of pipeline time samples per x-pixel in the plot.
    //   The 'ny' arg is the number of y-pixels (assumed to be the same for all plots in the group).
    //   The return value is the group_id arg needed in add_plot(), and group_ids always go 0,1,...
    //
    //   Note: add_plot_group() should be called in _start_pipeline(), not in the pipeline_object
    //   constructor or in _bind().  This is because plot_groups are "reset" between pipeline runs.
    //
    // add_plot(): Call just before writing a plot.
    //
    //   The range of time samples in the plot is [it0:it0+nt), where (it0,nt) are "upsampled" 
    //   time sample indices at the "native" pipeline time resolution.  The pixel dimensions of 
    //   the plot are (nx,ny).
    //
    // 	 Note that nx is redundant, since nt should always equal (nt_per_pix * nds * nx), where nt_per_pix
    //   was specified in add_plot_group() and nds is the current time downsampling factor in the pipeline.
    //   Similarly, ny is redundant since it must match the value specified in add_plot_group(). Currently,
    //   we require the caller to specify these values so that they can be used for error checking.
    //
    //   The return value is the full pathname ('basename' with the pipeline output_dir prepended)
    //
    // add_file(): 
    //
    //   Call just before writing a non-plot file, to check for filename collisions between transforms.
    //   The return value is the full pathname ('basename' with stream output_dir prepended)

    int add_plot_group(const std::string &name, int nt_per_pix, int ny);   // returns group id
    std::string add_plot(const std::string &basename, int64_t it0, int nt, int nx, int ny, int group_id=0);
    std::string add_file(const std::string &basename);


    // Json deserialization
    //
    // pipeline_object::from_json() is the "master" deserializer, which accepts a json-serialized
    // pipeline object, and returns a shared_ptr<pipeline_object>.
    //
    // In order for a new pipeline_object subclass to support json-serialization, it must
    // define the virtual function jsonize().  To support deserialization, it must define
    // a "deserializer" (a function f(x) whose argument x is a Json::Value&, and returns
    // a shared_ptr<pipeline_object>), and "register" the deserializer by calling
    // pipeline_object::register_deserializer().
    //
    // By convention, the deserializer for a new class X is usually a static member function
    //    shared_ptr<pipeline_objecT> X::from_json(Json::Value &j)
    //
    // and the boilerplate for registering the deserializer is:
    //    pipeline_object::register_json_deserializer("X", X::from_json);
    //
    // (This boilerplate can be wrapped in a class constructor in an anonymous namespace, see
    // pipeline.cpp for an example.)

    static std::shared_ptr<pipeline_object> from_json(const Json::Value &j);

    using json_deserializer_t = std::function<std::shared_ptr<pipeline_object> (const Json::Value &)>;
    static void register_json_deserializer(const std::string &class_name, const json_deserializer_t &f);

    // Disable copy/move constructors (I always use pipeline_objects via shared_ptr<>).
    pipeline_object(pipeline_object &&) = delete;
    pipeline_object(const pipeline_object &) = delete;
    pipeline_object &operator=(pipeline_object &&) = delete;
    pipeline_object &operator=(const pipeline_object &) = delete;

    virtual ~pipeline_object();

    // Helper: throws runtime_error with prefix "rf_pipelines: <name>: ..."
    void _throw(const std::string &msg) const;

    // For debugging or internal use.
    static void _show_registered_json_deserializers();
    static json_deserializer_t _find_json_deserializer(const std::string &class_name);  // can return NULL.
};


// -------------------------------------------------------------------------------------------------
//
// chunked_pipeline_object: corresponds to a pipeline_object which processes data in fixed-size chunks.
//
// This is a "semi-abstract base class": it defines some virtuals in its pipeline_object base class,
// but leaves others to be defined by an additional level of subclassing.
 

struct chunked_pipeline_object : public pipeline_object {
public:

    // Note: inherits 'name' member from pipeline_object base class.
    //
    // The name must be initialized at construction (possibly to something simple like the class name),
    // but can be changed later to something more verbose.
    //
    // The 'can_be_first' parameter should be true for stream-type objects which can be first in
    // a pipeline, and false for transform-type objects which process existing data.

    chunked_pipeline_object(const std::string &name, bool can_be_first);

    const bool can_be_first;

    // The 'nt_chunk' parameter is the chunk size, in time samples, with no downsampling factor applied.
    // It can either be initialized to something nonzero, or determined automatically by the pipeline.
    //
    // (Reminder: each ring buffer has its own downsampling factor 'nds', and the number of buffer samples
    // in each chunk will be nt_chunk/nds.)
    //
    // The value of nt_chunk must be a multiple of the ring buffer downsampling factor, for each
    // ring buffer which is used ("bound") by the chunked_pipeline_object.  This is checked in
    // bind(), and an exception will be thrown on failure.
    //
    // The value of nt_chunk is initially zero, but it will be initialized to a nonzero value
    // in one of three ways:
    //
    //   - The subclass may initialize nt_chunk by hand, either in its constructor or in the
    //     subclass-defined virtual function _bindc().  In this case, the subclass is responsible
    //     for ensuring that nt_chunk is valid (e.g. a multiple of all relevant 'nds' values)
    //
    //   - If the subclass does not initialize nt_chunk to something nonzero, then the pipeline
    //     will assign a default value, just after _bindc() is called, by calling the helper method
    //     finalize_nt_chunk().  This has the potential disadvantage that in _bindc(), nt_chunk will 
    //     not be initialized yet.
    //
    //   - The subclass may call finalize_nt_chunk() in _bindc(), after ring buffers are bound, but
    //     before the value of nt_chunk is used.  This option makes sense for chunked_pipeline_objects
    //     which want nt_chunk to be determined automatically, but also need to know its value in _bindc().
    //
    // Note 1: the value of 'nt_chunk' should not be modifed after _bindc() returns, or strange
    // things will happen!
    //
    // Note 2: if 'can_be_first' is true, then nt_chunk must be initialized in the constructor.  In
    // particular, this means that it can't be determined automatically, as in options #2 and #3 above.

    ssize_t nt_chunk = 0;

    // Helper method which automatically chooses nt_chunk, if has not already been initialized to something nonzero.
    // (This is virtual because wi_transform defines 'kernel_chunk_size', which needs to be incorporated.)
    virtual void finalize_nt_chunk();

    // These virtuals in the pipeline_object base class are defined by 'chunked_pipeline_object'.
    // We make them 'final', so that if e.g. a subclass erroneously overrides _bind() instead of _bindc(),
    // this error will be detected by the compiler.

    virtual void _bind(ring_buffer_dict &rb_dict, Json::Value &json_attrs) final override;
    virtual void _unbind() final override;
    virtual ssize_t _advance() final override;
    virtual ssize_t get_preferred_chunk_size() final override;

    // New pure virtuals, to be defined by subclass.
    //
    // _bindc(): responsible for calling get_buffer() or create_buffer(), for each ring bufffer which will
    //   be accessed by the chunked_pipeline_object.  Optionally, _bindc() may also read/create json attributes.
    //
    // _process_chunk(pos): responsible for processing data in the range [pos, pos+nt_chunk).
    //   Note that this range of sample indices should be used when accessing pipeline ring buffers.
    //   The return value is a boolean which is usually true, but stream-type classes which have reached 
    //   end-of-stream should return false.
    //
    // _unbindc(): any special logic needed to undo _bindc() should go here.

    virtual void _bindc(ring_buffer_dict &rb_dict, Json::Value &json_attrs) = 0;
    virtual bool _process_chunk(ssize_t pos) = 0;
    virtual void _unbindc();

    // prebind_nt_chunk: saved value of nt_chunk, before it is finalized in bind().
    // This is useful in jsonize() and _unbind().

    ssize_t _prebind_nt_chunk = 0;   // intended to be accessed through get_prebind_nt_chunk()
    ssize_t get_prebind_nt_chunk() const { return (state >= BOUND) ? _prebind_nt_chunk : nt_chunk; }

    // Internal helper function, assumes nt_chunk has been initialized.
    // (This is virtual because wi_transform defines 'kernel_chunk_size', which needs to be incorporated.)
    virtual void _check_nt_chunk() const;

    // Subclass can optionally override: jsonize(), _allocate(), _deallocate(), _start_pipeline(), _end_pipeline(), _reset().
};


// -------------------------------------------------------------------------------------------------
//
// wi_stream: "weighted intensity stream"
//
// Represents a stream-type object which generates weights and intensity arrays.  These are
// generated in regular chunks (i.e. wi_stream is a subclass of chunked_pipeline_object).
//
// This is a "semi-abstract base class": it defines some virtuals in its pipeline_object base class,
// but leaves others to be defined by an additional level of subclassing.


struct wi_stream : chunked_pipeline_object {
public:

    // Note: inherits 'name' member from base class.  The name must be initialized at construction 
    // (say to something simple, like the class name), but can be changed later to something more verbose.

    wi_stream(const std::string &name);

    // Note: inherits 'nt_chunk' member from base class.
    //
    // The values of 'nfreq' and 'nt_chunk' must be initialized to nonzero values, either
    // in the constructor, or in the _bind_stream() virtual.  (Note that nfreq and nt_chunk
    // can be initialized anywhere, in contrast with wi_transform where the rules are
    // more complicated, see below.)

    ssize_t nfreq = 0;

    // These virtuals in the chunked_pipeline_object base class are defined by 'wi_stream'.
    virtual void _bindc(ring_buffer_dict &rb_dict, Json::Value &json_attrs) final override;
    virtual bool _process_chunk(ssize_t pos) final override;
    virtual void _unbindc() final override;

    // New virtuals, to be defined by subclass.
    //
    // _bind_stream(json_attrs)
    //
    //    This is the "last chance" to initialize 'nfreq', 'nt_chunk', if these members
    //    were not already initialized in the constructor.
    //
    //    Optionally, _bind_stream() may read/create json attributes.  In CHIME, the stream
    //    is responsible for creating attributes 'freq_lo_MHz', 'freq_hi_MHz', 'dt_sample'.
    //
    // _unbind_stream()
    //
    //    Any code needed to undo _bind_stream() can go here.  (Usually not needed.)
    //
    // _fill_chunk(intensity, istride, weights, wstride, pos)
    //
    //    This is the "core" method which is responsible for filling the 'intensity' and 'weights'
    //    arrays.  These are arrays of shape (nfreq, nt_chunk), which must be filled with data
    //    corresponding to sample range [pos, pos+nt_chunk).  The memory strides are istride/wstride,
    //    i.e. the (i,j)-th element of the intensity array is intensity[i*istride+j], and similarly
    //    for the weights.
    //
    //    The return value should be 'true' normally, or 'false' if end-of-stream has been reached.

    virtual void _bind_stream(Json::Value &json_attrs);  // non-pure virtual (default does nothing)
    virtual bool _fill_chunk(float *intensity, ssize_t istride, float *weights, ssize_t wstride, ssize_t pos) = 0;
    virtual void _unbind_stream();  // non-pure virtual (default does nothing)

    // The wi_stream subclass shouldn't need to use these directly, 
    // since its _fill_chunk() method operates directly on pointers/strides.
    std::shared_ptr<ring_buffer> rb_intensity;
    std::shared_ptr<ring_buffer> rb_weights;

    // Subclass can optionally override: jsonize(), _allocate(), _deallocate(), _start_pipeline(), _end_pipeline(), _reset().
};


// -------------------------------------------------------------------------------------------------
//
// wi_transform: "weighted intensity transform"
// 
// Represents a transform-type object which reads and/or modifies weights and intensity arrays,
// usually as a processing step in a larger pipeline.  The processing is done in regular chunks
// (i.e. wi_transform is a subclass of chunked_pipeline_object).
//
// This is a "semi-abstract base class": it defines some virtuals in its pipeline_object base class,
// but leaves others to be defined by an additional level of subclassing.
//
// An important note!  The intensity and weights arrays which are processed by the wi_transform
// may be downsampled (in time) relative to the "native" pipeline resolution.  The level of 
// downsampling is determined by where the transform is placed into a larger pipeline, not by 
// the transform itself.  Currently, the only way downsampling can arise is if the transform 
// is placed in a downsampled "sub-pipeline" (see class wi_sub_pipeline).
//
// The 'nds' member is the downsampling factor, relative to the native resolution (i.e. nds=1 means
// there is no downsampling, and nds > 1 is the downsampled case).  The value of nds is initialized 
// just before the subclass-defined virtual _bind_transform() is called.
//
// In the downsampled case (i.e. nds > 1), the 'intensity' and 'weights' arrays which are processed 
// by the transform have length (nt_chunk/nds), not length nt_chunk.  This is consistent with a 
// general rf_pipelines convention that time indices (e.g. 'nt_chunk', or the 'pos' argument to 
// _process_chunk()) do not have downsampling factors applied, but array dimensions do have 
// downsampling factors applied.


struct wi_transform : chunked_pipeline_object {
public:

    // Note: inherits 'name' member from base class.  The name must be initialized at construction 
    // (say to something simple, like the class name), but can be changed later to something more verbose.

    wi_transform(const std::string &name);

    // The rules for initializing 'nfreq' and 'nds' are as follows:
    //
    //   - If nfreq and nds are changed from their default values, this should
    //     be done in the subclass constructor.  Changing their values in
    //     _bind_transform() is "too late"!
    //
    //   - If nfreq and/or nds are nonzero, this means "a specific value is required".
    //     For example, if nfreq is set to 1024, then an exception will be thrown in
    //     bind() if the number of frequency channels in the data is not equal to 1024.
    //
    //   - If nfreq and/or nds are zero, this means "any value is allowed, and initialize
    //     before calling _bind_transform()".  For example, if nds is zero, then the
    //     value of nds will be initialized to the actual downsampling factor of the
    //     data, just before _bind_transform() is called.
    //
    //   - Note that the default value of nfreq is 0 (i.e. any number of frequency channels
    //     is allowed), but the default value of nds is 1 (i.e. exception is thrown if the
    //     transform is downsampled).  Transforms which support downsampling should set
    //     nds to 0 in their subclass constructor, and take special care to ensure correctness
    //     in the downsampled case.  In particular, remember that time indices are at "native" 
    //     resolution but array dimensions are downsampled, e.g. the array arguments to
    //     _process_chunk() have length (nt_chunk/nds), not length nt_chunk.

    ssize_t nfreq = 0;  // Number of frequency channels.
    ssize_t nds = 1;    // Time downsampling factor, relative to "native" time resolution of pipeline.

    // Note: inherits 'nt_chunk' and finalize_nt_chunk() from base class.
    //
    // We also define 'kernel_chunk_size'.  If nonzero, then there is a requirement that nt_chunk
    // be a multiple of (kernel_chunk_size * nds).  A typical use case is a kernel which requires
    // that the array length be a multiple of 8 (say) for vectorization.
    //
    // The rules for initializing 'nt_chunk', are the same as in the chunked_pipeline_object base class, 
    // except that wi_transform::_bind_transform() plays the role of chunked_pipeline_object::_bindc().  
    //
    // For more details, please refer to the chunked_pipeline_object documentation.  The upshot is that 
    // if you need to use the value of nt_chunk in _bind_transform(), you need to either initialize it 
    // to something nonzero, or call finalize_nt_chunk() to determine it automatically.

    ssize_t kernel_chunk_size = 0;

    // These chunked_pipeline_object virtuals are defined here.
    virtual void _bindc(ring_buffer_dict &rb_dict, Json::Value &json_attrs) final override;
    virtual bool _process_chunk(ssize_t pos) final override;
    virtual void _unbindc() final override;

    // New virtuals, to be defined by subclass.
    //
    // _bind_transform(json_attrs)
    //
    //   The values of nfreq and nds are initialized just before _bind_transform() is called.
    //   Initializations which depend on nfreq and nds should go in _bind_transform(), not the 
    //   constructor.
    //
    //   Note that the value of nt_chunk is not necessarily initialized when _bind_transform()
    //   is called.  If you need the value of nt_chunk in _bind_transform(), you should either
    //   set it by hand, or call finalize_nt_chunk(), which initializes it automatically if it
    //   hasn't been initialized yet.
    //
    //    Optionally, _bind_transform() may read/create json attributes.  In CHIME, the attributes
    //    'freq_lo_MHz', 'freq_hi_MHz', 'dt_sample' should already exist.
    //
    // _unbind_transform()
    //
    //    Any code needed to undo _bind_transform() can go here.  (Usually not needed.)
    //
    // _process_chunk(intensity, istride, weights, wstride, pos)
    //
    //    This is the "core" method which is responsible for processing the 'intensity' and 'weights'
    //    arrays over timestamp range [pos,pos+nt_chunk).  The array strides are istride/wstride, 
    //    i.e. the (i,j)-th element of the intensity  array is intensity[i*istride+j], and similarly 
    //    for the weights.
    //
    //    Transforms which support downsampling (i.e. nds > 1) should note that 'pos'
    //    and 'nt_chunk' do not have the downsampling factor applied, but array dimensions
    //    do.  That is, the 'intensity' and 'weights' arrays have shape (nfreq, nt_chunk/nds),
    //    not shape (nfreq, nt_chunk), and 'pos' increases by nt_chunk (not nt_chunk/nds)
    //    in each call to _process_chunk();

    virtual void _bind_transform(Json::Value &json_attrs);  // non-pure virtual (default does nothing)
    virtual void _process_chunk(float *intensity, ssize_t istride, float *weights, ssize_t wstride, ssize_t pos) = 0;
    virtual void _unbind_transform();

    // prebind_nfreq, prebind_nds: saved values of nt_chunk, before it is finalized in bind().
    // This is useful in jsonize() and _unbind().

    // When _bindc() is called, it saves the values of 'nfreq' and 'nds' (before calling _bind_transform(), which can set them.)
    // This is useful in jsonize() and _unbind().

    ssize_t _prebind_nfreq = 0;  // intended to be accessed through get_prebind_nfreq
    ssize_t _prebind_nds = 0;    // intended to be accessed through get_prebind_nds

    // Helpers for jsonize(): originally-specified values of 'nfreq', 'nds' before bind() is called.
    // (Similar to chunked_pipeline_object::get_orig_nt_chunk(), which is inherited by wi_transform.)

    ssize_t get_prebind_nfreq() const { return (state >= BOUND) ? nfreq : _prebind_nfreq; }
    ssize_t get_prebind_nds() const   { return (state >= BOUND) ? nds : _prebind_nds; }

    // The wi_transform subclass shouldn't need to use these directly, 
    // since its _process_chunk() method operates directly on pointers/strides.
    
    std::shared_ptr<ring_buffer> rb_intensity;
    std::shared_ptr<ring_buffer> rb_weights;

    // We override chunked_pipeline::finalize_nt_chunk() and _check_nt_chunk(), in order to incorporate 'kernel_chunk_size'.
    virtual void finalize_nt_chunk() override;
    virtual void _check_nt_chunk() const override;

    // Subclass can optionally override: jsonize(), _allocate(), _deallocate(), _start_pipeline(), _end_pipeline(), _reset().
};


// -------------------------------------------------------------------------------------------------
//
// chime_file_stream_base
//
// The CHIME file-reading code is complex enough that it has its own base class!
// Eventually this may move to another library...
//
// Note: class chime_file_stream_base is currently not exported to python.


class chime_file_stream_base : public wi_stream {
public:
    chime_file_stream_base(const std::string &stream_name, const std::vector<std::string> &filename_list, ssize_t nt_chunk, ssize_t noise_source_align);
    virtual ~chime_file_stream_base() { }

    // Throughout the CHIMEFRB backend, we represent times in seconds, but the raw packets use timestamps
    // constructed from FPGA counts.  We convert by assuming that each FPGA count is exactly 2.56e-6 seconds.
    // The precise conversion matters (to machine precision!) when predicting the location of the noise
    // source edges from the timestamps.  Therefore, the 2.56e-6 "magic number" must be used consistently
    // throughout rf_pipelines and ch_vdif_assembler.

    static constexpr double chime_seconds_per_fpga_count = 2.56e-6;

protected:
    // Specified at construction.  For an explanation of the 'noise_source_align' field see rf_pipelines.hpp.
    // Note that the 'nt_chunk' constructor argument is used to initialize the base class member wi_stream::nt_maxwrite.
    const std::vector<std::string> filename_list;
    const ssize_t noise_source_align;   // if zero, no alignment will be performed

    int curr_ifile = -1;  // index of current file in filename_list

    // Only nonzero if noise source alignment is requested
    ssize_t initial_discard_count = 0;

    // Putting these here for now, they will move into a json object soon...
    double freq_lo_MHz = 0.0;
    double freq_hi_MHz = 0.0;
    double dt_sample = 0.0;

    // Parameters from current file.
    double time_lo = 0.0;
    double time_hi = 0.0;
    ssize_t nt_file = 0;
    bool frequencies_are_increasing = false;

    // Time index in current file.  Can be negative!  This means there is a time gap between files.
    ssize_t it_file = 0;

    // Devirtualize wi_stream base class.
    virtual void _bind_stream(Json::Value &json_attrs) override;
    virtual bool _fill_chunk(float *intensity, ssize_t istride, float *weights, ssize_t wstride, ssize_t pos) override;
    virtual void _end_pipeline(Json::Value &json_output) override;
    virtual void _unbind_stream() override;

    // Pure virtuals which follow must be defined by subclass!

    // Reads file from disk into a subclass-specific internal data structure.
    // This will be followed by calls to set_params_from_file() and/or check_file_consistency().
    virtual void load_file(const std::string& filename) = 0;

    // Responsible for initializing wi_stream members: { nfreq, freq_lo_MHz, freq_hi_MHz, dt_sample }
    // and chime_file_stream_base members: { time_lo, time_hi, nt, frequencies_are_increasing }
    virtual void set_params_from_file() = 0;
    
    // Checks consistency between file and data members mentioned above, but doesn't initialize them.
    virtual void check_file_consistency() const = 0;

    // Reads a shape-(nfreq,n) block from current file.
    // The 'dst_int' and 'dst_wt' arrays have shape (nfreq, n) and strides 'dst_istride', 'dst_wstride'.
    // The 'it_file' argument is the initial timestamp of the block, relative to the start of the current file.
    //
    // The subclass can choose whether the 'dst_int' and 'dst_wt' arrays are indexed from lowest frequency 
    // channel to highest, or from highest frequency channel to lowest (the rf_pipelines_default).  This is 
    // done by setting the 'frequencies_are_increasing' flag in set_params_from_file().  The way it is 
    // implemented is via logic in the base class which flips the sign of 'dst_stride' if necessary.

    virtual void read_data(float* dst_int, float* dst_wt, ssize_t it_file, ssize_t n, ssize_t dst_istride, ssize_t dst_wstride) const = 0;

    virtual void close_file() = 0;
};


}  // namespace rf_pipelines

#endif // _RF_PIPELINES_BASE_CLASSES_HPP
