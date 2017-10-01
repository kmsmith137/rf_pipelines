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


// Defined in rf_pipelines_internals.hpp
struct outdir_manager;
struct plot_group;


// -------------------------------------------------------------------------------------------------
//
// ring_buffer, ring_buffer_dict, ring_buffer_subarray


class ring_buffer {
public:
    // "Complementary" dimensions (all dimensions except time axis)
    const std::vector<ssize_t> cdims;
    const ssize_t csize;  // product of all cdims
    
    // Downsampling factor
    const ssize_t nds;                 

    ring_buffer(const std::vector<ssize_t> &cdims, ssize_t nds);

    void update_params(ssize_t nt_contig, ssize_t nt_maxlag);

    void allocate();
    void deallocate();
    void start();
    
    // The access_mode is optional, but enables some debug checks, and can
    // also help performance.  The numerical values are chosen for convenient
    // bitmasking.  The ACCESS_NONE value is a placeholder which throws an exception.

    static constexpr int ACCESS_NONE = 0;
    static constexpr int ACCESS_READ = 0x01;
    static constexpr int ACCESS_WRITE = 0x02;
    static constexpr int ACCESS_RW = 0x03;
    static constexpr int ACCESS_APPEND = 0x06;
    
    float *get(ssize_t pos0, ssize_t pos1, int access_mode);
    void put(float *p, ssize_t pos0, ssize_t pos1, int access_mode);

    ssize_t get_stride() const;

    static std::string access_mode_to_string(int access_mode);

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
    ssize_t curr_pos = 0;
    ssize_t first_valid_sample = 0;
    ssize_t last_valid_sample = 0;

    // Is there an active pointer?
    float *ap = nullptr;
    float ap_pos0 = 0;
    float ap_pos1 = 0;
    float ap_mode = ACCESS_NONE;

    // Helper functions
    void _mirror_initial(ssize_t it0);
    void _mirror_final(ssize_t it1);
    void _copy(ssize_t it_dst, ssize_t it_src, ssize_t n);

    friend struct ring_buffer_subarray;
};


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

    // For now, this is a lightweight class designed to live briefly on the call stack, with no copying.
    ring_buffer_subarray(const ring_buffer_subarray &) = delete;
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
// Some abstract base classes: pipeline_object, chunked_pipeline_object, wi_stream, wi_transform.


struct pipeline_object {
public:
    // We now require the name to be initialized at construction (usually to something simple like the class name).
    // It can be changed later (in the subclass constructor, or in _bind()) to something more verbose!
    std::string name;
    
    // These parameters control the flow of data into the pipeline_object.
    // They are set "externally", just before the virtual function _bind() is called.
    // Don't initialize or change them in the pipeline_object subclass, or strange things will happen!

    ssize_t nt_chunk_in = 0;   // step size in pos_hi
    ssize_t nt_maxlag = 0;     // max difference between pos_hi and pos_max, when _advance() is called

    // These parameters control the flow of data out of the pipeline object.
    // The subclass must initialize them, in the virtual function _bind().
    // Don't change them after _bind() returns, or strange things will happen!
    
    ssize_t nt_chunk_out = 0;  // step size in pos_lo
    ssize_t nt_maxgap = -1;    // max allowed value of (pos_hi-pos_lo), after _advance() returns
    ssize_t nt_contig = 0;     // max contiguous chunk size requested from ring buffers
    
    // Runtime state.
    ssize_t pos_lo = 0;   // always a multiple of nt_chunk_out
    ssize_t pos_hi = 0;   // always a multiple of nt_chunk_in
    ssize_t pos_max = 0;

    // Constructor for this abstract base class.
    explicit pipeline_object(const std::string &name);

    virtual ~pipeline_object() { }

    // High-level API: to run a pipeline, just call run().
    //
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
    // The meaning of the 'verbosity' argument is:
    //   0 = no output
    //   1 = high-level summary output (names of transforms, number of samples processed etc.)
    //   2 = show all output files
    //   3 = debug trace through pipeline

    Json::Value run(const std::string &outdir=".", int verbosity=2, bool clobber=true);

    // A more fine-grained high-level API.
    // bind() is the first step in pipeline: determines pipeline parameters such as ring buffer sizes
    void bind();  
    void allocate();
    void deallocate();    
    bool is_bound() const;

    // By default, this virtual function throws an exception ("jsonize() unimplemented...")
    // Objects which support jsonization should override this function, and also add
    // deserialization code to pipeline_object::from_json().

    virtual Json::Value jsonize() const;

    // A pipeline_object subclass which implements jsonize() will also want to implement a
    // from_json() static member function, and "register" it with register_json_constructor().
    //
    // The 'f' argument to register_json_converter() should be a function f(x) whose single
    // argument 'x' is a Json::Value, and returns a shared_ptr<pipeline_ohbject>.  By convention,
    // f() is usually chosen to be a static member function 'from_json' of the pipeline_object
    // subclass.  See pipeline.cpp for an example.

    using json_constructor_t = std::function<std::shared_ptr<pipeline_object> (const Json::Value &)>;
    static void register_json_constructor(const std::string &class_name, const json_constructor_t &f);

    // The implementation of from_json() in the pipeline_object base class will search the
    // registry for the matching class_name, and call the registered constructor.

    static std::shared_ptr<pipeline_object> from_json(const Json::Value &j);

    // You probably don't want to call the functions below, but they need to be public
    // (I think) in order to implement container-like subclasses of pipeline_object,
    // such as 'class pipeline' or 'class wi_sub_pipeline'.
    //
    // A note on the virtual function get_preferred_chunk_size().  This is only called
    // on the first pipeline_object in the pipeline (subsequent chunk sizes are determined
    // automatically).  By default, this function returns 0, which results in an
    // exception "...: this pipeline_object cannot appear first in pipeline".
    //
    // Stream-type pipeline_objects which can appear first in a pipeline should override
    // get_preferred_chunk_size() to return a nonzero value.
    
    void bind(ring_buffer_dict &rb_dict, ssize_t nt_chunk_in, ssize_t nt_maxlag, Json::Value &json_data);
    ssize_t advance(ssize_t pos_hi, ssize_t pos_max);
    
    void start_pipeline(const std::shared_ptr<outdir_manager> &mp, Json::Value &j);
    void end_pipeline(Json::Value &j);
    
    virtual ssize_t get_preferred_chunk_size();

    // Disable copy constructors
    pipeline_object(const pipeline_object &) = delete;
    pipeline_object &operator=(const pipeline_object &) = delete;

// FIXME 'protected' removed temporarily, until I figure out how to handle it in the python-wrapping code!
// protected:

    // Helper functions called by _bind().
    std::shared_ptr<ring_buffer> get_buffer(ring_buffer_dict &rb_dict, const std::string &key);
    std::shared_ptr<ring_buffer> create_buffer(ring_buffer_dict &rb_dict, const std::string &key, const std::vector<ssize_t> &cdims, ssize_t nds);

    // These helper functions are used by pipeline_objects which write output files (e.g. hdf5, png).
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

    // Pure virtuals
    //
    // _bind()
    //    - initializes nt_chunk_out, nt_contig, nt_maxgap
    //    - calls get_buffer() or create_buffer() for each ring buffer used by pipeline_object
    //
    //    - The (Json::Value &) argument contains "pipeline-global" parameters.  Each pipeline_object
    //      can either create new parameters, or read existing ones.  In CHIME, the pipeline-global
    //      parameters are: 'freq_lo_MHz', 'freq_hi_MHz', 'dt_sample', and are 
    //
    // _allocate()
    //    - should do nothing if already allocated.
    //
    // _start_pipeline()
    //
    //    - Note that start_pipeline() has a shared_ptr<outdir_manager> argument, but _start_pipeline()
    //      does not.  In _start_pipeline(), the outdir_manager is available as 'this->outdir_manager'.
    //
    //    - Container pipeline_objects must call 
    //         p.start_pipeline(this->outdir_manager, json_data)
    //      for each sub-object 'p'.
    //
    //    - Both _start_pipeline() and _end_pipeline() take a (Json::Value &) argument, but the meaning
    //      is different.  In _start_pipeline(), the Json::Value contains "pipeline-global" parameters
    //      (just like _bind()).  In _end_pipeline(), the Json::Value contains pipeline output, which
    //      each pipeline_object is free to define independently.
    //
    //    - Some pipeline-global parameters are passed to _bind(), and others are passed to _start_pipeline().
    //      The division is a little arbitrary, but we try to use _bind() if possible.  In CHIME, the only
    //      _start_pipeline() parameter is 'initial_fpga_count', which isn't known until the first packet
    //      is received.  As a consequence, _start_pipeline() should be very fast (e.g. avoid allocating
    //      memory) since it's part of the real-time pipeline.
    //
    // _advance()
    //
    //    From the perspective of a single pipeline_object, the pipeline run logic is as follows.
    //
    //    An external caller advances all input_buffers, and the index 'pos_hi', by a multiple of 
    //    nt_chunk_in, then calls _advance().  The pipeline_object is then responsible for advancing
    //    'pos_lo' by a multiple of nt_chunk_out, so that pos_lo >= pos_hi - nt_maxlag.
    //
    //    The return value from _advance() is a new nt_end.  Thus:
    //       - return SSIZE_MAX to continue processing.
    //       - return 0 to "hard-terminate" the pipeline, without throwing an exception.
    //       - if end-of-stream is encountered, return the number of samples in the stream.
    //
    // _end_pipeline()
    //
    //    - Both _start_pipeline() and _end_pipeline() take a (Json::Value &) argument, but the meaning
    //      is different.  In _start_pipeline(), the Json::Value contains "pipeline-global" parameters
    //      (just like _bind()).  In _end_pipeline(), the Json::Value contains pipeline output, which
    //      each pipeline_object is free to define independently.
    //
    //    - In addition to any subclass-dependent outputs defined in _end_pipeline(), the pipeline
    //      defines the following default json outputs: "name", "cpu_time", "plots".  For many pipeline_objects,
    //      these default outputs are sufficient, and defining _end_pipeline() isn't necessary.

    virtual void _bind(ring_buffer_dict &rb_dict, Json::Value &json_data) = 0;    
    virtual ssize_t _advance() = 0;

    // By default, these virtual functions do nothing.
    virtual void _allocate();
    virtual void _deallocate();
    virtual void _start_pipeline(Json::Value &j);
    virtual void _end_pipeline(Json::Value &j);

    // Helper: throws runtime_error with prefix "rf_pipelines: <name>: ..."
    void _throw(const std::string &msg) const;

    // For debugging or internal use
    static void _show_registered_json_constructors();
    static json_constructor_t _find_json_constructor(const std::string &class_name);  // can return NULL.

    // Internal state, managed automatically by calls to create_ring_buffer(), add_plot_group(), etc.

    // These fields are initialized in bind().
    std::vector<std::shared_ptr<ring_buffer>> new_ring_buffers;    // ring buffers created by this pipeline object
    std::vector<std::shared_ptr<ring_buffer>> all_ring_buffers;    // all ring buffers used by this pipeline object (including new_ring_buffers)

    // These fields are initialized in start_pipeline() and cleared in end_pipeline().
    // Note: if the 'outdir' argument to run() is an empty string (or python None), then
    // 'outdir_manager' will be a nonempty pointer, but calls to outdir_manager->add_file() will fail.

    std::shared_ptr<outdir_manager> out_mp;
    std::vector<plot_group> plot_groups;
    double time_spent_in_transform = 0.0;
};


struct chunked_pipeline_object : public pipeline_object {
public:
    const bool can_be_first;

    // If nt_chunk is zero, then it will be initialized (in bind()) to match the previous transform.
    // However, if 'can_be_first' is true, then nt_chunk must be initialized to a nonzero value.
    // Note that nt_chunk can be initialized anywhere in the subclass constructor, but not after the constructor exits.
    ssize_t nt_chunk = 0;

    // 'name' must be initialized at construction, but can be changed later (in subclass constructor, or _bind_chunked()).
    chunked_pipeline_object(const std::string &name, bool can_be_first, ssize_t nt_chunk=0);

// FIXME 'protected' removed temporarily, until I figure out how to handle it in the python-wrapping code!
// protected:

    // These pipeline_object virtuals are defined here.
    virtual void _bind(ring_buffer_dict &rb_dict, Json::Value &json_data) override;
    virtual ssize_t _advance() override;
    virtual ssize_t get_preferred_chunk_size() override;

    // New pure virtuals, to be defined by subclass.
    virtual void _bind_chunked(ring_buffer_dict &rb_dict, Json::Value &json_data) = 0;
    virtual bool _process_chunk(ssize_t pos) = 0;

    // Helper for jsonize(): this is the originally-specified value of 'nt_chunk' before bind() is called.
    // This can differ from 'nt_chunk', if nt_chunk is originally zero, and bind() initializes it to a nonzero value.
    ssize_t get_orig_nt_chunk() const;
    ssize_t _save_nt_chunk = 0;

    // Subclass can optionally override: jsonize(), _allocate(), _deallocate(), _start_pipeline(), _end_pipeline().
};


struct wi_stream : chunked_pipeline_object {
public:
    // Note: inherits 'name', 'nt_chunk' from base classes.

    // 'nfreq' and 'nt_chunk' must be initialized to nonzero values,
    // but this can be done any time before bind() is called.
    ssize_t nfreq = 0;

    std::shared_ptr<ring_buffer> rb_intensity;
    std::shared_ptr<ring_buffer> rb_weights;

    // 'name' must be initialized at construction, but can be changed later (in subclass constructor, or _bind_stream()).
    wi_stream(const std::string &name, ssize_t nfreq=0, ssize_t nt_chunk=0);

// FIXME 'protected' removed temporarily, until I figure out how to handle it in the python-wrapping code!
// protected:

    // These chunked pipeline_object virtuals are defined here.
    virtual void _bind_chunked(ring_buffer_dict &rb_dict, Json::Value &json_data) override;
    virtual bool _process_chunk(ssize_t pos) override;

    // New virtuals, to be defined by subclass.
    // Sometimes it's convenient to defer initialization of nfreq and nt_chunk to _bind_stream().
    // The _fill_chunk() return value should be 'false' if EOF occurred somewhere in the chunk.
    virtual void _bind_stream(Json::Value &json_data);  // non-pure virtual (default does nothing)
    virtual bool _fill_chunk(float *intensity, ssize_t istride, float *weights, ssize_t wstride, ssize_t pos) = 0;

    // Subclass can optionally override: jsonize(), _allocate(), _deallocate(), _start_pipeline(), _end_pipeline().
};


struct wi_transform : chunked_pipeline_object {
public:
    // Note: inherits 'name', 'nt_chunk' from base classes.
    
    std::shared_ptr<ring_buffer> rb_intensity;
    std::shared_ptr<ring_buffer> rb_weights;

    // If (nfreq, nds) are not initialized in the constructor, then they will
    // be initialized in bind().
    //
    // If (nfreq, nds) are initialized to nonzero values in the constructor,
    // then bind() will throw an exception if there is a mismatch with the
    // true (nfreq, nds) of the data.
  
    ssize_t nfreq = 0;
    ssize_t nds = 0;

    // 'name' must be initialized at construction, but can be changed later (in subclass constructor, or _bind_stream()).
    wi_transform(const std::string &name, ssize_t nt_chunk=0, ssize_t nfreq=0, ssize_t nds=0);

// FIXME 'protected' removed temporarily, until I figure out how to handle it in the python-wrapping code!
// protected:

    // These chunked_pipeline_object virtuals are defined here.
    virtual void _bind_chunked(ring_buffer_dict &rb_dict, Json::Value &json_data) override;
    virtual bool _process_chunk(ssize_t pos) override;

    // New virtuals, to be defined by subclass.
    // Note that _bind_transform() is called after 'nfreq' and 'nds' get initialized, so subclass-dependent
    // initializations which depend on their values should go there.
    virtual void _bind_transform(Json::Value &json_data);  // non-pure virtual (default does nothing)
    virtual void _process_chunk(float *intensity, ssize_t istride, float *weights, ssize_t wstride, ssize_t pos) = 0;

    // Helpers for jsonize(): originally-specified values of 'nfreq', 'nds' before bind() is called.
    // (Similar to chunked_pipeline_object::get_orig_nt_chunk(), which is inherited by wi_transform.)
    ssize_t get_orig_nfreq() const;
    ssize_t get_orig_nds() const;
    ssize_t _save_nfreq = 0;
    ssize_t _save_nds = 0;

    // Subclass can optionally override: jsonize(), _allocate(), _deallocate(), _start_pipeline(), _end_pipeline().
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

    // Devirtualize wi_stream::_fill_chunk()
    virtual void _bind_stream(Json::Value &json_data) override;
    virtual bool _fill_chunk(float *intensity, ssize_t istride, float *weights, ssize_t wstride, ssize_t pos) override;

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
