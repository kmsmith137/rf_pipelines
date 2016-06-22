#ifndef _RF_PIPELINES_HPP
#define _RF_PIPELINES_HPP

#if (__cplusplus < 201103) && !defined(__GXX_EXPERIMENTAL_CXX0X__)
#error "This source file needs to be compiled with C++0x support (g++ -std=c++0x)"
#endif

#include <vector>
#include <memory>

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif

struct wi_transform;
struct wi_run_state;


struct wi_stream {
    int nfreq;
    double freq_lo_MHz;
    double freq_hi_MHz;
    double dt_sample;

    // Max number of time samples per chunk emitted by stream
    // (Shouldn't be too large, since this determines an internal buffer size)
    int nt_maxwrite;

    // Subclass constructor can either call this base class constructor immediately...
    wi_stream(int nfreq, double freq_lo_MHz, double freq_hi_MHz, double dt_sample, int nt_maxwrite);

    // ...or use a fake default constructor, followed later by a call to construt()
    wi_stream();
    void construct(int nfreq, double freq_lo_MHz, double freq_hi_MHz, double dt_sample, int nt_maxwrite);

    // This run() method is intended to be the main interface for running a sequence of transforms on a stream
    void run(const std::vector<std::shared_ptr<wi_transform> > &transforms);
    
    //
    // The virtual run() method implemented by the subclass has a different interface which operates
    // on a wi_run_state (a helper class; see below).
    //
    // It should call chain.setup_write() and chain.finalize_write() in a loop, and return when done.
    //
    virtual void run(const wi_run_state &rstate) = 0;

    virtual ~wi_stream() { }
};


struct wi_transform {
    int nt_chunk;
    int nt_prepad;
    int nt_postpad;

    // Subclass constructor can either call this base class constructor immediately...
    explicit wi_transform(int nt_chunk, int nt_prepad=0, int nt_postpad=0);

    // ...or use a fake default constructor, followed later by a call to set_nt()
    wi_transform();

    void set_nt(int nt_chunk, int nt_prepad=0, int nt_postpad=0);
    void reset_nt();

    // Note: OK if call to set_nt() doesn't get called until initialize_transform().
    virtual void start_stream(int nfreq, double freq_lo_MHz, double freq_hi_MHz, double dt_sample) = 0;
    virtual void process_chunk(float *data, float *weight, int stride, float *pp_data, float *pp_weight, int pp_stride) = 0;
    virtual void end_stream() = 0;

    virtual ~wi_transform() { }
};


// -------------------------------------------------------------------------------------------------
//
// Helper classes which are probably not needed from the "outside world"
//
// Exception: implementations of 
//


struct wraparound_buf {
    // specified at construction
    int nfreq;
    int nt_contig;
    int nt_logical_size;

    // 2d arrays of shape (nfreq, nt_tot)
    std::vector<float> intensity;
    std::vector<float> weights;
    int nt_tot;

    int ipos;

    // Main constructor syntax
    wraparound_buf(int nfreq, int nt_contig, int nt_logical_size);

    // Alternate syntax: use default constuctor, then call construct()
    wraparound_buf();

    void construct(int nfreq, int nt_contig, int nt_logical_size);
    void reset();

    void setup_read(int it0, int nt, float* &intensityp, float* &weightp, int &stride);
    void setup_write(int it0, int nt, float* &intensityp, float* &weightp, int &stride, bool zero_flag);
    void finalize_write(int it0, int nt);

    void _copy(int it_dst, int it_src, int nt);
    void _check_integrity();

    static void run_unit_tests();
};


class wi_run_state {
protected:
    // make noncopyable
    wi_run_state(const wi_run_state &) = delete;
    wi_run_state& operator=(const wi_run_state &) = delete;

    // stream data
    const int nfreq;
    const double freq_lo_MHz;
    const double freq_hi_MHz;
    const double dt_sample;
    const int nt_stream_maxwrite;
    
    // transform data
    const int ntransforms;
    const std::vector<std::shared_ptr<wi_transform> > transforms;

    // timeline
    std::vector<int> transform_ipos;   // satisfies transform_ipos[0] >= transform_ipos[1] >= ...
    int stream_ipos;
    bool is_running;

    // buffers
    wraparound_buf main_buffer;
    std::vector<wraparound_buf> prepad_buffers;

public:
    wi_run_state(const wi_stream &stream, const std::vector<std::shared_ptr<wi_transform> > &transforms);

    // Called by wi_stream::run()
    void start_stream();
    void setup_write(int nt, float* &intensityp, float* &weightp, int &stride, bool zero_flag);
    void finalize_write(int nt);
    void end_stream();
};



}  // namespace rf_pipelines

#endif // _RF_PIPELINES_HPP
