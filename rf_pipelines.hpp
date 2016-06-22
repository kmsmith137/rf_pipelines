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
class wi_run_state;


struct wi_stream {
    int nfreq;
    double freq_lo_MHz;
    double freq_hi_MHz;
    double dt_sample;

    // Max number of time samples per chunk emitted by stream
    // (Shouldn't be too large, since it determines an internal buffer size)
    int nt_maxwrite;
    
    // If the default constructor is used, then the subclass is responsible for initializing the
    // fields { nfreq, freq_lo_MHz, freq_hi_MHz, dt_sample, nt_mawrite } before the first call to run().
    wi_stream();
    wi_stream(int nfreq, double freq_lo_MHz, double freq_hi_MHz, double dt_sample, int nt_maxwrite);

    virtual ~wi_stream() { }

    // Helper function
    void check_invariants() const;

    // run_transforms() is intended to be the main interface for running a sequence of transforms on a stream
    void run_transforms(const std::vector<std::shared_ptr<wi_transform> > &transforms);
    
    //
    // Schematically, the run_stream() method should look something like this:
    //
    //    rstate.start_stream()
    //    while (cond) {
    //        rstate.setup_write();
    //        rstate.finalize_write();
    //    }
    //    rstate.end_stream()
    //
    virtual void run_stream(wi_run_state &rstate) = 0;
};


struct wi_transform {
    int nt_chunk;
    int nt_prepad;
    int nt_postpad;

    // If the default constructor is used, then the subclass is responsible for initializing the
    // fields { nt_chunk, nt_prepad, nt_postpad } either in the constructor or start_stream().
    explicit wi_transform(int nt_chunk, int nt_prepad=0, int nt_postpad=0);
    wi_transform();

    virtual ~wi_transform() { }

    // Helper function
    void check_invariants() const;

    // Note: OK if call to set_nt() doesn't get called until initialize_transform().
    virtual void start_stream(int nfreq, double freq_lo_MHz, double freq_hi_MHz, double dt_sample) = 0;
    virtual void process_chunk(float *intensity, float *weight, int stride, float *pp_intensity, float *pp_weight, int pp_stride) = 0;
    virtual void end_stream() = 0;
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
    int nt_ring;

    // 2d arrays of shape (nfreq, nt_tot)
    std::vector<float> intensity;
    std::vector<float> weights;
    int nt_tot;

    int ipos;

    // Main constructor syntax
    wraparound_buf(int nfreq, int nt_contig, int nt_ring);

    // Alternate syntax: use default constuctor, then call construct()
    wraparound_buf();

    void construct(int nfreq, int nt_contig, int nt_ring);
    void reset();

    void setup_write(int it0, int nt, float* &intensityp, float* &weightp, int &stride);
    void setup_append(int nt, float* &intensityp, float* &weightp, int &stride, bool zero_flag);
    void append_zeros(int nt);

    void finalize_write(int it0, int nt);
    void finalize_append(int nt);

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

    bool is_running;
};



}  // namespace rf_pipelines

#endif // _RF_PIPELINES_HPP
