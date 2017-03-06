#include "rf_pipelines_internals.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
};  // pacify emacs c-mode
#endif

// The 'bb_dedisperser' wrapper class converts Ben Barsdell's 'dedisp' GPU code
// to the rf_pipelines API.  It does this by subclassing rf_pipelines::wi_transform,
// and defining the approporiate virtual functions.

struct bb_dedisperser : public wi_transform {
    // Constructor arguments
    double dm_start;
    double dm_end;
    double dm_tol;
    double pulse_width_ms;

    // Additional class members (for example, the dedisp_plan object) can be added here.

    bb_dedisperser(double dm_start_, double dm_end_, double dm_tol_, double pulse_width_ms_);

    // The behavior of bb_dedisperser (or any subclass of wi_transform) is defined by
    // implementing these virtual functions.  For documentation on what needs to be
    // implemented, see comments in rf_pipelines.hpp.

    virtual void set_stream(const wi_stream &stream) override;
    virtual void start_substream(int isubstream, double t0) override;
    virtual void process_chunk(double t0, double t1, float *intensity, float *weights, ssize_t stride, float *pp_intensity, float *pp_weights, ssize_t pp_stride) override;
    virtual void end_substream() override;
};


bb_dedisperser::bb_dedisperser(double dm_start_, double dm_end_, double dm_tol_, double pulse_width_ms_) :
    dm_start(dm_start_),
    dm_end(dm_end_),
    dm_tol(dm_tol_),
    pulse_width_ms(pulse_width_ms_)
{ 
    // Initialize this->name
    stringstream ss;
    ss << "bb_dedisperser(dm_start=" << dm_start << ",dm_end=" << dm_end
       << ",dm_tol=" << dm_tol << ",pulse_width_ms=" << pulse_width_ms << ")";

    this->name = ss.str();
}

// Placeholders for virtual functions follow.

void bb_dedisperser::set_stream(const wi_stream &stream)
{
    throw std::runtime_error("bb_dedisperser::set_stream() called but not implemented yet");
}

void bb_dedisperser::start_substream(int isubstream, double t0)
{
    throw std::runtime_error("bb_dedisperser::start_substream() called but not implemented yet");
}

void bb_dedisperser::process_chunk(double t0, double t1, float *intensity, float *weights, ssize_t stride, float *pp_intensity, float *pp_weights, ssize_t pp_stride)
{
    throw std::runtime_error("bb_dedisperser::process_chunk() called but not implemented yet");
}

void bb_dedisperser::end_substream()
{
    throw std::runtime_error("bb_dedisperser::end_substream() called but not implemented yet");

    // After the dedisperser has run to completion, you should have values of the FRB DM,
    // arrival time, and signal-to-noise.  In general, rf_pipelines transforms communicate
    // their outputs by setting fields in their JSON output (this will end up in a pipeline
    // output file named something like 'rf_pipelines_0.json').  The syntax below shows
    // how to do this!

    // Placeholder values
    double sn = 30.0;
    double dm = 100.0;
    double tarr = 2.0;

    // Initialize fields in JSON output.
    this->json_per_substream["frb_sn"] = sn;
    this->json_per_substream["frb_dm"] = dm;
    this->json_per_substream["frb_arrival_time"] = tarr;
}


// -------------------------------------------------------------------------------------------------
//
// make_bb_dedisperser()
//
// This interface (function which returns a shared pointer to the wi_transform base class) 
// allows the implementation of struct bb_dedisperser to be "hidden" from the rest of rf_pipelines.


std::shared_ptr<wi_transform> make_bb_dedisperser(double dm_start, double dm_end, double dm_tol, double pulse_width_ms)
{
    return std::make_shared<bb_dedisperser> (dm_start, dm_end, dm_tol, pulse_width_ms);
}


}  // namespace rf_pipelines
