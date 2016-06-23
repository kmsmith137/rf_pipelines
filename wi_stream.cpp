#include "rf_pipelines_internals.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif


void wi_stream::run(const std::vector<std::shared_ptr<wi_transform> > &transforms)
{
    if (transforms.size() == 0)
	throw runtime_error("wi_stream::run() called on empty transform list");

    // Lots of checks to make sure everything is initialized

    if (nfreq <= 0)
	throw runtime_error("wi_stream::nfreq is non-positive or uninitialized");
    if (freq_lo_MHz <= 0.0)
	throw runtime_error("wi_stream::frequency range is invalid or uninitialized");
    if (freq_lo_MHz >= freq_hi_MHz)
	throw runtime_error("wi_stream::frequency range is invalid or uninitialized");
    if (dt_sample <= 0.0)
	throw runtime_error("wi_stream::dt_sample is non-positive or uninitialized");	
    if (nt_maxwrite <= 0)
	throw runtime_error("wi_stream::nt_maxwrite is non-positive or uninitialized");

    for (unsigned int it = 0; it < transforms.size(); it++) {
	transforms[it]->set_stream(*this);

	if (transforms[it]->nfreq != this->nfreq)
	    throw runtime_error("wi_transform::nt_chunk does not match stream nfreq");
	if (transforms[it]->nt_chunk <= 0)
	    throw runtime_error("wi_transform::nt_chunk is non-positive or uninitialized");
	if (transforms[it]->nt_prepad < 0)
	    throw runtime_error("wi_transform::nt_prepad is negative");
	if (transforms[it]->nt_postpad < 0)
	    throw runtime_error("wi_transform::nt_postpad is negative");
    }

    // Delegate to stream_body() method implemented in subclass
    wi_run_state run_state(*this, transforms);
    this->stream_body(run_state);

    // Only state=4 is OK, otherwise we complain!
    if (run_state.state == 0)
	throw runtime_error("rf_pipelines: logic error in stream: run() returned without doing anything");
    if (run_state.state == 1)
	throw runtime_error("rf_pipelines: logic error in stream: run() returned after calling set_substream(), without writing data or calling end_substream()");
    if (run_state.state == 2)
	throw runtime_error("rf_pipelines: logic error in stream: run() returned after calling setup_write(), without calling finalize_write() or end_substream()");
    if (run_state.state == 3)
	throw runtime_error("rf_pipelines: logic error in stream: run() returned without calling end_substream()");
}


}  // namespace rf_pipelines
