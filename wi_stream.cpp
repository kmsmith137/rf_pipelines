// Note: I haven't systematically documented the C++ interface to rf_pipelines,
// so the level of documentation will be hit-or-miss.  Also please note that the
// python-wrapping in rf_pipelines_c.cpp is kind of a mess which I hope to improve
// soon.  In the meantime if you want to python-wrap a C++ class, just email me
// and I'll help navigate the mess!

#include "rf_pipelines_internals.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif

//
// A helper class which resets wi_transform::outdir_manager pointers in its destructor.
// This is convenient for ensuring that the pointer always gets reset, e.g. in the case where an
// exception is thrown.  We also use this class to detect transform reuse (currently treated as an error).
//
struct outdir_janitor {
    shared_ptr<outdir_manager> manager;
    vector<shared_ptr<wi_transform> > transform_list;

    outdir_janitor(const string &outdir, bool clobber)
    {
	this->manager = make_shared<outdir_manager> (outdir, clobber);
    }

    // noncopyable
    outdir_janitor(const outdir_janitor &) = delete;
    outdir_janitor &operator=(const outdir_janitor &) = delete;

    ~outdir_janitor() 
    { 
	for (unsigned int i = 0; i < transform_list.size(); i++)
	    transform_list[i]->outdir_manager.reset(); 
	transform_list.clear();
    }

    void set_outdir_manager(const shared_ptr<wi_transform> &transform)
    {
	rf_assert(transform);

	if (transform->outdir_manager) {
	    string name = transform->get_name();
	    throw runtime_error("Fatal: transform->outdir_manager initialized twice.  This probably means a transform is being reused (transform name=" + name + ")");
	}

	transform->outdir_manager = manager;
	this->transform_list.push_back(transform);
    }
};


void wi_stream::run(const vector<shared_ptr<wi_transform> > &transforms, const string &outdir, bool noisy, bool clobber)
{
    int ntransforms = transforms.size();

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
    if (ntransforms == 0)
	throw runtime_error("wi_stream::run() called on empty transform list");

    for (const auto &transform: transforms)
	if (!transform)
	    throw runtime_error("rf_pipelines: empty transform pointer passed to wi_stream::run()");

    outdir_janitor janitor(outdir, clobber);

    for (int it = 0; it < ntransforms; it++) {
	janitor.set_outdir_manager(transforms[it]);
	
	transforms[it]->set_stream(*this);

	if (transforms[it]->nfreq != this->nfreq)
	    throw runtime_error("rf_pipelines: transform's value of 'nfreq' does not match stream's value of 'nfreq'");
	if (transforms[it]->nt_chunk <= 0)
	    throw runtime_error("rf_pipelines: transform's value of 'nt_chunk' is non-positive or uninitialized");
	if (transforms[it]->nt_prepad < 0)
	    throw runtime_error("rf_pipelines: wi_transform::nt_prepad is negative");
	if (transforms[it]->nt_postpad < 0)
	    throw runtime_error("rf_pipelines: wi_transform::nt_postpad is negative");
    }

    // Delegate to stream_body() method implemented in subclass
    wi_run_state run_state(*this, transforms, janitor.manager, noisy);
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
