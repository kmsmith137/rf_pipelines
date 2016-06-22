#include "rf_pipelines_internals.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif


wi_stream::wi_stream() :
    nfreq(0), freq_lo_MHz(0.), freq_hi_MHz(0.), dt_sample(0.), nt_maxwrite(0)
{ }


wi_stream::wi_stream(int nfreq_, double freq_lo_MHz_, double freq_hi_MHz_, double dt_sample_, int nt_maxwrite_) :
    nfreq(nfreq_), freq_lo_MHz(freq_lo_MHz_), freq_hi_MHz(freq_hi_MHz_), dt_sample(dt_sample_), nt_maxwrite(nt_maxwrite_)
{
    this->check_invariants();
}


void wi_stream::check_invariants() const
{
    if (nfreq <= 0)
	throw runtime_error("wi_stream: nfreq is non-positive or uninitialized");
    if (freq_lo_MHz <= 0.0)
	throw runtime_error("wi_stream: frequency range is invalid or uninitialized");
    if (freq_lo_MHz >= freq_hi_MHz)
	throw runtime_error("wi_stream: frequency range is invalid or uninitialized");
    if (dt_sample <= 0.0)
	throw runtime_error("wi_stream: dt_sample is non-positive or uninitialized");	
    if (nt_maxwrite <= 0)
	throw runtime_error("wi_stream: nt_maxwrite is non-positive or uninitialized");
}


void wi_stream::run_transforms(const std::vector<std::shared_ptr<wi_transform> > &transforms)
{
    this->check_invariants();

    // Delegate to run_stream() method implemented in subclass
    wi_run_state rstate(*this, transforms);
    this->run_stream(rstate);

    if (rstate.is_running)
	throw runtime_error("rf_pipelines: stream is still in \"running\" state after run_stream() returns (maybe you forgot to call wi_run_state::end_stream() in wi_stream::run_stream()?)");
}


}  // namespace rf_pipelines
