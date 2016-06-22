#include "rf_pipelines_internals.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif


wi_stream::wi_stream()
{
    this->nfreq = 0;
    this->freq_lo_MHz = 0.0;
    this->freq_hi_MHz = 0.0;
    this->dt_sample = 0.0;
    this->nt_maxwrite = 0;
}


wi_stream::wi_stream(int nfreq_, double freq_lo_MHz_, double freq_hi_MHz_, double dt_sample_, int nt_maxwrite_) :
    wi_stream()
{
    this->construct(nfreq_, freq_lo_MHz_, freq_hi_MHz_, dt_sample_, nt_maxwrite_);
}


void wi_stream::construct(int nfreq_, double freq_lo_MHz_, double freq_hi_MHz_, double dt_sample_, int nt_maxwrite_)
{
    this->nfreq = nfreq_;
    this->freq_lo_MHz = freq_lo_MHz_;
    this->freq_hi_MHz = freq_hi_MHz_;
    this->dt_sample = dt_sample_;
    this->nt_maxwrite = nt_maxwrite_;

    if (nfreq <= 0)
	throw runtime_error("wi_stream::construct(): expected nfreq > 0");
    if (freq_lo_MHz <= 0.0)
	throw runtime_error("wi_stream::construct(): expected freq_lo_MHz > 0.0");
    if (freq_lo_MHz >= freq_hi_MHz)
	throw runtime_error("wi_stream::construct(): expected freq_lo_MHz < freq_hi_MHz");
    if (dt_sample <= 0.0)
	throw runtime_error("wi_stream::construct(): expected dt_sample > 0.0");
    if (nt_maxwrite <= 0)
	throw runtime_error("wi_stream::construct(): expected nt_maxwrite > 0");
}


void wi_stream::run(const std::vector<std::shared_ptr<wi_transform> > &transforms)
{
    wi_run_state rstate(*this, transforms);

    // Delegate to pure virtual run() method implemented in subclass
    this->run(rstate);
}


}  // namespace rf_pipelines
