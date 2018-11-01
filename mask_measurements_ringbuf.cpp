#include "rf_pipelines_internals.hpp"
#include "rf_pipelines_inventory.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif

typedef lock_guard<mutex> ulock;


mask_measurements::mask_measurements(ssize_t pos_, int nf_, int nt_)
{
    this->pos = pos_;
    this->nf = nf_;
    this->nt = nt_;
    this->nsamples = nf_ * nt_;
    this->freqs_unmasked = make_sptr<int> (nf_);
}


mask_measurements_ringbuf::mask_measurements_ringbuf(int nhistory) :
    next(0),
    maxsize(nhistory)
{
    if (nhistory <= 0)
	throw runtime_error("rf_pipelines::mask_measurements_ringbuf constructor called with nhistory <= 0");

    // Allocate ring buffer
    ringbuf.resize(nhistory);
}

void mask_measurements_ringbuf::add(rf_pipelines::mask_measurements& meas) {
    ulock l(mutex);
    ringbuf[next % maxsize] = meas;
    next++;
}
    
std::vector<rf_pipelines::mask_measurements>
mask_measurements_ringbuf::get_all_measurements() {
    std::vector<rf_pipelines::mask_measurements> copy;
    {
        ulock l(mutex);
        // The returned vector has the chunks listed in time order
        int start;
        int end;
        if (next <= maxsize) {
            start = 0;
            end = next;
        } else {
            start = next;
            end = next + maxsize;
        }
        for (int i=start; i<end; i++)
            copy.push_back(ringbuf[i % maxsize]);
    }
    return copy;
}

std::unordered_map<std::string, float> 
mask_measurements_ringbuf::get_stats(int nchunks) {
    unordered_map<string, float> stats;
    float totsamp = 0;
    float totunmasked = 0;
    if (nchunks > maxsize)
        nchunks = maxsize;
    {
        ulock l(mutex);
        int start = next - nchunks;
        if (start < 0)
            start = 0;
        for (int i=start; i<next; i++) {
            int reali = (i % maxsize);
            totsamp     += ringbuf[reali].nsamples;
            totunmasked += ringbuf[reali].nsamples_unmasked;
        }
    }
    stats["rfi_mask_pct_masked"]   = 100. * (totsamp - totunmasked) / max(totsamp, 1.f);
    return stats;
}

}  // namespace rf_pipelines
