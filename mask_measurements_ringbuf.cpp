#include "rf_pipelines_internals.hpp"
#include "rf_pipelines_inventory.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif

typedef lock_guard<mutex> ulock;

mask_measurements_ringbuf::mask_measurements_ringbuf(int nhistory) :
    current(0),
    maxsize(nhistory)
{}

void mask_measurements_ringbuf::add(rf_pipelines::mask_measurements& meas) {

    cout << "mask_measurements_ringbuf::add: got pos " << meas.pos 
	 << ": N samples masked: " << (meas.nsamples - meas.nsamples_unmasked)
	 << "/" << meas.nsamples << endl;

    ulock l(mutex);
    if (current < maxsize)
        ringbuf.push_back(meas);
    else
        // ring buffer
        ringbuf[current % maxsize] = meas;
    current++;
}

    
std::vector<rf_pipelines::mask_measurements>
mask_measurements_ringbuf::get_all_measurements() {
    std::vector<rf_pipelines::mask_measurements> copy;
    {
        ulock l(mutex);
        // Reorder the ring buffer.
        int n = ringbuf.size();
        for (int off=0; off<n; off++)
            copy.push_back(ringbuf[(current + 1 + off) % n]);
    }
    return copy;
}

std::unordered_map<std::string, float> 
mask_measurements_ringbuf::get_stats(float period) {
    unordered_map<string, float> stats;
    // FIXME -- assume one sample per second!
    int nsteps = (int)period;
    float totsamp = 0;
    float totunmasked = 0;
    {
        ulock l(mutex);
        int n = ringbuf.size();
        //cout << "Ringbuf n: " << n << ", current " << current << endl;
        if (nsteps > n)
            nsteps = n;
        // avoid  (i % 0)
        int istart = (n ? (current - nsteps + n) % n : 0);
        for (int offset=0; offset<nsteps; offset++) {
            int i = (istart + offset) % n;
            //cout << "offset " << offset << " of " << nsteps << " -> i " << i << endl;
            totsamp += ringbuf[i].nsamples;
            totunmasked += ringbuf[i].nsamples_unmasked;
        }
    }
    stats["rfi_mask_pct_masked"]   = 100. * (totsamp - totunmasked) / max(totsamp, 1.f);
    return stats;
}



}  // namespace rf_pipelines
