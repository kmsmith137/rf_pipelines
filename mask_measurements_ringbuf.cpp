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

    cout << "mask_measurements_ringbuf::add: got pos " << meas.pos << ": N samples masked: " << meas.nsamples_masked << "/" << meas.nsamples << "; n times " << meas.nt_masked << "/" << meas.nt << "; n freqs " << meas.nf_masked << "/" << meas.nf << endl;

    ulock l(mutex);
    if (current < maxsize)
        ringbuf.push_back(meas);
    else
        // ring buffer
        ringbuf[current] = meas;
    current = (current + 1) % maxsize;
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
    float totmasked = 0;
    float tot_t = 0;
    float tot_tmasked = 0;
    float tot_fmasked = 0;
    float tot_f = 0;
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
            totmasked += ringbuf[i].nsamples_masked;
            tot_t += ringbuf[i].nt;
            tot_tmasked += ringbuf[i].nt_masked;
            tot_f += ringbuf[i].nf;
            tot_fmasked += ringbuf[i].nf_masked;
        }
    }
    stats["rfi_mask_pct_masked"]   = 100. * totmasked / max(totsamp, 1.f);
    stats["rfi_mask_pct_t_masked"] = 100. * tot_tmasked / max(tot_t, 1.f);
    stats["rfi_mask_pct_f_masked"] = 100. * tot_fmasked / max(tot_f, 1.f);
    return stats;
}



}  // namespace rf_pipelines
