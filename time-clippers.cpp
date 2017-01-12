#include "rf_pipelines_internals.hpp"
#include "kernels/polyfit.hpp"

#include <mutex>
#include <cassert>
#include <condition_variable>

using namespace std;
using namespace rf_pipelines;

// Semi-placeholder version: just times the wi_transform objects (via transform_timing_thread base class).
// Will be expanded later to time kernels in a more granular way.

struct clipper_timing_thread : public transform_timing_thread
{
    const int Df;
    const int Dt;
    const int niter;

    clipper_timing_thread(const shared_ptr<timing_thread_pool> &pool_, int nfreq_, int nt_chunk_, int stride_, int Df_, int Dt_, int niter_) :
	transform_timing_thread{ pool_, nfreq_, nt_chunk_, stride_,
	    { make_intensity_clipper(nt_chunk_, AXIS_FREQ, 1.0e10, niter_, 1.0e10, Df_, Dt_),
	      make_intensity_clipper(nt_chunk_, AXIS_TIME, 1.0e10, niter_, 1.0e10, Df_, Dt_),
	      make_intensity_clipper(nt_chunk_, AXIS_NONE, 1.0e10, niter_, 1.0e10, Df_, Dt_),
	      make_std_dev_clipper(nt_chunk_, AXIS_FREQ, 1.0e10, Df_, Dt_),
	      make_std_dev_clipper(nt_chunk_, AXIS_TIME, 1.0e10, Df_, Dt_)
	    }
        },
	Df(Df_), 
	Dt(Dt_), 
	niter(niter_)
    { }

    virtual void thread_top() override
    {
	if (thread_id == 0) {
            cout << "time-clippers: nfreq=" << nfreq << ", nt_chunk=" << nt_chunk 
		 << ", stride=" << stride  << ", Df=" << Df << ", Dt=" << Dt 
		 << ", niter=" << niter << endl;
	}
    }

    virtual void thread_bottom() override { }
};


int main(int argc, char **argv)
{
    if (argc != 8) {
	cerr << "usage: time-clippers <nthreads> <nfreq> <nt_chunk> <stride> <Df> <Dt> <niter>\n";
	exit(2);
    }

    int nthreads = atoi(argv[1]);
    int nfreq = atoi(argv[2]);
    int nt_chunk = atoi(argv[3]);
    int stride = atoi(argv[4]);
    int Df = atoi(argv[5]);
    int Dt = atoi(argv[6]);
    int niter = atoi(argv[7]);

    assert(nfreq > 0);
    assert((nt_chunk > 0) && (nt_chunk % 8 == 0));
    assert(stride >= nt_chunk);
    assert(Df >= 1 && is_power_of_two(Df));
    assert(Dt >= 1 && is_power_of_two(Dt));
    assert(niter >= 1);
    
    cout << "nthreads = " << nthreads << endl;
    auto pool = make_shared<timing_thread_pool> (nthreads);

    vector<std::thread> threads(nthreads);
    for (int i = 0; i < nthreads; i++)
        threads[i] = spawn_timing_thread<clipper_timing_thread> (pool, nfreq, nt_chunk, stride, Df, Dt, niter);
    for (int i = 0; i < nthreads; i++)
        threads[i].join();

    return 0;
}
