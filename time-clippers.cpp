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

    clipper_timing_thread(const shared_ptr<timing_thread_pool> &pool_, int nfreq_, int nt_chunk_, int stride_, int Df_, int Dt_, int niter_, bool twopass_) :
	transform_timing_thread{ pool_, nfreq_, nt_chunk_, stride_,
	    { make_intensity_clipper(nt_chunk_, AXIS_FREQ, 1.0e10, niter_, 1.0e10, Df_, Dt_, twopass_),
	      make_intensity_clipper(nt_chunk_, AXIS_TIME, 1.0e10, niter_, 1.0e10, Df_, Dt_, twopass_),
	      make_intensity_clipper(nt_chunk_, AXIS_NONE, 1.0e10, niter_, 1.0e10, Df_, Dt_, twopass_),
	      make_std_dev_clipper(nt_chunk_, AXIS_FREQ, 1.0e10, Df_, Dt_, twopass_),
	      make_std_dev_clipper(nt_chunk_, AXIS_TIME, 1.0e10, Df_, Dt_, twopass_)
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


static void usage()
{
    cerr << "usage: time-clippers [-t] <nthreads> <nfreq> <nt_chunk> <stride> <Df> <Dt> <niter>\n"
	 << "   The -t flag uses \"two-pass\" clippers\n";

    exit(2);
}


int main(int argc, char **argv)
{
    bool two_pass = false;

    std::vector<char *> args;

    for (int i = 1; i < argc; i++) {
	if (!strcmp(argv[i], "-t"))
	    two_pass = true;
	else
	    args.push_back(argv[i]);
    }

    if (args.size() != 7)
	usage();

    int nthreads = atoi(args[0]);
    int nfreq = atoi(args[1]);
    int nt_chunk = atoi(args[2]);
    int stride = atoi(args[3]);
    int Df = atoi(args[4]);
    int Dt = atoi(args[5]);
    int niter = atoi(args[6]);

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
        threads[i] = spawn_timing_thread<clipper_timing_thread> (pool, nfreq, nt_chunk, stride, Df, Dt, niter, two_pass);
    for (int i = 0; i < nthreads; i++)
        threads[i].join();

    return 0;
}
