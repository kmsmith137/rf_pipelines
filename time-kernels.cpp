#include "rf_pipelines_internals.hpp"
#include "kernels/polyfit.hpp"

#include <mutex>
#include <thread>
#include <cassert>
#include <condition_variable>

using namespace std;
using namespace rf_pipelines;


// -------------------------------------------------------------------------------------------------
//
// General-purpose timing thread, part 1.  
// This can go in rf_pipelines_internals.hpp if it needs to be made available elsewhere.


class timing_thread_pool {
public:
    const int nthreads;

    timing_thread_pool(int nthreads);

    typedef struct timeval time_t;

    time_t start_timer();
    double stop_timer(const time_t &start_time);
    
    // Helper function called by timing_thread.
    int get_and_increment_thread_id();

protected:
    std::mutex lock;
    std::condition_variable cond0;
    std::condition_variable cond1;
    std::condition_variable cond2;

    double total_dt = 0.0;
    int threads_so_far = 0;

    int ix0 = 0;
    int ix1 = 0;
    int ix2 = 0;
};


class timing_thread {
public:
    const std::shared_ptr<timing_thread_pool> pool;
    const bool pinned_to_core;
    const int thread_id;
    const int nthreads;

    static void _thread_main(timing_thread *t);

    virtual ~timing_thread() { }

protected:
    timing_thread(const std::shared_ptr<timing_thread_pool> &pool, bool pin_to_core);

    virtual void thread_body() = 0;

    timing_thread_pool::time_t start_time;
    bool timer_is_running = false;

    void start_timer();

    // If 'name' is non-null, then timing will be announced on thread ID zero.
    double stop_timer(const char *name=nullptr);
};


template<typename T, typename... Args>
std::thread spawn_timing_thread(Args... args)
{
    timing_thread *t = new T(args...);
    return std::thread(timing_thread::_thread_main, t);
}


// -------------------------------------------------------------------------------------------------
//
// General-purpose timing thread, part 2.
// This can go in its own source file if it needs to be made available elsewhere.



static void pin_current_thread_to_core(int core_id)
{
#ifdef __APPLE__
    if (core_id == 0)
	cerr << "warning: pinning threads to cores is not implemented in osx\n";
    return;
#else
    int hwcores = std::thread::hardware_concurrency();
    
    if ((core_id < 0) || (core_id >= hwcores))
	throw runtime_error("pin_thread_to_core: core_id=" + to_string(core_id) + " is out of range (hwcores=" + to_string(hwcores) + ")");

    pthread_t thread = pthread_self();

    cpu_set_t cs;
    CPU_ZERO(&cs);
    CPU_SET(core_id, &cs);

    int err = pthread_setaffinity_np(thread, sizeof(cs), &cs);
    if (err)
        throw runtime_error("pthread_setaffinity_np() failed");
#endif
}


timing_thread_pool::timing_thread_pool(int nthreads_) :
    nthreads(nthreads_)
{ 
    if (nthreads <= 0)
	throw runtime_error("timing_thread_pool constructor called with nthreads <= 0");
}


timing_thread_pool::time_t timing_thread_pool::start_timer()
{
    unique_lock<mutex> l(lock);

    ix0++;
    if (ix0 == nthreads) {
	ix1 = ix2 = 0;
	cond0.notify_all();
    }
    while (ix0 < nthreads)
	cond0.wait(l);

    time_t start_time;
    gettimeofday(&start_time, NULL);
    return start_time;
}


double timing_thread_pool::stop_timer(const time_t &start_time)
{
    time_t end_time;
    gettimeofday(&end_time, NULL);

    double local_dt = (end_time.tv_sec - start_time.tv_sec) + 1.0e-6 * (end_time.tv_usec - start_time.tv_usec);

    unique_lock<mutex> l(this->lock);
    total_dt += local_dt;

    ix1++;
    if (ix1 == nthreads) {
	ix0 = ix2 = 0;
	cond1.notify_all();
    }
    while (ix1 < nthreads)
	cond1.wait(l);

    double ret = total_dt / nthreads;

    ix2++;
    if (ix2 == nthreads) {
	total_dt = 0.0;
	ix0 = ix1 = 0;
	cond2.notify_all();
    }
    while (ix2 < nthreads)
	cond2.wait(l);

    return ret;
}


int timing_thread_pool::get_and_increment_thread_id()
{
    lock_guard<mutex> l(lock);
    return threads_so_far++;
}


timing_thread::timing_thread(const shared_ptr<timing_thread_pool> &pool_, bool pin_to_core) :
    pool(pool_), 
    pinned_to_core(pin_to_core),
    thread_id(pool_->get_and_increment_thread_id()),
    nthreads(pool_->nthreads)
{ }


// static member function
void timing_thread::_thread_main(timing_thread *t)
{
    // Ensure delete(t) is called
    auto p = unique_ptr<timing_thread> (t);

    if (t->pinned_to_core)
	pin_current_thread_to_core(t->thread_id);

    t->thread_body();
}


void timing_thread::start_timer()
{
    if (timer_is_running)
	throw runtime_error("double call to timing_thread::start_timer(), without call to stop_timer() in between");

    this->timer_is_running = true;
    this->start_time = pool->start_timer();
}


double timing_thread::stop_timer(const char *name)
{
    if (!timer_is_running)
	throw runtime_error("timing_thread::stop_timer() called without calling start_timer(), or double call to stop_timer()");

    this->timer_is_running = false;

    double ret = pool->stop_timer(start_time);

    if (name && !thread_id)
	cout << (string(name) + ": " + to_string(ret) + " seconds\n");

    return ret;
}


// -------------------------------------------------------------------------------------------------
//
// Kernel timing code starts here


// The template parameter N is (polydeg + 1)
template<typename T, unsigned int S, unsigned int N>
struct kernel_timing_thread : public timing_thread
{
    const int nfreq;
    const int nt_chunk;
    const int stride;
    const int niter = 16;

    float *intensity = nullptr;
    float *weights = nullptr;

    // A place to write dummy results, to keep the compiler from optimizing things out
    float *dummyp = nullptr;

    kernel_timing_thread(const shared_ptr<timing_thread_pool> &pool_, int nfreq_, int nt_chunk_, int stride_) :
	timing_thread(pool_, true),    // pin_to_core=true
	nfreq(nfreq_), nt_chunk(nt_chunk_), stride(stride_)
    { 
	assert(nfreq > 0);
	assert(nt_chunk > 0 && (nt_chunk % 8 == 0));
	assert(stride >= nt_chunk);

	intensity = aligned_alloc<float> (nfreq * stride);
	weights = aligned_alloc<float> (nfreq * stride);
	dummyp = aligned_alloc<float> (16);

	for (int i = 0; i < nfreq*stride; i++)
	    weights[i] = 1.0;
    }

    virtual void thread_body() override
    {
	simd_trimatrix<T,S,N> smat;
	simd_ntuple<T,S,N> svec;
	simd_t<T,S> dummy(0.0);

	if (thread_id == 0) {
            cout << "nfreq=" << nfreq << ", nt_chunk=" << nt_chunk 
		 << ", stride=" << stride  << ", polydeg=" << (N-1)
		 << ", niter=" << niter << endl;
	}

	this->start_timer();
        for (int iter = 0; iter < niter; iter++)
	    _kernel_detrend_t<T,S,N> (nfreq, nt_chunk, intensity, weights, stride);
        this->stop_timer("kernel_detrend_t");

#if 0
	this->start_timer();
        for (int iter = 0; iter < niter; iter++) {
	    for (int ifreq = 0; ifreq < nfreq; ifreq++) {
		_kernel_detrend_t_pass1<T,S,N> (smat, svec, nt_chunk, intensity + ifreq*stride, weights + ifreq*stride);
		dummy += smat.vertical_sum();
		dummy += svec.vertical_sum();
	    }
	}
        this->stop_timer("kernel_detrend_t_pass1");
#endif

	this->start_timer();
        for (int iter = 0; iter < niter; iter++)
	    _kernel_detrend_f<T,S,N> (nfreq, nt_chunk, intensity, weights, stride);
        this->stop_timer("kernel_detrend_f");

#if 0
	this->start_timer();
        for (int iter = 0; iter < niter; iter++) {
	    for (int it = 0; it < nt_chunk; it += S) {
		_kernel_detrend_f_pass1<T,S,N> (smat, svec, nfreq, intensity + it, weights + it, stride);
		dummy += smat.vertical_sum();
		dummy += svec.vertical_sum();
	    }
	}
        this->stop_timer("kernel_detrend_f_pass1");
#endif

	dummy.store(dummyp);
    }
};


template<typename T, unsigned int S, unsigned int Nmax, typename std::enable_if<(Nmax>0),int>::type = 0>
inline std::thread make_kernel_timing_thread(const shared_ptr<timing_thread_pool> &pool, int polydeg, int nfreq, int nt_chunk, int stride)
{
    if (Nmax == polydeg + 1)
	return spawn_timing_thread< kernel_timing_thread<T,S,Nmax> > (pool, nfreq, nt_chunk, stride);
    return make_kernel_timing_thread<T,S,(Nmax-1)> (pool, polydeg, nfreq, nt_chunk, stride);
}

template<typename T, unsigned int S, unsigned int Nmax, typename std::enable_if<(Nmax==0),int>::type = 0>
inline std::thread make_kernel_timing_thread(const shared_ptr<timing_thread_pool> &pool, int polydeg, int nfreq, int nt_chunk, int stride)
{
    throw runtime_error("internal error in make_kernel_timing_thread()");
}


int main(int argc, char **argv)
{
    // (max polynomial degree) + 1
    static const int Nmax = 17;

    if (argc != 6) {
	cerr << "usage: time-kernels <nfreq> <nt_chunk> <stride> <polydeg> <nthreads>\n";
	exit(2);
    }

    int nfreq = atoi(argv[1]);
    int nt_chunk = atoi(argv[2]);
    int stride = atoi(argv[3]);
    int polydeg = atoi(argv[4]);
    int nthreads = atoi(argv[5]);

    assert(nfreq > 0);
    assert((nt_chunk > 0) && (nt_chunk % 8 == 0));
    assert(stride >= nt_chunk);
    assert(polydeg >= 0 && polydeg < Nmax);
    assert(nthreads > 0 && nthreads <= 20);
    
    cout << "nthreads = " << nthreads << endl;
    auto pool = make_shared<timing_thread_pool> (nthreads);

    vector<std::thread> threads(nthreads);
    for (int i = 0; i < nthreads; i++)
        threads[i] = make_kernel_timing_thread<float,8,Nmax> (pool, polydeg, nfreq, nt_chunk, stride);
    for (int i = 0; i < nthreads; i++)
        threads[i].join();

    return 0;
}
