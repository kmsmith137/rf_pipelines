#include <random>
#include "rf_pipelines_internals.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif


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


// -------------------------------------------------------------------------------------------------
//
// timing_thread_pool


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


// -------------------------------------------------------------------------------------------------
//
// timing_thread


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
// transform_timing_thread (+ helper class fake_stream)


struct fake_stream : wi_stream {
    fake_stream(ssize_t nfreq_, ssize_t nt_maxwrite_)
    {
	// Initialize base class members
	this->nfreq = nfreq_;
	this->freq_lo_MHz = 400.0;
	this->freq_hi_MHz = 800.0;
	this->dt_sample = 1.0e-3;
	this->nt_maxwrite = nt_maxwrite_;
    }

    virtual void stream_body(wi_run_state &run_state) override
    {
	throw runtime_error("fake_stream::stream_body() should never be called");
    }
};


transform_timing_thread::transform_timing_thread(const shared_ptr<timing_thread_pool> &pool_, int nfreq_, int nt_chunk_, int stride_, 
						 const vector<shared_ptr<wi_transform>> &transform_list_) :
    timing_thread(pool_, true),    // pin_to_core=true
    nfreq(nfreq_),
    nt_chunk(nt_chunk_),
    stride(stride_),
    transform_list(transform_list_),
    ntransforms(transform_list_.size())
{ 
    rf_assert(nfreq > 0);
    rf_assert(nt_chunk > 0 && (nt_chunk % 8 == 0));
    rf_assert(stride >= nt_chunk);
    rf_assert(ntransforms > 0);
    
    intensity = aligned_alloc<float> (nfreq * stride);
    weights = aligned_alloc<float> (nfreq * stride);
}


transform_timing_thread::~transform_timing_thread()
{
    free(intensity);
    free(weights);
    intensity = weights = nullptr;
}


void transform_timing_thread::thread_body()
{
    std::random_device rd;
    std::mt19937 rng(rd());
    
    if (thread_id == 0) {
	cout << "nfreq=" << nfreq << ", nt_chunk=" << nt_chunk 
	     << ", stride=" << stride << ", niter=" << niter << endl;
    }

    fake_stream s(nfreq, nt_chunk);

    for (int itr = 0; itr < ntransforms; itr++) {
	
	// This is a little awkward: we don't want the buffer-resets (see below) 
	// to contribute to the reported time, so we start and stop the timer in
	// every iteration of the loop.
	
	double total_time = 0.0;
	
	transform_list[itr]->set_stream(s);
	transform_list[itr]->start_substream(0, 0.0);
	
	// A few sanity checks
	if (transform_list[itr]->nfreq != nfreq)
	    throw runtime_error("rf_pipelines::transform_timing_thread: nfreq mismatch");
	if (transform_list[itr]->nt_chunk != nt_chunk)
	    throw runtime_error("rf_pipelines::transform_timing_thread: nt_chunk mismatch");
	if ((transform_list[itr]->nt_prepad != 0) || (transform_list[itr]->nt_postpad != 0))
	    throw runtime_error("rf_pipelines::transform_timing_thread expects nt_prepad=nt_prepad=0 (this would be easy to fix if needed)");
	
	for (int iter = 0; iter < niter; iter++) {

	    // We reset buffers between calls to wi_transform::process_chunk().
	    // This is to avoid unexpected situations, such as a clipper which eventually
	    // clips the entire buffer if it is iterated many times.

	    for (int i = 0; i < nfreq * stride; i++) {
		intensity[i] = uniform_real_distribution<>()(rng);
		weights[i] = 1.0;
	    }

	    double tchunk = s.dt_sample * nt_chunk;

	    this->start_timer();
	    transform_list[itr]->process_chunk(iter*tchunk, (iter+1)*tchunk, intensity, weights, stride, nullptr, nullptr, 0);
	    total_time += this->stop_timer();
	}

	transform_list[itr]->end_substream();
    }

    this->timing_thread_body();
}


}   // namespace rf_pipelines
