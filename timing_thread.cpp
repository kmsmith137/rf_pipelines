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


}   // namespace rf_pipelines
