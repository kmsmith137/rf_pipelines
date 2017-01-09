#ifndef _RF_PIPELINES_INTERNALS_HPP
#define _RF_PIPELINES_INTERNALS_HPP

#if (__cplusplus < 201103) && !defined(__GXX_EXPERIMENTAL_CXX0X__)
#error "This source file needs to be compiled with C++0x support (g++ -std=c++0x)"
#endif

#include <cmath>
#include <cstring>
#include <sstream>
#include <stdexcept>
#include <sys/time.h>

#include <thread>
#include <condition_variable>

#include "rf_pipelines.hpp"


// Branch predictor hint
#ifndef _unlikely
#define _unlikely(cond)  (__builtin_expect(cond,0))
#endif

// rf_assert(): like assert, but throws an exception in order to work smoothly with python.
#define rf_assert(cond) rf_assert2(cond, __LINE__)

#define rf_assert2(cond,line) \
    do { \
        if (_unlikely(!(cond))) { \
	    const char *msg = "rf_pipelines: assertion '" __STRING(cond) "' failed (" __FILE__ ":" __STRING(line) ")\n"; \
	    throw std::runtime_error(msg); \
	} \
    } while (0)


namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif


//
// A note for the future: if rf_pipelines is ever made multithreaded, then there are
// some race conditions related to the output_tracker which will need to be fixed.  The
// basename_set should be protected by a lock, and we also probably want a lock to
// protect the wi_transform::output_output_tracker pointer itself.  We may also want the
// output_tracker to create a lockfile in the output directory.
//
struct outdir_manager {
    std::string outdir;  // can be an empty string, otherwise includes trailing slash
    bool clobber_ok = true;

    std::set<std::string> basename_set;

    // Constructor creates the output directory.
    outdir_manager(const std::string &outdir, bool clobber_ok);

    // Returns the full pathname, throws exception if filename has already been written in this pipeline run.
    std::string add_file(const std::string &basename);

    void write_per_substream_json_file(int isubstream, const Json::Value &data, bool noisy);

    static bool is_json_basename(const std::string &basename);
};


struct plot_group {
    std::string name;
    int nt_per_pix = 0;
    int ny = 0;
    
    bool is_empty = true;
    int64_t curr_it0 = 0;
    int64_t curr_it1 = 0;
    Json::Value files;

    plot_group(const std::string &name_, int nt_per_pix_, int ny_) :
	name(name_), nt_per_pix(nt_per_pix_), ny(ny_) 
    { 
	if (nt_per_pix < 1)
	    throw std::runtime_error("rf_pipelines::plot_group: nt_per_pix must be >= 1");
	if (ny < 1)
	    throw std::runtime_error("rf_pieplines::plot_group: ny must be >= 1");
    }
};


// -------------------------------------------------------------------------------------------------
//
// timing_thread (general-purpose timing thread), and transform_timing_thread (subclass for timing wi_transforms).


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

    // Thread-collective: all threads wait at a barrier, then initialize their local timers.
    void start_timer();

    // Thread-collective: the returned time is the average taken over all threads.
    // If 'name' is non-null, then timing will be announced on thread ID zero.
    double stop_timer(const char *name=nullptr);
};


struct transform_timing_thread : public timing_thread
{
    const int nfreq;
    const int nt_chunk;
    const int stride;
    const int niter = 16;

    float *intensity = nullptr;
    float *weights = nullptr;

    std::vector<std::shared_ptr<wi_transform>> transform_list;
    int ntransforms = 0;

    transform_timing_thread(const std::shared_ptr<timing_thread_pool> &pool, int nfreq, int nt_chunk, int stride, 
			    const std::vector<std::shared_ptr<wi_transform>> &transform_list);

    ~transform_timing_thread();

    // Noncopyable
    transform_timing_thread(const transform_timing_thread &) = delete;
    transform_timing_thread &operator=(const transform_timing_thread &) = delete;

    virtual void thread_top();             // default thread_top(): prints nfreq, nt_chunk, stride on thread_id 0
    virtual void thread_body() override;   // overrides timing_thread::thread_body(), times transforms in transform_list
    virtual void thread_bottom() {}        // optional: if the timing thread should do anything else, it can go here
};


// Can be called for any subclass T of timing_thread (including T=transform_timing_thread)
template<typename T, typename... Args>
std::thread spawn_timing_thread(Args... args)
{
    timing_thread *t = new T(args...);
    return std::thread(timing_thread::_thread_main, t);
}


// -------------------------------------------------------------------------------------------------


// Non-inline helper functions (more to come?)
extern bool file_exists(const std::string &filename);
extern void makedirs(const std::string &dirname);
extern std::vector<std::string> listdir(const std::string &dirname);


// Inlines follow...

inline bool is_power_of_two(int n)
{
    rf_assert(n >= 1);
    return (n & (n-1)) == 0;
}

inline double uniform_rand()
{
    return (rand() + 0.5) / (RAND_MAX + 1.0);
}

inline double uniform_rand(double lo, double hi)
{
    return lo + (hi-lo)*uniform_rand();
}

inline ssize_t randint(ssize_t lo, ssize_t hi)
{
    rf_assert(lo < hi);
    ssize_t ret = lo + (ssize_t)((hi-lo)*uniform_rand());
    ret = std::max(ret, lo);    // should be redundant
    ret = std::min(ret, hi-1);  // should be redundant
    return ret;
}

inline double dist(double x, double y)
{
    return fabs(x-y);
}

inline double reldist(double x, double y)
{
    return fabs(x-y) / (fabs(x) + fabs(y));
}

// round up m to nearest multiple of n 
inline ssize_t round_up(ssize_t m, ssize_t n)
{
    rf_assert(m >= 0);
    rf_assert(n > 0);
    return ((m+n-1)/n) * n;
}

inline ssize_t gcd(ssize_t m, ssize_t n)
{
    if (m < n)
	std::swap(m, n);
    if (n < 0)
	throw std::runtime_error("gcd() called with negative argument");

    while (n > 0) {
	ssize_t d = m % n;
	m = n;
	n = d;
    }

    return m;
}

inline bool startswith(const std::string &str, const std::string &prefix)
{
    return std::equal(prefix.begin(), prefix.end(), str.begin());
}

inline bool endswith(const std::string &str, const std::string &suffix)
{
    return std::equal(suffix.rbegin(), suffix.rend(), str.rbegin());
}

template<typename T, typename... Args>
std::unique_ptr<T> make_unique(Args&& ...args)
{
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

template<typename T>
inline T *aligned_alloc(size_t nelts)
{
    if (nelts == 0)
	return NULL;

    // align to 64-byte cache lines
    void *p = NULL;
    if (posix_memalign(&p, 64, nelts * sizeof(T)) != 0)
	throw std::runtime_error("couldn't allocate memory");

    memset(p, 0, nelts * sizeof(T));
    return reinterpret_cast<T *> (p);
}

// std::vector doesn't provide a member function which guarantees deallocation!
template<typename T> static inline void deallocate(std::vector<T> &v)
{
    std::vector<T> w;
    v.swap(w);
}

inline double time_diff(const struct timeval &tv1, const struct timeval &tv2)
{
    return (tv2.tv_sec - tv1.tv_sec) + 1.0e-6 * (tv2.tv_usec - tv1.tv_usec);
}

inline struct timeval get_time()
{
    struct timeval ret;
    if (gettimeofday(&ret, NULL) < 0)
	throw std::runtime_error("gettimeofday() failed");
    return ret;
}


}  // namespace rf_pipelines

#endif // _RF_PIPELINES_INTERNALS_HPP
