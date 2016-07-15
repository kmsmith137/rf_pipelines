#ifndef _RF_PIPELINES_INTERNALS_HPP
#define _RF_PIPELINES_INTERNALS_HPP

#if (__cplusplus < 201103) && !defined(__GXX_EXPERIMENTAL_CXX0X__)
#error "This source file needs to be compiled with C++0x support (g++ -std=c++0x)"
#endif

#include <cmath>
#include <cstring>
#include <stdexcept>

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


// Non-inline helper functions (more to come?)
extern bool file_exists(const std::string &filename);
extern void listdir(std::vector<std::string> &filenames, const std::string &dirname);


// Inlines follow...

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


}  // namespace rf_pipelines

#endif // _RF_PIPELINES_INTERNALS_HPP
