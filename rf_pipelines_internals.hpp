#ifndef _RF_PIPELINES_INTERNALS_HPP
#define _RF_PIPELINES_INTERNALS_HPP

#if (__cplusplus < 201103) && !defined(__GXX_EXPERIMENTAL_CXX0X__)
#error "This source file needs to be compiled with C++0x support (g++ -std=c++0x)"
#endif

#include <cmath>
#include <cstring>
#include <iostream>
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

#define rf_assert_close(x,y,epsilon) rf_assert2(fabs((x)-(y)) <= epsilon, __LINE__)


namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif


inline double uniform_rand()
{
    return (rand() + 0.5) / (RAND_MAX + 1.0);
}

inline double uniform_rand(double lo, double hi)
{
    return lo + (hi-lo)*uniform_rand();
}

inline int randint(int lo, int hi)
{
    rf_assert(lo < hi);
    int ret = lo + (int)((hi-lo)*uniform_rand());
    ret = std::max(ret, lo);    // should be redundant
    ret = std::min(ret, hi-1);  // should be redundant
    return ret;
}

inline bool is_close(double x, double y)
{
    return fabs(x-y) <= 1.0e-5 * (fabs(x) + fabs(y));
}

// round up m to nearest multiple of n 
inline int round_up(int m, int n)
{
    rf_assert(m >= 0);
    rf_assert(n > 0);
    return ((m+n-1)/n) * n;
}

inline int gcd(int m, int n)
{
    if (m < n)
	std::swap(m, n);
    if (n < 0)
	throw std::runtime_error("gcd() called with negative argument");

    while (n > 0) {
	int d = m % n;
	m = n;
	n = d;
    }

    return m;
}

// std::vector doesn't provide a member function which guarantees deallocation!
template<typename T> static inline void deallocate(std::vector<T> &v)
{
    std::vector<T> w;
    v.swap(w);
}


}  // namespace rf_pipelines

#endif // _RF_PIPELINES_INTERNALS_HPP
