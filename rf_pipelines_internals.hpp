#ifndef _RF_PIPELINES_INTERNALS_HPP
#define _RF_PIPELINES_INTERNALS_HPP

#if (__cplusplus < 201103) && !defined(__GXX_EXPERIMENTAL_CXX0X__)
#error "This source file needs to be compiled with C++11 support (g++ -std=c++11)"
#endif

#include <random>
#include <cstring>
#include <sstream>
#include <algorithm>
#include <unordered_set>
#include <sys/time.h>

#include "rf_kernels/core.hpp"
#include "rf_pipelines_base_classes.hpp"
#include "rf_pipelines_inventory.hpp"

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
}  // emacs pacifier
#endif


// -------------------------------------------------------------------------------------------------
//
// outdir_manager, plot_group.


struct outdir_manager {
    // If 'outdir' is an empty string, then the pipeline has no output directory,
    // and calls to outdir_manager::add_file() or pipeline_object::add_plot() will fail.

    std::string outdir;  // if nonempty string, includes trailing slash.
    bool clobber_ok = true;

    std::unordered_set<std::string> basenames;

    // Constructor creates the output directory.
    outdir_manager(const std::string &outdir, bool clobber_ok);

    // Called before file is written, to test for filename collisions.
    // Returns the full pathname, throws exception if filename has already been written in this pipeline run.
    std::string add_file(const std::string &basename);
};


struct plot_group {
    std::string name;
    int nt_per_pix = 0;
    int ny = 0;
    
    bool is_empty = true;
    int64_t curr_it0 = 0;
    int64_t curr_it1 = 0;
    Json::Value files;
};


// -------------------------------------------------------------------------------------------------


// file_utils.cpp
extern bool file_exists(const std::string &filename);
extern std::vector<std::string> listdir(const std::string &dirname);
extern void makedirs(const std::string &dirname);


// json_utils.cpp
extern Json::Value array_from_json(const Json::Value &x, const std::string &k);
extern std::string string_from_json(const Json::Value &x, const std::string &k);
extern rf_kernels::axis_type axis_type_from_json(const Json::Value &x, const std::string &k);
extern double double_from_json(const Json::Value &x, const std::string &k);
extern int int_from_json(const Json::Value &x, const std::string &k);
extern bool bool_from_json(const Json::Value &j, const std::string &k);
extern void add_json_object(Json::Value &dst, const Json::Value &src);


// -------------------------------------------------------------------------------------------------
//
// lexical_cast


// Utility routine: converts a string to type T (only a few T's are defined; see lexical_cast.cpp)
// Returns true on success, false on failure
template<typename T> extern bool lexical_cast(const std::string &x, T &ret);

// Also defined in lexical_cast.cpp (for the same values of T)
template<typename T> extern const char *typestr();

// Version of lexical_cast() which throws exception on failure.
template<typename T> inline T lexical_cast(const std::string &x, const char *name="string")
{
    T ret;
    if (lexical_cast(x, ret))
	return ret;
    throw std::runtime_error("couldn't convert " + std::string(name) + "='" + x + "' to " + typestr<T>());
}


// -------------------------------------------------------------------------------------------------
//
// RNG helpers


inline double uniform_rand(std::mt19937 &rng)
{
    return std::uniform_real_distribution<>()(rng);
}

inline double uniform_rand(std::mt19937 &rng, double lo, double hi)
{
    return lo + (hi-lo) * uniform_rand(rng);
}

inline ssize_t randint(std::mt19937 &rng, ssize_t lo, ssize_t hi)
{
    if (lo >= hi)
	throw std::runtime_error("rf_pipelines internal error: expected (lo < hi) in randint()");

    return std::uniform_int_distribution<>(lo,hi-1)(rng);   // note hi-1 here!
}

inline double dist(double x, double y)
{
    return fabs(x-y);
}

inline double reldist(double x, double y)
{
    double num = fabs(x-y);
    double den = fabs(x) + fabs(y);
    return (den > 0.0) ? (num/den) : 0.0;
}

inline std::vector<float> uniform_randvec(std::mt19937 &rng, ssize_t n, double lo, double hi)
{
    rf_assert(n > 0);
    
    std::vector<float> ret(n);
    for (ssize_t i = 0; i < n; i++)
	ret[i] = uniform_rand(rng, lo, hi);

    return ret;
}


// -------------------------------------------------------------------------------------------------
//
// Misc inlines


template<typename T>
inline T *aligned_alloc(size_t nelts, size_t nalign=128, bool zero=true)
{
    if (nelts == 0)
        return NULL;

    void *p = NULL;
    if (posix_memalign(&p, nalign, nelts * sizeof(T)) != 0)
        throw std::runtime_error("couldn't allocate memory");

    if (zero)
	memset(p, 0, nelts * sizeof(T));

    return reinterpret_cast<T *> (p);
}

template<typename T, typename... Args>
std::unique_ptr<T> make_unique(Args&& ...args)
{
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
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

// Returns (m/n), in a situation where we want to assert that n evenly divides m.
inline ssize_t xdiv(ssize_t m, ssize_t n)
{
    rf_assert(m >= 0);
    rf_assert(n > 0);
    rf_assert(m % n == 0);
    return m / n;
}

// Returns (m % n), in a situation where we want to assert that the % operation makes sense
inline ssize_t xmod(ssize_t m, ssize_t n)
{
    rf_assert(m >= 0);
    rf_assert(n > 0);
    return m % n;
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

// round up m to nearest multiple of n
inline ssize_t round_up(ssize_t m, ssize_t n)
{
    rf_assert(m >= 0);
    rf_assert(n > 0);
    return ((m+n-1)/n) * n;
}

template<typename T> inline T prod(const std::vector<T> &v)
{
    T ret = 1;
    for (unsigned int i = 0; i < v.size(); i++)
        ret *= v[i];
    return ret;
}

template<typename T, typename U, typename V>
inline bool has_key(const std::unordered_map<T,U> &m, const V &key)
{
    return (m.find(T(key)) != m.end());
}


template<typename T>
inline void stringify(std::stringstream &ss, const T &x)
{
    ss << x;
}

template<typename T>
inline void stringify(std::stringstream &ss, const std::vector<T> &v)
{
    ss << "[";
    for (unsigned int i = 0; i < v.size(); i++) {
	if (i > 0)
	    ss << ",";
	ss << v[i];
    }
    ss << "]";
}

template<typename T>
inline std::string stringify(const T &x)
{
    std::stringstream ss;
    stringify(ss, x);
    return ss.str();
}

inline bool startswith(const std::string &str, const std::string &prefix)
{
    return std::equal(prefix.begin(), prefix.end(), str.begin());
}

inline bool endswith(const std::string &str, const std::string &suffix)
{
    return std::equal(suffix.rbegin(), suffix.rend(), str.rbegin());
}


template<typename T>
inline T median(std::vector<T> &v)
{
    rf_assert(v.size() > 0);

    int n = v.size();
    int m = n/2;

    std::nth_element(v.begin(), v.begin()+m, v.end());
    
    if (n == 2*m+1)
	return v[m];

    return 0.5 * (v[m] + *std::max_element(v.begin(), v.begin()+m));
}


}  // namespace rf_pipelines

#endif  // _RF_PIPELINES_INTERNALS_HPP
