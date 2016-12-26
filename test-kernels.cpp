// Note: this file should only include kernel headers (kernels/*.hpp), not toplevel headers (./*.hpp)
// (This is because it has a Makefile dependency on the former but not the latter.)

#include "kernels/clip2d.hpp"
#include "simd_helpers/simd_debug.hpp"

using namespace std;
using namespace rf_pipelines;


// -------------------------------------------------------------------------------------------------
//
// A few things cut-and-pasted from rf_pipelines_internals.hpp, since we don't want to depend on it.


#ifndef _unlikely
#define _unlikely(cond)  (__builtin_expect(cond,0))
#endif

#define rf_assert(cond) rf_assert2(cond, __LINE__)

#define rf_assert2(cond,line) \
    do { \
        if (_unlikely(!(cond))) { \
	    const char *msg = "rf_pipelines: assertion '" __STRING(cond) "' failed (" __FILE__ ":" __STRING(line) ")\n"; \
	    throw std::runtime_error(msg); \
	} \
    } while (0)


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


// -------------------------------------------------------------------------------------------------


struct random_chunk {
    const int nfreq;
    const int nt;
    const int stride;
    
    float *intensity = nullptr;
    float *weights = nullptr;

    random_chunk(std::mt19937 &rng, int nfreq, int nt);
    ~random_chunk();
    
    // noncopyable
    random_chunk(const random_chunk &) = delete;
    random_chunk &operator=(const random_chunk &) = delete;
};


random_chunk::random_chunk(std::mt19937 &rng, int nfreq_, int nt_)
    : nfreq(nfreq_), nt(nt_), stride(nt + std::uniform_int_distribution<>(0,4)(rng))
{
    rf_assert(nfreq > 0);
    rf_assert(nt > 0);

    intensity = aligned_alloc<float> (nfreq * stride);
    weights = aligned_alloc<float> (nfreq * stride);

    std::normal_distribution<> gdist;

    for (int i = 0; i < nfreq * stride; i++) {
	intensity[i] = gdist(rng) + 1.0;
	weights[i] = std::uniform_real_distribution<>()(rng);
    }
}


random_chunk::~random_chunk()
{
    free(intensity);
    free(weights);
    intensity = weights = nullptr;
}


// -------------------------------------------------------------------------------------------------


template<typename T>
static void reference_clip2d_wrms(T &mean, T &rms, const T *intensity, const T *weights, int nfreq, int nt, int stride, int nds_f, int nds_t)
{
    rf_assert(nfreq % nds_f == 0);
    rf_assert(nt % nds_t == 0);

    // double-precision here
    double acc0 = 0.0;
    double acc1 = 0.0;
    double acc2 = 0.0;
    
    for (int ifreq = 0; ifreq < nfreq; ifreq += nds_f) {
	for (int it = 0; it < nt; it += nds_t) {
	    T ival = 0.0;
	    T wval = 0.0;

	    for (int jfreq = ifreq; jfreq < ifreq + nds_f; jfreq++) {
		for (int jt = it; jt < it + nds_t; jt++) {
		    ival += intensity[jfreq*stride + jt];
		    wval += weights[jfreq*stride + jt];
		}
	    }

	    acc0 += double(wval);
	    acc1 += double(wval) * double(ival);
	    acc2 += double(wval) * double(ival) * double(ival);
	}
    }

    // FIXME case of invalid entries not tested
    mean = acc1/acc0;
    rms = sqrt(acc2/acc0 - mean*mean);
}


template<typename T, unsigned int S, unsigned int Df, unsigned int Dt>
static void test_kernel_clip2d_wrms(std::mt19937 &rng)
{
    int nfreq = Df * std::uniform_int_distribution<>(10,20)(rng);
    int nt = Dt * S * std::uniform_int_distribution<>(10,20)(rng);

    random_chunk rc(rng, nfreq, nt);
    
    simd_t<T,S> mean, rms;
    _kernel_clip2d_wrms<T,S,Df,Dt> (mean, rms, rc.intensity, rc.weights, nfreq, nt, rc.stride);

    float ref_mean, ref_rms;
    reference_clip2d_wrms(ref_mean, ref_rms, rc.intensity, rc.weights, nfreq, nt, rc.stride, Df, Dt);

    cout << ref_mean << " | " << mean << endl;
    cout << ref_rms << " | " << rms << endl;
}


int main(int argc, char **argv)
{
    std::random_device rd;
    std::mt19937 rng(rd());
    
    test_kernel_clip2d_wrms<float,8,1,1> (rng);
    return 0;
}
