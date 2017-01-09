#include "rf_pipelines_internals.hpp"
#include "kernels/polyfit.hpp"

#include <mutex>
#include <cassert>
#include <condition_variable>

using namespace std;
using namespace rf_pipelines;


// The template parameter N is (polydeg + 1)
template<typename T, unsigned int S, unsigned int N>
struct detrender_timing_thread : public timing_thread
{
    const int nfreq;
    const int nt_chunk;
    const int stride;
    const int niter = 16;

    float *intensity = nullptr;
    float *weights = nullptr;

    // A place to write dummy results, to keep the compiler from optimizing things out
    float *dummyp = nullptr;

    detrender_timing_thread(const shared_ptr<timing_thread_pool> &pool_, int nfreq_, int nt_chunk_, int stride_) :
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

	this->start_timer();
        for (int iter = 0; iter < niter; iter++)
	    _kernel_detrend_f<T,S,N> (nfreq, nt_chunk, intensity, weights, stride);
        this->stop_timer("kernel_detrend_f");

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


template<typename T, unsigned int S, unsigned int Nmax, typename std::enable_if<(Nmax==0),int>::type = 0>
inline std::thread make_detrender_timing_thread(const shared_ptr<timing_thread_pool> &pool, int polydeg, int nfreq, int nt_chunk, int stride)
{
    throw runtime_error("internal error in make_detrender_timing_thread()");
}

template<typename T, unsigned int S, unsigned int Nmax, typename std::enable_if<(Nmax>0),int>::type = 0>
inline std::thread make_detrender_timing_thread(const shared_ptr<timing_thread_pool> &pool, int polydeg, int nfreq, int nt_chunk, int stride)
{
    if (Nmax == polydeg + 1)
	return spawn_timing_thread< detrender_timing_thread<T,S,Nmax> > (pool, nfreq, nt_chunk, stride);
    return make_detrender_timing_thread<T,S,(Nmax-1)> (pool, polydeg, nfreq, nt_chunk, stride);
}


int main(int argc, char **argv)
{
    // (max polynomial degree) + 1
    static const int Nmax = 17;

    if (argc != 6) {
	cerr << "usage: time-detrenders <nfreq> <nt_chunk> <stride> <polydeg> <nthreads>\n";
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
        threads[i] = make_detrender_timing_thread<float,8,Nmax> (pool, polydeg, nfreq, nt_chunk, stride);
    for (int i = 0; i < nthreads; i++)
        threads[i].join();

    return 0;
}
