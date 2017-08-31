// FIXME kernel unit tests in this file are currently incomplete!
//
// They were pretty complete at one point, but then stopped getting maintained,
// since test-cpp-python-equivalence.py became the "primary" test.  I still think
// the standalone tests are useful, so I'd like to bring them up-to-date with the
// current kernel code.  (This is an item in TODO_kernels.md.)

#include <memory>
#include <cassert>
#include <stdexcept>

#include <simd_helpers/simd_debug.hpp>
#include "rf_pipelines_internals.hpp"

#include "kernels/polyfit.hpp"

using namespace std;
using namespace rf_pipelines;


// -------------------------------------------------------------------------------------------------
//
// General-purpose helpers


// Generates a random number in the range [-2,2], but not too close to (+/- 1).
// This is useful when testing clippers, to avoid spurious roundoff-indiced 
// differences between reference code and fast code.
inline double clip_rand(std::mt19937 &rng)
{
    for (;;) {
	double t = std::uniform_real_distribution<>(-2.,2.)(rng);
	double u = fabs(fabs(t)-1.0);
	if (u > 1.0e-3)
	    return t;
    }
}


// Fills length-n 1D strided array with values of a randomly generated polynomial.
template<typename T>
inline void randpoly(T *dst, std::mt19937 &rng, int deg, int n, int stride)
{
    vector<T> coeffs = simd_helpers::gaussian_randvec<T> (rng, deg+1);

    for (int i = 0; i < n; i++) {
	T t = T(i) / T(n);
	T tp = 1.0;
	T y = 0.0;

	for (int p = 0; p <= deg; p++) {
	    y += coeffs[p] * tp;
	    tp *= t;
	}

	dst[i*stride] = y;
    }
}


// Makes length-n strided weights array badly conditioned, assuming a polynomial fit of degree 'deg'.
template<typename T>
inline void make_weights_badly_conditioned(T *dst, std::mt19937 &rng, int deg, int n, int stride)
{
    for (int i = 0; i < n; i++)
	dst[i*stride] = 0;

    for (int i = 0; i < deg; i++) {
	int j = std::uniform_int_distribution<>(0,n-1)(rng);
	dst[j*stride] = 1.0;
    }
}


template<typename T>
inline bool is_all_zero(const T *vec, int n, int stride)
{
    for (int i = 0; i < n; i++)
	if (vec[i*stride] != 0.0)
	    return false;
    return true;
}


template<typename T>
inline bool is_all_positive(const T *vec, int n, int stride)
{
    for (int i = 0; i < n; i++)
	if (vec[i*stride] <= 0.0)
	    return false;
    return true;
}


// hconst_simd_ntuple<T,S,N> (const T *p):    constructs "horizontally constant" N-tuple from a length N array
// hconst_simd_trimatrix<T,S,N> (const T *p): constructs "horizontally constant" N-tuple from a length N(N+1)/2 array
//
// "Horizontally constant" means "constant within each simd_t".


template<typename T, int S, int N, typename std::enable_if<(N==0),int>::type = 0>
inline simd_ntuple<T,S,N> hconst_simd_ntuple(const T *p) { return simd_ntuple<T,S,0> (); }

template<typename T, int S, int N, typename std::enable_if<(N==0),int>::type = 0>
inline simd_trimatrix<T,S,N> hconst_simd_trimatrix(const T *p) { return simd_trimatrix<T,S,0> (); }

template<typename T, int S, int N, typename std::enable_if<(N>0),int>::type = 0>
inline simd_ntuple<T,S,N> hconst_simd_ntuple(const T *p)
{ 
    return simd_ntuple<T,S,N> (hconst_simd_ntuple<T,S,N-1>(p), simd_t<T,S>(p[N-1]));
}

template<typename T, int S, int N, typename std::enable_if<(N>0),int>::type = 0>
inline simd_trimatrix<T,S,N> hconst_simd_trimatrix(const T *p) 
{
    return simd_trimatrix<T,S,N> (hconst_simd_trimatrix<T,S,N-1>(p), hconst_simd_ntuple<T,S,N>(p + (N*(N-1))/2));
}


// -------------------------------------------------------------------------------------------------


struct random_chunk {
    const int nfreq;
    const int nt;
    const int stride;
    
    float *intensity = nullptr;
    float *weights = nullptr;

    random_chunk(std::mt19937 &rng, int nfreq, int nt, int stride);
    random_chunk(std::mt19937 &rng, int nfreq, int nt);
    ~random_chunk();
    
    // noncopyable
    random_chunk(const random_chunk &) = delete;
    random_chunk &operator=(const random_chunk &) = delete;
};


random_chunk::random_chunk(std::mt19937 &rng, int nfreq_, int nt_, int stride_) :
    nfreq(nfreq_), nt(nt_), stride(stride_)
{
    assert(nfreq > 0);
    assert(nt > 0);
    assert(stride >= nt);

    intensity = aligned_alloc<float> (nfreq * stride);
    weights = aligned_alloc<float> (nfreq * stride);

    std::normal_distribution<> gdist;

    for (int i = 0; i < nfreq * stride; i++) {
	intensity[i] = gdist(rng) + 1.0;
	weights[i] = std::uniform_real_distribution<>()(rng);
    }
}


random_chunk::random_chunk(std::mt19937 &rng, int nfreq_, int nt_) :
    random_chunk(rng, nfreq_, nt_, nt_ + std::uniform_int_distribution<>(0,4)(rng))
{ }


random_chunk::~random_chunk()
{
    free(intensity);
    free(weights);
    intensity = weights = nullptr;
}


// -------------------------------------------------------------------------------------------------
//
// Test _kernel_legpoly_eval().


template<typename T>
static vector<T> reference_legpoly_eval(int npl, const vector<T> &zvec)
{
    assert(npl > 0);
    assert(zvec.size() > 0);

    int nz = zvec.size();
    vector<T> out_pl(npl * nz);

    for (int iz = 0; iz < nz; iz++)
	out_pl[iz] = 1.0;

    if (npl <= 1)
	return out_pl;

    for (int iz = 0; iz < nz; iz++)
	out_pl[nz+iz] = zvec[iz];

    for (int l = 2; l < npl; l++) {
	T a = (2*l-1) / T(l);
	T b = -(l-1) / T(l);
	
	for (int iz = 0; iz < nz; iz++)
	    out_pl[l*nz + iz] = a * zvec[iz] * out_pl[(l-1)*nz + iz] + b * out_pl[(l-2)*nz + iz];
    }

    return out_pl;
}


template<typename T, int S, int N>
static void test_legpoly_eval(std::mt19937 &rng)
{
    simd_t<T,S> z = simd_helpers::uniform_random_simd_t<T,S> (rng, -1.0, 1.0);

    simd_ntuple<T,S,N> pl;
    _kernel_legpoly_eval(pl, z);

    vector<T> pl0 = reference_legpoly_eval(N, vectorize(z));

#if 0
    for (int iz = 0; iz < S; iz++) {
	cout << z0[iz] << ": ";
	for (int l = 0; l < N; l++)
	    cout << " " << pl0[l*S+iz];
	cout << "\n";
    }
#endif

    T epsilon = simd_helpers::compare(vectorize(pl), pl0);
    assert(epsilon < 1.0e-5);
}


// -------------------------------------------------------------------------------------------------
//
// Test _kernel_detrend_t_pass1()


template<typename T>
static void reference_detrend_t_pass1(T *outm, T *outv, int npl, int nt, const T *ivec, const T *wvec)
{
    vector<T> tmp_z(nt);
    for (int it = 0; it < nt; it++)
	tmp_z[it] = 2 * (it+0.5) / T(nt) - 1;

    vector<T> tmp_pl = reference_legpoly_eval(npl, tmp_z);

    vector<T> tmp_wp(npl * nt);
    for (int l = 0; l < npl; l++)
	for (int it = 0; it < nt; it++)
	    tmp_wp[l*nt+it] = wvec[it] * tmp_pl[l*nt+it];

    for (int l = 0; l < npl; l++) {
	for (int l2 = 0; l2 <= l; l2++) {
	    T t = 0.0;
	    for (int it = 0; it < nt; it++)
		t += tmp_wp[l*nt+it] * tmp_pl[l2*nt+it];

	    outm[(l*(l+1))/2 + l2] = t;
	}

	T t = 0.0;
	for (int it = 0; it < nt; it++)
	    t += tmp_wp[l*nt+it] * ivec[it];

	outv[l] = t;
    }
}


template<typename T, int S, int N>
static void test_detrend_t_pass1(std::mt19937 &rng, int nt)
{
    constexpr int NN = (N*(N+1))/2;

    vector<T> ivec = simd_helpers::uniform_randvec<T> (rng, nt, 0.0, 1.0);
    vector<T> wvec = simd_helpers::uniform_randvec<T> (rng, nt, 0.0, 1.0);

    simd_trimatrix<T,S,N> outm;
    simd_ntuple<T,S,N> outv;

    _kernel_detrend_t_pass1(outm, outv, nt, &ivec[0], &wvec[0]);

    vector<T> outm0(NN);
    vector<T> outv0(N);

    reference_detrend_t_pass1(&outm0[0], &outv0[0], N, nt, &ivec[0], &wvec[0]);

    simd_trimatrix<T,S,N> outm1 = hconst_simd_trimatrix<T,S,N> (&outm0[0]);
    simd_ntuple<T,S,N> outv1 = hconst_simd_ntuple<T,S,N> (&outv0[0]);

    T epsilon_m = simd_helpers::compare(vectorize(outm), vectorize(outm1));
    T epsilon_v = simd_helpers::compare(vectorize(outv), vectorize(outv1));
    
    assert(epsilon_m < 1.0e-5);    
    assert(epsilon_v < 1.0e-5);    
}



// -------------------------------------------------------------------------------------------------
//
// Test _kernel_detrend_t_pass2()


template<typename T>
static void reference_detrend_t_pass2(T *ivec, int npl, int nt, const T *coeffs)
{
    vector<T> tmp_z(nt);
    for (int it = 0; it < nt; it++)
	tmp_z[it] = 2 * (it+0.5) / T(nt) - 1;

    vector<T> tmp_pl = reference_legpoly_eval(npl, tmp_z);

    for (int l = 0; l < npl; l++)
	for (int it = 0; it < nt; it++)
	    ivec[it] -= coeffs[l] * tmp_pl[l*nt + it];
}


template<typename T, int S, int N>
static void test_detrend_t_pass2(std::mt19937 &rng, int nt)
{
    vector<T> coeffs = simd_helpers::gaussian_randvec<T> (rng, N);
    vector<T> ivec = simd_helpers::gaussian_randvec<T> (rng, nt);
    vector<T> ivec2 = ivec;
    
    _kernel_detrend_t_pass2<T,S,N> (&ivec[0], nt, hconst_simd_ntuple<T,S,N> (&coeffs[0]));
    reference_detrend_t_pass2(&ivec2[0], N, nt, &coeffs[0]);

    T epsilon = simd_helpers::compare(ivec, ivec2);
    assert(epsilon < 1.0e-5);
}


// -------------------------------------------------------------------------------------------------
//
// Some general tests on kernel_detrend_t: 
//   "nulling": detrending a polynmomial should give zero,
//   "idempotency": detrending twice should be the same as detrending once
//
// Note: the idempotency test also tests masking, by randomly choosing some
// rows to make badly conditioned.


template<typename T, int S, int N>
static void test_detrend_t_nulling(std::mt19937 &rng, int nfreq, int nt, int stride)
{
    vector<T> intensity(nfreq * stride, 0.0);
    vector<T> weights = simd_helpers::uniform_randvec<T> (rng, nfreq * stride, 0.0, 1.0);

    for (int ifreq = 0; ifreq < nfreq; ifreq++)
	randpoly(&intensity[ifreq*stride], rng, N-1, nt, 1);

    _kernel_detrend_t<T,S,N> (nfreq, nt, &intensity[0], &weights[0], stride);

    double epsilon = simd_helpers::maxabs(intensity);

    if (epsilon > 1.0e-3) {
	cerr << "test_detrend_t_nulling failed (N=" << N << "): epsilon=" << epsilon << "\n";
	exit(1);
    }
}


template<typename T, int S, int N>
static void test_detrend_t_idempotency(std::mt19937 &rng, int nfreq, int nt, int stride)
{
    vector<T> intensity = simd_helpers::uniform_randvec<T> (rng, nfreq * stride, 0.0, 1.0);
    vector<T> weights = simd_helpers::uniform_randvec<T> (rng, nfreq * stride, 0.0, 1.0);

    _kernel_detrend_t<T,S,N> (nfreq, nt, &intensity[0], &weights[0], stride);

    // Give each row a 50% chance of being well-conditioned.
    vector<bool> well_conditioned(nfreq, true);
    for (int ifreq = 0; ifreq < nfreq; ifreq++) {
	if (std::uniform_real_distribution<>()(rng) > 0.5) {
	    well_conditioned[ifreq] = false;
	    make_weights_badly_conditioned(&weights[ifreq*stride], rng, N-1, nt, 1);
	}
    }

    vector<T> intensity2 = intensity;
    _kernel_detrend_t<T,S,N> (nfreq, nt, &intensity2[0], &weights[0], stride);

    for (int ifreq = 0; ifreq < nfreq; ifreq++) {
	if (well_conditioned[ifreq] && !is_all_positive(&weights[ifreq*stride], nt, 1)) {
	    cerr << "test_detrend_t_idempotency failed(N=" << N << ": well-conditioned weights were incorrectly masked\n";
	    exit(1);
	}

	if (!well_conditioned[ifreq] && !is_all_zero(&weights[ifreq*stride], nt, 1)) {
	    cerr << "test_detrend_t_idempotency failed(N=" << N << ": poorly conditioned weights did not get masked\n";
	    exit(1);
	}
    }

    double epsilon = simd_helpers::maxdiff(intensity, intensity2);
    
    if (epsilon > 1.0e-3) {
	cerr << "test_detrend_t_idempotency failed (N=" << N << "): epsilon=" << epsilon << endl;
	exit(1);
    }
}


// -------------------------------------------------------------------------------------------------
//
// "Nulling" and "idempotency" tests for kernel_detrend_f (analogous to tests for kernel_detrend_t above)
//
// Reminder: the idempotency test also tests masking, by randomly choosing some
// rows to make badly conditioned.


template<typename T, int S, int N>
static void test_detrend_f_nulling(std::mt19937 &rng, int nfreq, int nt, int stride)
{
    vector<T> intensity(nfreq * stride, 0.0);
    vector<T> weights = simd_helpers::uniform_randvec<T> (rng, nfreq * stride, 0.0, 1.0);

    for (int it = 0; it < nt; it++)
	randpoly(&intensity[it], rng, N-1, nfreq, stride);

    _kernel_detrend_f<T,S,N> (nfreq, nt, &intensity[0], &weights[0], stride);

    double epsilon = simd_helpers::maxabs(intensity);

    if (epsilon > 1.0e-3) {
	cerr << "test_detrend_f_nulling failed (N=" << N << "): epsilon=" << epsilon << "\n";
	exit(1);
    }
}


template<typename T, int S, int N>
static void test_detrend_f_idempotency(std::mt19937 &rng, int nfreq, int nt, int stride)
{
    vector<T> intensity = simd_helpers::uniform_randvec<T> (rng, nfreq * stride, 0.0, 1.0);
    vector<T> weights = simd_helpers::uniform_randvec<T> (rng, nfreq * stride, 0.0, 1.0);

    _kernel_detrend_f<T,S,N> (nfreq, nt, &intensity[0], &weights[0], stride);

    vector<bool> well_conditioned(nt, true);

    // Assign badly conditioned columns, by looping over S-column blocks
    for (int it = 0; it < nt; it += S) {
	// No badly conditioned columns in this block
	if (std::uniform_real_distribution<>()(rng) < 0.33)
	    continue;

	// If this flag is set, the block will be all badly conditioned
	bool all_badly_conditioned = (std::uniform_real_distribution<>()(rng) < 0.5);

	for (int jt = it; jt < it+S; jt++) {
	    if (all_badly_conditioned || (std::uniform_real_distribution<>()(rng) < 0.5)) {
		well_conditioned[jt] = false;
		make_weights_badly_conditioned(&weights[jt], rng, N-1, nfreq, stride);
	    }
	}
    }

    vector<T> intensity2 = intensity;
    _kernel_detrend_f<T,S,N> (nfreq, nt, &intensity2[0], &weights[0], stride);

    for (int it = 0; it < nt; it++) {
	if (well_conditioned[it] && !is_all_positive(&weights[it], nfreq, stride)) {
	    cerr << "test_detrend_f_idempotency failed(N=" << N << "): well-conditioned weights were incorrectly masked\n";
	    exit(1);
	}

	if (!well_conditioned[it] && !is_all_zero(&weights[it], nfreq, stride)) {
	    cerr << "test_detrend_f_idempotency failed(N=" << N << "): poorly conditioned weights did not get masked\n";
	    exit(1);
	}
    }
    
    double epsilon = simd_helpers::maxdiff(intensity, intensity2);

    if (epsilon > 1.0e-3) {
	cerr << "test_detrend_f_idempotency failed (N=" << N << "): epsilon=" << epsilon << endl;
	exit(1);
    }
}


// -------------------------------------------------------------------------------------------------


template<typename T, int S, int N>
static void test_detrend_transpose(std::mt19937 &rng, int n1, int n2, int stride1, int stride2)
{
    vector<T> intensity12(n1 * stride2, 0.0);
    vector<T> intensity21(n2 * stride1, 0.0);

    vector<T> weights12(n1 * stride2, 0.0);
    vector<T> weights21(n2 * stride1, 0.0);

    std::normal_distribution<> dist;

    for (int i = 0; i < n1; i++) {
	for (int j = 0; j < n2; j++) {
	    intensity12[i*stride2+j] = intensity21[j*stride1+i] = dist(rng);
	    weights12[i*stride2+j] = weights21[j*stride1+i] = std::uniform_real_distribution<>()(rng);
	}
    }

    _kernel_detrend_t<T,S,N> (n1, n2, &intensity12[0], &weights12[0], stride2);
    _kernel_detrend_f<T,S,N> (n2, n1, &intensity21[0], &weights21[0], stride1);

    T epsilon = 0;

    for (int i = 0; i < n1; i++) {
	for (int j = 0; j < n2; j++) {
	    T x = intensity12[i*stride2+j];
	    T y = intensity21[j*stride1+i];
	    epsilon = std::max(epsilon, std::fabs(x-y));
	}
    }

    if (epsilon > 1.0e-4) {
	cerr << "test_detrend_transpose failed: T=" << simd_helpers::type_name<T>() << ", S=" << S
	     << ", N=" << N << ", n1=" << n1 << ", n2=" << n2 << ", stride1=" << stride1
	     << ", stride2=" << stride2 << ": epsilon=" << epsilon << endl;
	exit(1);
    }
}


// -------------------------------------------------------------------------------------------------


template<typename T, int S, int Nmax, typename std::enable_if<(Nmax==0),int>::type = 0>
static void run_all_polynomial_detrender_tests(std::mt19937 &rng)
{
    return;
}


template<typename T, int S, int Nmax, typename std::enable_if<(Nmax>0),int>::type = 0>
static void run_all_polynomial_detrender_tests(std::mt19937 &rng)
{
    run_all_polynomial_detrender_tests<T,S,(Nmax-1)> (rng);

    for (int iter = 0; iter < 10; iter++) {
	int nfreq = std::uniform_int_distribution<>(10*Nmax,20*Nmax)(rng);
	int nt = S * std::uniform_int_distribution<>((10*Nmax)/S,(20*Nmax)/S)(rng);
	int stride = nt + S * std::uniform_int_distribution<>(0,4)(rng);

	test_legpoly_eval<T,S,Nmax> (rng);

	test_detrend_t_pass1<T,S,Nmax> (rng, nt);
	test_detrend_t_pass2<T,S,Nmax> (rng, nt);
	test_detrend_t_nulling<T,S,Nmax> (rng, nfreq, nt, stride);
	test_detrend_t_idempotency<T,S,Nmax> (rng, nfreq, nt, stride);

	test_detrend_f_nulling<T,S,Nmax> (rng, nfreq, nt, stride);
	test_detrend_f_idempotency<T,S,Nmax> (rng, nfreq, nt, stride);

	// n2, stride2 only used in test_detrend_transpose()
	int n2 = S * std::uniform_int_distribution<>((10*Nmax)/S,(20*Nmax)/S)(rng);
	int stride2 = n2 + S * std::uniform_int_distribution<>(0,4)(rng);
	test_detrend_transpose<T,S,Nmax> (rng, nt, n2, stride, stride2);
    }
}


// -------------------------------------------------------------------------------------------------


int main(int argc, char **argv)
{
    static const int num_iterations = 20;

    std::random_device rd;
    std::mt19937 rng(rd());

    for (int iter = 1; iter <= num_iterations; iter++) {
	cout << "iteration " << iter << "/" << num_iterations << "\n";

	// FIXME we only go to polynomial degree 8 because the detrender tests can
	// fail for higher degrees ("poorly conditioned weights did not get masked").
	//
	// Maybe it makes sense to implement double-precision polynomial detrending?
	// Or maybe there is a better way for the Cholesky inversion to detect numerical
	// instability?  (See TODO_kernels.md)

	run_all_polynomial_detrender_tests<float,8,9> (rng);   // max degree 8
    }

    cout << "test-kernels: all tests passed\n";
    return 0;
}
