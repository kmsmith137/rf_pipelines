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
#include "kernels/intensity_clippers.hpp"

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
//
// clippers_wrms_vops: contains wrms kernels for a fixed tuple (T,S,Df,Dt,Iflag,Wflag).
//
// This extra level of virtual functions was introduced in order to reduce compile time!


template<typename T>
struct clipper_ops_base {
    const int S;
    const int Df;
    const int Dt;
    const bool Iflag;
    const bool Wflag;
    
    clipper_ops_base(int S_, int Df_, int Dt_, bool Iflag_, bool Wflag_) :
	S(S_), Df(Df_), Dt(Dt_), Iflag(Iflag_), Wflag(Wflag_) { }

    virtual ~clipper_ops_base() { }

    // fast_kernel2d
    //   - takes an input array of shape (nfreq, nt)
    //   - computes a scalar mean/rms, which is repeated redundantly S times
    //   - outputs downsampled arrays of shape (nfreq/Df, nt/Dt)
    //
    // fast_kernel1d_f
    //   - takes an input array of shape (nfreq, Dt*S)
    //   - computes mean/rms arrays of length S
    //   - outputs downsampled arrays of shape (nfreq/Df, S)

    virtual void apply_fast_wrms_kernel_2d(T *mean, T *rms, const T *intensity, const T *weights, int nfreq, int nt, int stride, T *ds_int, T *ds_wt) = 0;
    virtual void apply_fast_wrms_kernel_1d_f(T *mean, T *rms, const T *intensity, const T *weights, int nfreq, int stride, T *ds_int, T *ds_w) = 0;
    virtual void apply_fast_mask_kernel_2d(T *weights, const T *ds_intensity, T mean, T thresh, int nfreq, int nt, int stride, int ds_stride) = 0;
    virtual void apply_fast_mask_kernel_1d_f(T *weights, const T *ds_intensity, T mean, T thresh, int nfreq, int stride, int ds_stride) = 0;
    virtual void apply_fast_iterate_kernel_2d(T *mean, T *rms, const T *intensity, const T *weights, T in_mean, T in_thresh, int nfreq, int nt, int stride) = 0;
    
    void apply_reference_wrms_kernel_2d(T &mean, T &rms, const T *intensity, const T *weights, int nfreq, int nt, int stride, T *ds_int, T *ds_wt);
    void apply_reference_wrms_kernel_1d_f(T *mean, T *rms, const T *intensity, const T *weights, int nfreq, int stride, T *ds_int, T *ds_wt);    
    void apply_reference_mask_kernel_2d(T *weights, const T *ds_intensity, T mean, T thresh, int nfreq, int nt, int stride, int ds_stride);
    void apply_reference_mask_kernel_1d_f(T *weights, const T *ds_intensity, T mean, T thresh, int nfreq, int stride, int ds_stride);
    void apply_reference_iterate_kernel_2d(T &out_mean, T &out_rms, const T *intensity, const T *weights, T in_mean, T in_thresh, int nfreq, int nt, int stride);

    void test_wrms_kernel_2d(std::mt19937 &rng, int nfreq, int nt, int stride);
    void test_wrms_kernel_1d_f(std::mt19937 &rng, int nfreq, int nt, int stride);
    void test_mask_kernel_2d(std::mt19937 &rng, int nfreq, int nt, int stride, int ds_stride);
    void test_mask_kernel_1d_f(std::mt19937 &rng, int nfreq, int stride, int ds_stride);
    void test_iterate_kernel_2d(std::mt19937 &rng, int nfreq, int nt, int stride);

    void run_tests(std::mt19937 &rng);

    void _reference_wrms_kernel_helper(T &mean, T &rms, const T *intensity, const T *weights, int nfreq, int nt, int stride, T *ds_int, T *ds_wt, int ds_stride);
};


template<typename T>
void clipper_ops_base<T>::_reference_wrms_kernel_helper(T &mean, T &rms, const T *intensity, const T *weights, int nfreq, int nt, int stride, T *ds_int, T *ds_wt, int ds_stride)
{
    assert(nfreq % Df == 0);
    assert(nt % Dt == 0);
    assert(ds_stride >= nt/Dt);

    // double-precision here
    double acc0 = 0.0;
    double acc1 = 0.0;
    double acc2 = 0.0;

    for (int ifreq = 0; ifreq < nfreq; ifreq += Df) {
	for (int it = 0; it < nt; it += Dt) {
	    T wival = 0.0;
	    T wval = 0.0;

	    for (int jfreq = ifreq; jfreq < ifreq + Df; jfreq++) {
		for (int jt = it; jt < it + Dt; jt++) {
		    int s = jfreq*stride + jt;
		    wival += weights[s] * intensity[s];
		    wval += weights[s];
		}
	    }

	    acc0 += double(wval);
	    acc1 += double(wival);
	    acc2 += double(wival) * double(wival) / double(wval);
	    
	    *(ds_int++) = wival / wval;
	    *(ds_wt++) = wval;
	}

	// A little awkward but that's OK
	ds_int += (ds_stride - (nt/Dt));
	ds_wt += (ds_stride - (nt/Dt));
    }

    // FIXME case of invalid entries not tested
    mean = acc1/acc0;
    rms = sqrt(acc2/acc0 - mean*mean);
}


template<typename T>
void clipper_ops_base<T>::apply_reference_wrms_kernel_2d(T &mean, T &rms, const T *intensity, const T *weights, int nfreq, int nt, int stride, T *ds_int, T *ds_wt)
{
    // Call _reference_wrms_kernel_helper() with ds_stride = nt/Dt.
    _reference_wrms_kernel_helper(mean, rms, intensity, weights, nfreq, nt, stride, ds_int, ds_wt, nt/Dt);
}

template<typename T>
void clipper_ops_base<T>::apply_reference_wrms_kernel_1d_f(T *mean, T *rms, const T *intensity, const T *weights, int nfreq, int stride, T *ds_int, T *ds_wt)
{
    // Call _reference_wrms_kernel_helper() S times, with nt=Dt and ds_stride=S.
    for (int s = 0; s < S; s++)
	_reference_wrms_kernel_helper(mean[s], rms[s], intensity + s*Dt, weights + s*Dt, nfreq, Dt, stride, ds_int+s, ds_wt+s, S);
}


template<typename T>
void clipper_ops_base<T>::apply_reference_mask_kernel_2d(T *weights, const T *ds_intensity, T mean, T thresh, int nfreq, int nt, int stride, int ds_stride)
{
    assert(nfreq % Df == 0);
    assert(nt % Dt == 0);

    int nfreq_ds = nfreq / Df;
    int nt_ds = nt / Dt;

    for (int ifreq_ds = 0; ifreq_ds < nfreq_ds; ifreq_ds++) {
	for (int it_ds = 0; it_ds < nt_ds; it_ds++) {
	    T ival = ds_intensity[ifreq_ds*ds_stride + it_ds];

	    if (fabs(ival-mean) < thresh)
		continue;

	    for (int ifreq = ifreq_ds*Df; ifreq < (ifreq_ds+1)*Df; ifreq++)
		for (int it = it_ds*Dt; it < (it_ds+1)*Dt; it++)
		    weights[ifreq*stride+it] = 0.0;
	}
    }    
}


template<typename T>
void clipper_ops_base<T>::apply_reference_mask_kernel_1d_f(T *weights, const T *ds_intensity, T mean, T thresh, int nfreq, int stride, int ds_stride)
{
    apply_reference_mask_kernel_2d(weights, ds_intensity, mean, thresh, nfreq, Dt*S, stride, ds_stride);
}


template<typename T>
void clipper_ops_base<T>::apply_reference_iterate_kernel_2d(T &out_mean, T &out_rms, const T *intensity, const T *weights, T in_mean, T in_thresh, int nfreq, int nt, int stride)
{
    // double-precision here
    double acc0 = 0.0;
    double acc1 = 0.0;
    double acc2 = 0.0;

    for (int ifreq = 0; ifreq < nfreq; ifreq++) {
	for (int it = 0; it < nt; it++) {
	    T ival = intensity[ifreq*stride + it];
	    T wval = weights[ifreq*stride + it];
	    
	    if (fabs(ival - in_mean) >= in_thresh)
		continue;

	    acc0 += double(wval);
	    acc1 += double(wval) * double(ival);
	    acc2 += double(wval) * double(ival) * double(ival);
	}
    }

    // FIXME case of invalid entries not tested
    out_mean = acc1/acc0;
    out_rms = sqrt(acc2/acc0 - out_mean*out_mean);
}


template<typename T>
void clipper_ops_base<T>::test_wrms_kernel_2d(std::mt19937 &rng, int nfreq, int nt, int stride)
{
    assert(nfreq % Df == 0);
    assert(nt % (Dt*S) == 0);
    assert(stride >= nt);

    random_chunk rc(rng, nfreq, nt, stride);

    float ref_mean, ref_rms;
    vector<T> ref_ds_int((nfreq/Df) * (nt/Dt), -1.0);
    vector<T> ref_ds_wt((nfreq/Df) * (nt/Dt), -1.0);

    apply_reference_wrms_kernel_2d(ref_mean, ref_rms, rc.intensity, rc.weights, nfreq, nt, stride, &ref_ds_int[0], &ref_ds_wt[0]);

    vector<T> mean(S), rms(S);
    vector<T> ds_intv((nfreq/Df) * (nt/Dt), -1.0);
    vector<T> ds_wtv((nfreq/Df) * (nt/Dt), -1.0);

    T *ds_int = Iflag ? &ds_intv[0] : nullptr;
    T *ds_wt = Wflag ? &ds_wtv[0] : nullptr;

    apply_fast_wrms_kernel_2d(&mean[0], &rms[0], rc.intensity, rc.weights, nfreq, nt, stride, ds_int, ds_wt);
    
    T delta_mean = simd_helpers::compare(mean, vector<T>(S, ref_mean));
    T delta_rms = simd_helpers::compare(rms, vector<T>(S, ref_rms));

    if ((std::abs(delta_mean) > 1.0e-3 * Df*Dt) || (std::abs(delta_rms) > 1.0e-3 * sqrt(Df*Dt))) {
	cerr << "kernel_clip2d_wrms mean/rms mismatch:"
	     << " T=" << simd_helpers::type_name<T>() << ", S=" << S << ", Df=" << Df << ", Dt=" << Dt << ", Iflag=" << Iflag 
	     << ", Wflag=" << Wflag << ", nfreq=" << nfreq << ", nt=" << nt << ", stride=" << stride << "\n"
	     << "  mean: " << ref_mean << ", " << simd_helpers::vecstr(mean) << "\n"
	     << "  rms: " << ref_rms << ", " << simd_helpers::vecstr(rms) << "\n";
	exit(1);
    }

    if (Iflag && (simd_helpers::maxdiff(ref_ds_int, ds_intv) > 1.0e-3 * Df*Dt)) {
	cerr << "kernel_clip2d_wrms ds_int mismatch:"
	     << " T=" << simd_helpers::type_name<T>() << ", S=" << S << ", Df=" << Df << ", Dt=" << Dt << ", Iflag=" << Iflag 
	     << ", Wflag=" << Wflag << ", nfreq=" << nfreq << ", nt=" << nt << ", stride=" << stride << "\n";
	exit(1);
    }

    if (Wflag && (simd_helpers::maxdiff(ref_ds_wt, ds_wtv) > 1.0e-3 * Df*Dt)) {
	cerr << "kernel_clip2d_wrms ds_wt mismatch:"
	     << " T=" << simd_helpers::type_name<T>() << ", S=" << S << ", Df=" << Df << ", Dt=" << Dt << ", Iflag=" << Iflag 
	     << ", Wflag=" << Wflag << ", nfreq=" << nfreq << ", nt=" << nt << ", stride=" << stride << "\n";
	exit(1);
    }
}


template<typename T>
void clipper_ops_base<T>::test_wrms_kernel_1d_f(std::mt19937 &rng, int nfreq, int nt, int stride)
{
    assert(nfreq % Df == 0);
    assert(nt % (Dt*S) == 0);
    assert(stride >= nt);

    random_chunk rc(rng, nfreq, nt, stride);

    vector<T> ref_mean(S), ref_rms(S);
    vector<T> ref_ds_int((nfreq/Df) * S);
    vector<T> ref_ds_wt((nfreq/Df) * S);

    vector<T> mean(S), rms(S);
    vector<T> ds_intv((nfreq/Df) * S);
    vector<T> ds_wtv((nfreq/Df) * S);

    T *ds_int = Iflag ? &ds_intv[0] : nullptr;
    T *ds_wt = Wflag ? &ds_wtv[0] : nullptr;

    for (int it = 0; it < nt; it += Dt*S) {
	const T *icol = rc.intensity + it;
	const T *wcol = rc.weights + it;

	apply_reference_wrms_kernel_1d_f(&ref_mean[0], &ref_rms[0], icol, wcol, nfreq, stride, &ref_ds_int[0], &ref_ds_wt[0]);
	apply_fast_wrms_kernel_1d_f(&mean[0], &rms[0], icol, wcol, nfreq, stride, ds_int, ds_wt);
    
	T delta_mean = simd_helpers::compare(mean, ref_mean);
	T delta_rms = simd_helpers::compare(rms, ref_rms);
	
	if ((std::abs(delta_mean) > 1.0e-3 * Df*Dt) || (std::abs(delta_rms) > 2.0e-3 * sqrt(Df*Dt))) {
	    cerr << "kernel_clip1d_f_wrms mean/rms mismatch:"
		 << " T=" << simd_helpers::type_name<T>() << ", S=" << S << ", Df=" << Df << ", Dt=" << Dt << ", Iflag=" << Iflag 
		 << ", Wflag=" << Wflag << ", nfreq=" << nfreq << ", nt=" << nt << ", stride=" << stride << "\n"
		 << "  mean: " << simd_helpers::vecstr(ref_mean) << ", " << simd_helpers::vecstr(mean) << "\n"
		 << "  rms: " << simd_helpers::vecstr(ref_rms) << ", " << simd_helpers::vecstr(rms) << "\n";
	    exit(1);
	}
	
	if (Iflag && (simd_helpers::maxdiff(ref_ds_int, ds_intv) > 1.0e-3 * Df*Dt)) {
	    cerr << "kernel_clip1d_f_wrms ds_int mismatch:"
		 << " T=" << simd_helpers::type_name<T>() << ", S=" << S << ", Df=" << Df << ", Dt=" << Dt << ", Iflag=" << Iflag 
		 << ", Wflag=" << Wflag << ", nfreq=" << nfreq << ", nt=" << nt << ", stride=" << stride << "\n";
	    exit(1);
	}
	
	if (Wflag && (simd_helpers::maxdiff(ref_ds_wt, ds_wtv) > 1.0e-3 * Df*Dt)) {
	    cerr << "kernel_clip1d_f_wrms ds_wt mismatch:"
		 << " T=" << simd_helpers::type_name<T>() << ", S=" << S << ", Df=" << Df << ", Dt=" << Dt << ", Iflag=" << Iflag 
		 << ", Wflag=" << Wflag << ", nfreq=" << nfreq << ", nt=" << nt << ", stride=" << stride << "\n";
	    exit(1);
	}
    }
}


template<typename T>
void clipper_ops_base<T>::test_mask_kernel_2d(std::mt19937 &rng, int nfreq, int nt, int stride, int ds_stride)
{
    assert(nfreq % Df == 0);
    assert(nt % Dt == 0);
    assert(stride >= nt);
    assert(ds_stride >= (nt/Dt));

    int nfreq_ds = nfreq / Df;
    int nt_ds = nt / Dt;

    T mean = std::uniform_real_distribution<>()(rng);
    T thresh = std::uniform_real_distribution<>()(rng);

    vector<T> ds_intensity(nfreq_ds * ds_stride, 0.0);
    vector<T> weights = simd_helpers::uniform_randvec<T> (rng, nfreq * stride, 0.0, 1.0);
    vector<T> weights2 = weights;

    for (int ifreq_ds = 0; ifreq_ds < nfreq_ds; ifreq_ds++)
	for (int it_ds = 0; it_ds < nt_ds; it_ds++)
	    ds_intensity[ifreq_ds*ds_stride + it_ds] = mean + thresh * clip_rand(rng);

    apply_reference_mask_kernel_2d(&weights[0], &ds_intensity[0], mean, thresh, nfreq, nt, stride, ds_stride);
    apply_fast_mask_kernel_2d(&weights2[0], &ds_intensity[0], mean, thresh, nfreq, nt, stride, ds_stride);

    for (int ifreq = 0; ifreq < nfreq; ifreq++) {
	for (int it = 0; it < nt; it++) {
	    if (weights[ifreq*stride+it] == weights2[ifreq*stride+it])
		continue;

	    int ifreq_ds = ifreq / Df;
	    int it_ds = it / Dt;

	    cerr << "test_clip2d_mask failed:"
		 << " T=" << simd_helpers::type_name<T>() << ", S=" << S << ", Df=" << Df << ", Dt=" << Dt 
		 << ", nfreq=" << nfreq << ", nt=" << nt << ", stride=" << stride << ", ds_stride=" << ds_stride << "\n"
		 << "   at (ifreq,it)=(" << ifreq << "," << it << "): "
		 << " wt_ref=" << weights[ifreq*stride+it] 
		 << ", wt_fast=" << weights2[ifreq*stride+it] << "\n"
		 << "   mean=" << mean << ", thresh=" << thresh 
		 << ", ds_int=" << ds_intensity[ifreq_ds*ds_stride + it_ds] << "\n";

	    exit(1);
	}
    }
}


template<typename T>
void clipper_ops_base<T>::test_mask_kernel_1d_f(std::mt19937 &rng, int nfreq, int stride, int ds_stride)
{
    assert(nfreq % Df == 0);
    assert(stride >= (Dt*S));
    assert(ds_stride >= S);

    int nfreq_ds = nfreq / Df;
 
    T mean = std::uniform_real_distribution<>()(rng);
    T thresh = std::uniform_real_distribution<>()(rng);

    vector<T> weights = simd_helpers::uniform_randvec<T> (rng, nfreq * stride, 0.0, 1.0);
    vector<T> weights2 = weights;

    vector<T> ds_intensity(nfreq_ds * ds_stride, 0.0);

    for (int ifreq_ds = 0; ifreq_ds < nfreq_ds; ifreq_ds++)
	for (int s = 0; s < S; s++)
	    ds_intensity[ifreq_ds*ds_stride + s] = mean + thresh * clip_rand(rng);

    apply_reference_mask_kernel_1d_f(&weights[0], &ds_intensity[0], mean, thresh, nfreq, stride, ds_stride);
    apply_fast_mask_kernel_1d_f(&weights2[0], &ds_intensity[0], mean, thresh, nfreq, stride, ds_stride);

    for (int ifreq = 0; ifreq < nfreq; ifreq++) {
	for (int it = 0; it < Dt*S; it++) {
	    if (weights[ifreq*stride+it] == weights2[ifreq*stride+it])
		continue;

	    int ifreq_ds = ifreq / Df;
	    int it_ds = it / Dt;

	    cerr << "test_make_kernel_1d_f failed:"
		 << " T=" << simd_helpers::type_name<T>() << ", S=" << S << ", Df=" << Df << ", Dt=" << Dt 
		 << ", nfreq=" << nfreq << ", stride=" << stride << ", ds_stride=" << ds_stride << "\n"
		 << "   at (ifreq,it)=(" << ifreq << "," << it << "): "
		 << " wt_ref=" << weights[ifreq*stride+it] 
		 << ", wt_fast=" << weights2[ifreq*stride+it] << "\n"
		 << "   mean=" << mean << ", thresh=" << thresh 
		 << ", ds_int=" << ds_intensity[ifreq_ds*ds_stride + it_ds] << "\n";

	    exit(1);
	}
    }
}

		     
template<typename T>
void clipper_ops_base<T>::test_iterate_kernel_2d(std::mt19937 &rng, int nfreq, int nt, int stride)
{
    random_chunk rc(rng, nfreq, nt, stride);

    T in_mean = std::uniform_real_distribution<>()(rng);
    T in_thresh = std::uniform_real_distribution<>(1.0, 2.0)(rng);

    for (int ifreq = 0; ifreq < nfreq; ifreq++)
	for (int it = 0; it < nt; it++)
	    rc.intensity[ifreq*stride + it] = in_mean + in_thresh * clip_rand(rng);
    
    T ref_mean, ref_rms;
    apply_reference_iterate_kernel_2d(ref_mean, ref_rms, rc.intensity, rc.weights, in_mean, in_thresh, nfreq, nt, stride);

    vector<T> mean(S), rms(S);
    apply_fast_iterate_kernel_2d(&mean[0], &rms[0], rc.intensity, rc.weights, in_mean, in_thresh, nfreq, nt, stride);

    T delta_mean = simd_helpers::compare(mean, vector<T>(S, ref_mean));
    T delta_rms = simd_helpers::compare(rms, vector<T>(S, ref_rms));

    if ((std::abs(delta_mean) > 1.0e-3) || (std::abs(delta_rms) > 1.0e-3)) {
	cerr << "test_iterate_kernel2d failed:"
	     << " T=" << simd_helpers::type_name<T>() << ", S=" << S
	     << ", nfreq=" << nfreq << ", nt=" << nt << ", stride=" << stride << "\n"
	     << "  mean: " << ref_mean << ", " << simd_helpers::vecstr(mean) << "\n"
	     << "  rms: " << ref_rms << ", " << simd_helpers::vecstr(rms) << "\n";

	exit(1);
    }
}


template<typename T>
void clipper_ops_base<T>::run_tests(std::mt19937 &rng)
{
    int nfreq = Df * std::uniform_int_distribution<>(50,100)(rng);
    int nt = Dt * S * std::uniform_int_distribution<>(10,20)(rng);
    int stride = nt + std::uniform_int_distribution<>(0,4)(rng);

    test_wrms_kernel_2d(rng, nfreq, nt, stride);
    test_wrms_kernel_1d_f(rng, nfreq, nt, stride);
    test_mask_kernel_2d(rng, nfreq, nt, stride, (nt/Dt) + std::uniform_int_distribution<>(0,4)(rng));
    test_mask_kernel_1d_f(rng, nfreq, Dt*S + std::uniform_int_distribution<>(0,4)(rng), S + std::uniform_int_distribution<>(0,4)(rng));
    test_iterate_kernel_2d(rng, nfreq, nt, stride);
}


// -------------------------------------------------------------------------------------------------


template<typename T, int S_>
struct clipper_ops_base2 : public clipper_ops_base<T>
{
    clipper_ops_base2(int Df_, int Dt_, bool Iflag_, bool Wflag_)
	: clipper_ops_base<T> (S_, Df_, Dt_, Iflag_, Wflag_) { }

    // FIXME: this test was designed for a different version of the kernel API, and
    // looks strange when expressed in terms of the new kernels.

    virtual void apply_fast_iterate_kernel_2d(T *mean, T *rms, const T *intensity, const T *weights, T in_mean, T in_thresh, int nfreq, int nt, int stride) override
    {
	simd_t<T,S_> mean_x = in_mean;
	simd_t<T,S_> rms_x = in_thresh;

	_kernel_wrms_iterate_2d(mean_x, rms_x, intensity, weights, nfreq, nt, stride, 2, 1.0);

	mean_x.storeu(mean);
	rms_x.storeu(rms);
    }
};


template<typename T, int S_, int Df_, int Dt_>
struct clipper_ops_base3 : public clipper_ops_base2<T,S_>
{
    clipper_ops_base3(bool Iflag_, bool Wflag_) : clipper_ops_base2<T,S_>(Df_, Dt_, Iflag_, Wflag_) { }

    virtual void apply_fast_mask_kernel_2d(T *weights, const T *ds_intensity, T mean, T thresh, int nfreq, int nt, int stride, int ds_stride) override
    {
	_kernel_intensity_mask_2d<T,S_,Df_,Dt_> (weights, ds_intensity, simd_t<T,S_>(mean), simd_t<T,S_>(thresh), nfreq, nt, stride, ds_stride);
    }

    virtual void apply_fast_mask_kernel_1d_f(T *weights, const T *ds_intensity, T mean, T thresh, int nfreq, int stride, int ds_stride) override
    {
	_kernel_intensity_mask_1d_f<T,S_,Df_,Dt_> (weights, ds_intensity, simd_t<T,S_> (mean), simd_t<T,S_> (thresh), nfreq, stride, ds_stride);
    }
};


template<typename T, int S_, int Df_, int Dt_, bool Iflag_, bool Wflag_>
struct clipper_ops : clipper_ops_base3<T,S_,Df_,Dt_>
{
    clipper_ops() : clipper_ops_base3<T,S_,Df_,Dt_> (Iflag_, Wflag_) { }
    
    virtual void apply_fast_wrms_kernel_2d(T *mean, T *rms, const T *intensity, const T *weights, int nfreq, int nt, int stride, T *ds_int, T *ds_wt) override
    {
	simd_t<T,S_> mean_x;
	simd_t<T,S_> rms_x;

	// FIXME only TwoPass=false tested for now
	_kernel_noniterative_wrms_2d<T,S_,Df_,Dt_,Iflag_,Wflag_,false> (mean_x, rms_x, intensity, weights, nfreq, nt, stride, ds_int, ds_wt);

	mean_x.storeu(mean);
	rms_x.storeu(rms);
    }

    virtual void apply_fast_wrms_kernel_1d_f(T *mean, T *rms, const T *intensity, const T *weights, int nfreq, int stride, T *ds_int, T *ds_wt) override
    {
	simd_t<T,S_> mean_x;
	simd_t<T,S_> rms_x;

	// FIXME only TwoPass=false tested for now
	_kernel_noniterative_wrms_1d_f<T,S_,Df_,Dt_,Iflag_,Wflag_,false> (mean_x, rms_x, intensity, weights, nfreq, stride, ds_int, ds_wt);

	mean_x.storeu(mean);
	rms_x.storeu(rms);
    }
};


template<typename T, int S, int Df, int Dt, typename enable_if<((Df>1) || (Dt>1)),int>::type = 0>
inline void populate3(vector<shared_ptr<clipper_ops_base<T>>> &testsuite)
{
    testsuite.push_back(make_shared<clipper_ops<T,S,Df,Dt,false,false>> ());
    testsuite.push_back(make_shared<clipper_ops<T,S,Df,Dt,false,true>> ());
    testsuite.push_back(make_shared<clipper_ops<T,S,Df,Dt,true,false>> ());
    testsuite.push_back(make_shared<clipper_ops<T,S,Df,Dt,true,true>> ());    
}

template<typename T, int S, int Df, int Dt, typename enable_if<((Df==1) && (Dt==1)),int>::type = 0>
inline void populate3(vector<shared_ptr<clipper_ops_base<T>>> &testsuite)
{
    testsuite.push_back(make_shared<clipper_ops<T,S,Df,Dt,false,false>> ());
}

template<typename T, int S, int Df, int MaxDt, typename std::enable_if<(MaxDt==0),int>::type = 0>
inline void populate2(vector<shared_ptr<clipper_ops_base<T>>> &testsuite) { return; }

template<typename T, int S, int Df, int MaxDt, typename std::enable_if<(MaxDt > 0),int>::type = 0>
inline void populate2(vector<shared_ptr<clipper_ops_base<T>>> &testsuite)
{
    populate2<T,S,Df,(MaxDt/2)> (testsuite);
    populate3<T,S,Df,MaxDt> (testsuite);
}


template<typename T, int S, int MaxDf, int MaxDt, typename std::enable_if<(MaxDf==0),int>::type = 0>
inline void populate(vector<shared_ptr<clipper_ops_base<T>>> &testsuite) { return; }

template<typename T, int S, int MaxDf, int MaxDt, typename std::enable_if<(MaxDf > 0),int>::type = 0>
inline void populate(vector<shared_ptr<clipper_ops_base<T>>> &testsuite)
{
    populate<T,S,(MaxDf/2),MaxDt> (testsuite);
    populate2<T,S,MaxDf,MaxDt> (testsuite);
}


template<typename T, int S, int MaxDf, int MaxDt>
inline void run_all_clipper_tests(std::mt19937 &rng)
{
    vector<shared_ptr<clipper_ops_base<T>>> testsuite;
    populate<T,S,MaxDf,MaxDt> (testsuite);
    
    for (const auto &test: testsuite)
	test->run_tests(rng);
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

	// Run clipper tests to max (Df,Dt)=(16,16).
	// I wanted to use (32,32) but that ended up being pretty slow!

	run_all_clipper_tests<float,8,16,16> (rng);
    }

    cout << "test-kernels: all tests passed\n";
    return 0;
}
