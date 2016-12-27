// Note: this file should only include kernel headers (kernels/*.hpp), not toplevel headers (./*.hpp)
// (This is because it has a Makefile dependency on the former but not the latter.)

#include <cassert>
#include "simd_helpers/simd_debug.hpp"
#include "kernels/clip2d.hpp"

using namespace std;
using namespace rf_pipelines;


// -------------------------------------------------------------------------------------------------
//
// General-purpose helpers


// Cut-and-pasted from rf_pipelines_internals.hpp, since we don't want to #include it.
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


// Frequently used in conjunction with simd_helpers::vectorize()
template<typename T> 
inline double maxabs(const std::vector<T> &v)
{
    assert(v.size() > 0);

    double ret = fabs(v[0]);
    for (unsigned int i = 1; i < v.size(); i++)
	ret = max(ret, fabs(double(v[i])));

    return ret;
}


template<typename T> 
inline double maxdiff(int n, const T *v1, const T *v2)
{
    assert(n > 0);

    double ret = fabs(v1[0]-v2[0]);
    for (int i = 1; i < n; i++)
	ret = max(ret, fabs(double(v1[i]-v2[i])));

    return ret;
}


// Generates a random number in the range [-2,2], but not too close to (+/- 1).
// This is useful when 
inline double clip_rand(std::mt19937 &rng)
{
    for (;;) {
	double t = std::uniform_real_distribution<>(-2.,2.)(rng);
	double u = fabs(fabs(t)-1.0);
	if (u > 1.0e-3)
	    return t;
    }
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


template<typename T>
static void reference_clip2d_wrms(T &mean, T &rms, const T *intensity, const T *weights, int nfreq, int nt, int stride, int nds_f, int nds_t, T *ds_int, T *ds_wt)
{
    assert(nfreq % nds_f == 0);
    assert(nt % nds_t == 0);

    // double-precision here
    double acc0 = 0.0;
    double acc1 = 0.0;
    double acc2 = 0.0;
    
    for (int ifreq = 0; ifreq < nfreq; ifreq += nds_f) {
	for (int it = 0; it < nt; it += nds_t) {
	    T wival = 0.0;
	    T wval = 0.0;

	    for (int jfreq = ifreq; jfreq < ifreq + nds_f; jfreq++) {
		for (int jt = it; jt < it + nds_t; jt++) {
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
    }

    // FIXME case of invalid entries not tested
    mean = acc1/acc0;
    rms = sqrt(acc2/acc0 - mean*mean);
}


template<typename T, unsigned int S>
void test_clip2d_wrms_postmortem(int Df, int Dt, int nfreq, int nt, int stride, T ref_mean, 
				 T ref_rms, T *ref_ds_int, T *ref_ds_wt, simd_t<T,S> mean, 
				 simd_t<T,S> rms, T *ds_int, T *ds_wt)
{
    stringstream ss;
    ss << "test_clip2d_wrms failed:"
       << " T=" << simd_helpers::type_name<T>()
       << ", S=" << S << ", Df=" << Df << ", Dt=" << Dt 
       << ", nfreq=" << nfreq << ", nt=" << nt << ", stride=" << stride;

    vector<float> delta1 = vectorize(mean - simd_t<T,S> (ref_mean));
    vector<float> delta2 = vectorize(rms - simd_t<T,S> (ref_rms));    
    int n = (nfreq * nt) / (Df * Dt);

    if ((maxabs(delta1) > 1.0e-3 * Df*Dt) || (maxabs(delta2) > 1.0e-3 * sqrt(Df*Dt))) {
	cerr << ss.str() << "\n"
	     << "  mean: " << ref_mean << ", " << mean << "\n"
	     << "  rms: " << ref_rms << ", " << rms << "\n";
	exit(1);
    }

    if (ds_int && (maxdiff(n,ref_ds_int,ds_int) > 1.0e-3 * Df*Dt)) {
	cerr << ss.str() << ": ds_int arrays differ\n";
	exit(1);
    }

    if (ds_wt && (maxdiff(n,ref_ds_wt,ds_wt) > 1.0e-3 * Df*Dt)) {
	cerr << ss.str() << ": ds_wt arrays differ\n";
	exit(1);
    }
}


template<typename T, unsigned int S, unsigned int Df, unsigned int Dt>
static void test_clip2d_wrms(std::mt19937 &rng, int nfreq, int nt, int stride)
{
    assert(nfreq % Df == 0);
    assert(nt % (Dt*S) == 0);
    assert(stride >= nt);

    random_chunk rc(rng, nfreq, nt, stride);

    float ref_mean, ref_rms;
    vector<T> ref_ds_int((nfreq/Df) * (nt/Dt), -1.0);
    vector<T> ref_ds_wt((nfreq/Df) * (nt/Dt), -1.0);

    reference_clip2d_wrms(ref_mean, ref_rms, rc.intensity, rc.weights, nfreq, nt, rc.stride, Df, Dt, &ref_ds_int[0], &ref_ds_wt[0]);
    
    simd_t<T,S> mean, rms;
    vector<T> ds_int((nfreq/Df) * (nt/Dt), -1.0);
    vector<T> ds_wt((nfreq/Df) * (nt/Dt), -1.0);
    
    _kernel_clip2d_wrms<T,S,Df,Dt,false,false> (mean, rms, rc.intensity, rc.weights, nfreq, nt, rc.stride, nullptr, nullptr);
    test_clip2d_wrms_postmortem(Df, Dt, nfreq, nt, stride, ref_mean, ref_rms, &ref_ds_int[0], &ref_ds_wt[0], mean, rms, (T *)nullptr, (T *)nullptr);

    _kernel_clip2d_wrms<T,S,Df,Dt,false,true> (mean, rms, rc.intensity, rc.weights, nfreq, nt, rc.stride, nullptr, &ds_wt[0]);
    test_clip2d_wrms_postmortem(Df, Dt, nfreq, nt, stride, ref_mean, ref_rms, &ref_ds_int[0], &ref_ds_wt[0], mean, rms, (T *)nullptr, &ds_wt[0]);

    _kernel_clip2d_wrms<T,S,Df,Dt,true,false> (mean, rms, rc.intensity, rc.weights, nfreq, nt, rc.stride, &ds_int[0], nullptr);
    test_clip2d_wrms_postmortem(Df, Dt, nfreq, nt, stride, ref_mean, ref_rms, &ref_ds_int[0], &ref_ds_wt[0], mean, rms, &ds_int[0], (T *)nullptr);

    _kernel_clip2d_wrms<T,S,Df,Dt,true,true> (mean, rms, rc.intensity, rc.weights, nfreq, nt, rc.stride, &ds_int[0], &ds_wt[0]);
    test_clip2d_wrms_postmortem(Df, Dt, nfreq, nt, stride, ref_mean, ref_rms, &ref_ds_int[0], &ref_ds_wt[0], mean, rms, &ds_int[0], &ds_wt[0]);
}


template<typename T, unsigned int S, unsigned int Df, unsigned int Dt>
static void test_clip2d_wrms(std::mt19937 &rng)
{
    int nfreq = Df * std::uniform_int_distribution<>(10,20)(rng);
    int nt = Dt * S * std::uniform_int_distribution<>(10,20)(rng);
    int stride = nt + std::uniform_int_distribution<>(0,4)(rng);

    test_clip2d_wrms<T,S,Df,Dt> (rng, nfreq, nt, stride);
}


// Fixed Df, many Dt
template<typename T, unsigned int S, unsigned int Df, unsigned int MaxDt, typename std::enable_if<(MaxDt==1),int>::type = 0>
static void test_clip2d_wrms_varying_dt(std::mt19937 &rng)
{
    test_clip2d_wrms<T,S,Df,1> (rng);
}

// Fixed Df, many Dt
template<typename T, unsigned int S, unsigned int Df, unsigned int MaxDt, typename std::enable_if<(MaxDt>1),int>::type = 0>
static void test_clip2d_wrms_varying_dt(std::mt19937 &rng)
{
    test_clip2d_wrms_varying_dt<T,S,Df,MaxDt/2> (rng);
    test_clip2d_wrms<T,S,Df,MaxDt> (rng);
}

// Many Df, many Dt
template<typename T, unsigned int S, unsigned int MaxDf, unsigned int MaxDt, typename std::enable_if<(MaxDf==1),int>::type = 0>
static void test_clip2d_wrms_all(std::mt19937 &rng)
{
    test_clip2d_wrms_varying_dt<T,S,1,MaxDt> (rng);
}

// Many Df, many Dt
template<typename T, unsigned int S, unsigned int MaxDf, unsigned int MaxDt, typename std::enable_if<(MaxDf>1),int>::type = 0>
static void test_clip2d_wrms_all(std::mt19937 &rng)
{
    test_clip2d_wrms_all<T,S,MaxDf/2,MaxDt> (rng);
    test_clip2d_wrms_varying_dt<T,S,MaxDf,MaxDt> (rng);
}


// -------------------------------------------------------------------------------------------------


template<typename T>
static void reference_clip2d_mask(T *weights, const T *ds_intensity, T mean, T thresh, int nfreq, int nt, int stride, int Df, int Dt, int ds_stride)
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


template<typename T, unsigned int S, unsigned int Df, unsigned int Dt>
static void test_clip2d_mask(std::mt19937 &rng, int nfreq, int nt, int stride, int ds_stride)
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

    reference_clip2d_mask(&weights[0], &ds_intensity[0], mean, thresh, nfreq, nt, stride, Df, Dt, ds_stride);
    _kernel_clip2d_mask<T,S,Df,Dt> (&weights2[0], &ds_intensity[0], simd_t<T,S>(mean), simd_t<T,S>(thresh), nfreq, nt, stride, ds_stride);

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


template<typename T, unsigned int S, unsigned int Df, unsigned int Dt>
static void test_clip2d_mask(std::mt19937 &rng)
{
    int nfreq = Df * std::uniform_int_distribution<>(10,20)(rng);
    int nt = Dt * S * std::uniform_int_distribution<>(10,20)(rng);
    int stride = nt + std::uniform_int_distribution<>(0,4)(rng);
    int ds_stride = (nt/Dt) + std::uniform_int_distribution<>(0,4)(rng);

    test_clip2d_mask<T,S,Df,Dt> (rng, nfreq, nt, stride, ds_stride);
}


// Fixed Df, many Dt
template<typename T, unsigned int S, unsigned int Df, unsigned int MaxDt, typename std::enable_if<(MaxDt==1),int>::type = 0>
static void test_clip2d_mask_varying_dt(std::mt19937 &rng)
{
    test_clip2d_mask<T,S,Df,1> (rng);
}

// Fixed Df, many Dt
template<typename T, unsigned int S, unsigned int Df, unsigned int MaxDt, typename std::enable_if<(MaxDt>1),int>::type = 0>
static void test_clip2d_mask_varying_dt(std::mt19937 &rng)
{
    test_clip2d_mask_varying_dt<T,S,Df,MaxDt/2> (rng);
    test_clip2d_mask<T,S,Df,MaxDt> (rng);
}

// Many Df, many Dt
template<typename T, unsigned int S, unsigned int MaxDf, unsigned int MaxDt, typename std::enable_if<(MaxDf==1),int>::type = 0>
static void test_clip2d_mask_all(std::mt19937 &rng)
{
    test_clip2d_mask_varying_dt<T,S,1,MaxDt> (rng);
}

// Many Df, many Dt
template<typename T, unsigned int S, unsigned int MaxDf, unsigned int MaxDt, typename std::enable_if<(MaxDf>1),int>::type = 0>
static void test_clip2d_mask_all(std::mt19937 &rng)
{
    test_clip2d_mask_all<T,S,MaxDf/2,MaxDt> (rng);
    test_clip2d_mask_varying_dt<T,S,MaxDf,MaxDt> (rng);
}



// -------------------------------------------------------------------------------------------------


int main(int argc, char **argv)
{
    std::random_device rd;
    std::mt19937 rng(rd());

    test_clip2d_wrms_all<float,8,32,32> (rng);
    test_clip2d_mask_all<float,8,32,32> (rng);

    cout << "test-kernels: all tests passed\n";
    return 0;
}
