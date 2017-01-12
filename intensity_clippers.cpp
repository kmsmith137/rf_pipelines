// #include <simd_helpers/simd_debug.hpp>

#include <array>
#include "rf_pipelines_internals.hpp"
#include "kernels/intensity_clippers.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif


// f_nds(nds_int, nds_wt, nfreq, nt)
// f_clip(intensity, weights, nfreq, nt, stride, niter, sigma, iter_sigma, ds_int, ds_wt)

using kernel_nds_t = void (*)(int &, int &, int, int);
using kernel_clip_t = void (*)(float *, float *, int, int, int, int, double, double, float *, float *);


// -------------------------------------------------------------------------------------------------
//
// clipper_transform_base
//
// The compile-time arguments Df,Dt are the frequency/time downsampling factors, and the
// compile-time boolean argument IterFlag should be set to 'true' if and only if niter > 1.
//
// FIXME: currently we need to compile a new kernel for every (Df,Dt) pair.  Eventually I'd
// like to improve this by having special kernels to handle the large-Df and large-Dt cases.

struct clipper_transform_base : public wi_transform 
{
    // (Frequency, time) downsampling factors and axis.
    const int nds_f;
    const int nds_t;
    const axis_type axis;
    
    // Clipping thresholds.
    const int niter;
    const double sigma;
    const double iter_sigma;

    // Allocated in set_stream()
    float *ds_intensity = nullptr;
    float *ds_weights = nullptr;

    // Kernels
    kernel_nds_t f_nds;
    kernel_clip_t f_clip;

    // Noncopyable
    clipper_transform_base(const clipper_transform_base &) = delete;
    clipper_transform_base &operator=(const clipper_transform_base &) = delete;

    clipper_transform_base(int nds_f_, int nds_t_, axis_type axis_, int nt_chunk_, double sigma_, int niter_, double iter_sigma_, kernel_nds_t f_nds_, kernel_clip_t f_clip_)
	: nds_f(nds_f_), nds_t(nds_t_), axis(axis_), niter(niter_), sigma(sigma_), iter_sigma(iter_sigma_ ? iter_sigma_ : sigma_), f_nds(f_nds_), f_clip(f_clip_)
    {
	stringstream ss;
	ss << "intensity_clipper_transform(Df=" << nds_f << ",Dt=" << nds_t << ",axis=" << axis
	   << ",nt_chunk=" << nt_chunk_ << ",sigma=" << sigma << ",niter=" << niter 
	   << ",iter_sigma=" << iter_sigma << ")";

	this->name = ss.str();
	this->nt_chunk = nt_chunk_;
	this->nt_prepad = 0;
	this->nt_postpad = 0;

	// No need to make these asserts "verbose", since they should have been checked in make_intensity_clipper().
	rf_assert(sigma >= 1.0);
	rf_assert(iter_sigma >= 1.0);
	rf_assert(nt_chunk > 0);
	rf_assert(nt_chunk % nds_t == 0);
    }

    virtual ~clipper_transform_base()
    {
	free(ds_intensity);
	free(ds_weights);
	ds_intensity = ds_weights = nullptr;
    }

    virtual void set_stream(const wi_stream &stream) override
    {
	rf_assert(stream.nfreq % nds_f == 0);

	int nds_int, nds_wt;
	this->f_nds(nds_int, nds_wt, nfreq, nt_chunk);

	this->nfreq = stream.nfreq;
	this->ds_intensity = aligned_alloc<float> (nds_int);
	this->ds_weights = aligned_alloc<float> (nds_wt);
    }

    virtual void process_chunk(double t0, double t1, float *intensity, float *weights, ssize_t stride, float *pp_intensity, float *pp_weights, ssize_t pp_stride) override
    {
	f_clip(intensity, weights, nfreq, nt_chunk, stride, niter, sigma, iter_sigma, ds_intensity, ds_weights);
    }

    virtual void start_substream(int isubstream, double t0) override { }
    virtual void end_substream() override { }
};


// -------------------------------------------------------------------------------------------------
//
// clipper_transform_2d


template<unsigned int S, unsigned int Df, unsigned int Dt, bool IterFlag>
static void kernel_nds_2d(int &nds_int, int &nds_wt, int nfreq, int nt)
{
    static constexpr bool DsiFlag = (Df > 1) || (Dt > 1);
    static constexpr bool DswFlag = IterFlag && ((Df > 1) || (Dt > 1));

    nds_int = DsiFlag ? ((nfreq/Df) * (nt/Dt)) : 0;
    nds_wt = DswFlag ? ((nfreq/Df) * (nt/Dt)) : 0;
}


template<unsigned int S, unsigned int Df, unsigned int Dt, bool IterFlag>
static void kernel_clip_2d(float *intensity, float *weights, int nfreq, int nt, int stride, int niter, double sigma, double iter_sigma, float *ds_int, float *ds_wt)
{
    static constexpr bool DsiFlag = (Df > 1) || (Dt > 1);
    static constexpr bool DswFlag = IterFlag && ((Df > 1) || (Dt > 1));

    simd_t<float,S> mean, rms;
    _kernel_clip2d_wrms<float,S,Df,Dt,DsiFlag,DswFlag,float,S> (mean, rms, intensity, weights, nfreq, nt, stride, ds_int, ds_wt);

    const float *s_intensity = DsiFlag ? ds_int : intensity;
    const float *s_weights = DswFlag ? ds_wt : weights;
    int s_stride = DsiFlag ? (nt/Dt) : stride;   // must use DsiFlag here, not DswFlag
	
    for (int iter = 1; iter < niter; iter++) {
	// (s_intensity, s_weights, iter_sigma) here
	simd_t<float,S> thresh = simd_t<float,S>(iter_sigma) * rms;
	_kernel_clip2d_iterate<float,S> (mean, rms, s_intensity, s_weights, mean, thresh, nfreq/Df, nt/Dt, s_stride);
    }

    // needs simd_debug.hpp
    // cerr << "kernel_clip_2d: mean=" << mean << endl;
    // cerr << "kernel_clip_2d: rms=" << rms << endl;

    // (s_intensity, weights, sigma) here
    simd_t<float,S> thresh = simd_t<float,S>(sigma) * rms;
    _kernel_clip2d_mask<float,S,Df,Dt> (weights, s_intensity, mean, thresh, nfreq, nt, stride, s_stride);
}


template<unsigned int S, unsigned int Df, unsigned int Dt, bool IterFlag>
struct clipper_transform_2d : public clipper_transform_base
{
    clipper_transform_2d(int nt_chunk_, double sigma_, int niter_, double iter_sigma_)
	: clipper_transform_base(Df, Dt, AXIS_NONE, nt_chunk_, sigma_, niter_, iter_sigma_, kernel_nds_2d<S,Df,Dt,IterFlag>, kernel_clip_2d<S,Df,Dt,IterFlag>)
    { 
	if (IterFlag) rf_assert(niter > 1);
	if (!IterFlag) rf_assert(niter == 1);
    }
};


// -------------------------------------------------------------------------------------------------
//
// clipper_transform_time_axis
//
// Currently implemented by calling the 2d kernels many times with nfreq=1.
//
// FIXME there is a little extra overhead here, should improve by writing real 1D kernels.


template<unsigned int S, unsigned int Df, unsigned int Dt, bool IterFlag>
static void kernel_nds_1d_t(int &nds_int, int &nds_wt, int nfreq, int nt)
{
    static constexpr bool DsiFlag = (Df > 1) || (Dt > 1);
    static constexpr bool DswFlag = IterFlag && ((Df > 1) || (Dt > 1));

    nds_int = DsiFlag ? (nt/Dt) : 0;
    nds_wt = DswFlag ? (nt/Dt) : 0;
}


template<unsigned int S, unsigned int Df, unsigned int Dt, bool IterFlag>
static void kernel_clip_1d_t(float *intensity, float *weights, int nfreq, int nt, int stride, int niter, double sigma, double iter_sigma, float *ds_int, float *ds_wt)
{
    static constexpr bool DsiFlag = (Df > 1) || (Dt > 1);
    static constexpr bool DswFlag = IterFlag && ((Df > 1) || (Dt > 1));

    for (int ifreq = 0; ifreq < nfreq; ifreq += Df) {
	float *irow = intensity + ifreq * stride;
	float *wrow = weights + ifreq * stride;
	
	// We pass nfreq=Df to _kernel_clip2d_wrms, not the "true" nfreq
	simd_t<float,S> mean, rms;
	_kernel_clip2d_wrms<float,S,Df,Dt,DsiFlag,DswFlag,float,S> (mean, rms, irow, wrow, Df, nt, stride, ds_int, ds_wt);
	
	const float *irow2 = DsiFlag ? ds_int : irow;
	const float *wrow2 = DswFlag ? ds_wt : wrow;
	
	for (int iter = 1; iter < niter; iter++) {
	    // Here we pass nfreq=1 and stride=0
	    simd_t<float,S> thresh = simd_t<float,S>(iter_sigma) * rms;
	    _kernel_clip2d_iterate<float,S> (mean, rms, irow2, wrow2, mean, thresh, 1, nt/Dt, 0);    // (irow2, wrow2, iter_sigma)
	}

	// Here we pass nfreq=Df.  Setting both strides to 'stride' is OK but this isn't completely obvious.
	simd_t<float,S> thresh = simd_t<float,S>(sigma) * rms;
	_kernel_clip2d_mask<float,S,Df,Dt> (wrow, irow2, mean, thresh, Df, nt, stride, stride);    // (wrow, irow2, sigma)
    }
}


template<unsigned int S, unsigned int Df, unsigned int Dt, bool IterFlag>
struct clipper_transform_time_axis : clipper_transform_base
{
    clipper_transform_time_axis(int nt_chunk_, double sigma_, int niter_, double iter_sigma_) :
	clipper_transform_base(Df, Dt, AXIS_TIME, nt_chunk_, sigma_, niter_, iter_sigma_, kernel_nds_1d_t<S,Df,Dt,IterFlag>, kernel_clip_1d_t<S,Df,Dt,IterFlag>)
    { 
	if (IterFlag) rf_assert(niter > 1);
	if (!IterFlag) rf_assert(niter == 1);
    }
};


// -------------------------------------------------------------------------------------------------
//
// clipper_transform_freq_axis
//
// Currently implemented by calling the 2d kernels many times with nfreq=1.
//
// FIXME there is a little extra overhead here, should improve by writing real 1D kernels.


template<unsigned int S, unsigned int Df, unsigned int Dt, bool IterFlag>
static void kernel_nds_1d_f(int &nds_int, int &nds_wt, int nfreq, int nt)
{
    static constexpr bool DsiFlag = (Df > 1) || (Dt > 1);
    static constexpr bool DswFlag = IterFlag && ((Df > 1) || (Dt > 1));

    nds_int = DsiFlag ? ((nfreq/Df) * S) : 0;
    nds_wt = DswFlag ? ((nfreq/Df) * S) : 0;
}


template<unsigned int S, unsigned int Df, unsigned int Dt, bool IterFlag>
static void kernel_clip_1d_f(float *intensity, float *weights, int nfreq, int nt, int stride, int niter, double sigma, double iter_sigma, float *ds_int, float *ds_wt)
{
    static constexpr bool DsiFlag = (Df > 1) || (Dt > 1);
    static constexpr bool DswFlag = IterFlag && ((Df > 1) || (Dt > 1));

    for (int it = 0; it < nt; it += Dt*S) {
	float *icol = intensity + it;
	float *wcol = weights + it;
	
	simd_t<float,S> mean, rms;	
	_kernel_clip1d_f_wrms<float,S,Df,Dt,DsiFlag,DswFlag> (mean, rms, icol, wcol, nfreq, stride, ds_int, ds_wt);
	
	const float *icol2 = DsiFlag ? ds_int : icol;
	const float *wcol2 = DswFlag ? ds_wt : wcol;
	int stride2 = DsiFlag ? S : stride;   // must use DsiFlag here, not DswFlag
	
	for (int iter = 1; iter < niter; iter++) {
	    // Here we can use _kernel_clip2d_iterate() with (nfreq,nt) replaced by (nfreq/Df,S)
	    simd_t<float,S> thresh = simd_t<float,S>(iter_sigma) * rms;
	    _kernel_clip2d_iterate<float,S> (mean, rms, icol2, wcol2, mean, thresh, nfreq/Df, S, stride2);    // (irow2, wrow2, iter_sigma)
	}

	simd_t<float,S> thresh = simd_t<float,S>(sigma) * rms;
	_kernel_clip1d_f_mask<float,S,Df,Dt> (wcol, icol2, mean, thresh, nfreq, stride, stride2);
    }
}


template<unsigned int S, unsigned int Df, unsigned int Dt, bool IterFlag>
struct clipper_transform_freq_axis : clipper_transform_base
{
    clipper_transform_freq_axis(int nt_chunk_, double sigma_, int niter_, double iter_sigma_) :
	clipper_transform_base(Df, Dt, AXIS_FREQ, nt_chunk_, sigma_, niter_, iter_sigma_, kernel_nds_1d_f<S,Df,Dt,IterFlag>, kernel_clip_1d_f<S,Df,Dt,IterFlag>)
    { 
	if (IterFlag) rf_assert(niter > 1);
	if (!IterFlag) rf_assert(niter == 1);
    }
};


// -------------------------------------------------------------------------------------------------
//
// Boilerplate needed to instantiate templates and export factory functions.


template<unsigned int S, unsigned int Df, unsigned int Dt, unsigned int IterFlag>
inline shared_ptr<clipper_transform_base> _make_intensity_clipper5(axis_type axis, int nt_chunk, double sigma, int niter, double iter_sigma)
{
    if (axis == AXIS_FREQ)
	return make_shared<clipper_transform_freq_axis<S,Df,Dt,IterFlag>> (nt_chunk, sigma, niter, iter_sigma);
    if (axis == AXIS_TIME)
	return make_shared<clipper_transform_time_axis<S,Df,Dt,IterFlag>> (nt_chunk, sigma, niter, iter_sigma);
    if (axis == AXIS_NONE)
	return make_shared<clipper_transform_2d<S,Df,Dt,IterFlag>> (nt_chunk, sigma, niter, iter_sigma);

    throw runtime_error("rf_pipelines::make_intensity_clipper(): axis='" + to_string(axis) + "' is not a valid value");
}


template<unsigned int S, unsigned int Df, unsigned int Dt>
inline shared_ptr<clipper_transform_base> _make_intensity_clipper4(axis_type axis, int nt_chunk, double sigma, int niter, double iter_sigma)
{
    if (niter > 1)
	return _make_intensity_clipper5<S,Df,Dt,true> (axis, nt_chunk, sigma, niter, iter_sigma);
    else
	return _make_intensity_clipper5<S,Df,Dt,false> (axis, nt_chunk, sigma, niter, iter_sigma);
}


template<unsigned int S, unsigned int Df, unsigned int MaxDt, typename std::enable_if<(MaxDt==0),int>::type = 0>
inline shared_ptr<clipper_transform_base> _make_intensity_clipper3(int Dt, axis_type axis, int nt_chunk, double sigma, int niter, double iter_sigma)
{
    throw runtime_error("rf_pipelines internal error: Dt=" + to_string(Dt) + " not found in intensity_clipper template chain");
}

template<unsigned int S, unsigned int Df, unsigned int MaxDt, typename std::enable_if<(MaxDt > 0),int>::type = 0>
inline shared_ptr<clipper_transform_base> _make_intensity_clipper3(int Dt, axis_type axis, int nt_chunk, double sigma, int niter, double iter_sigma)
{
    if (Dt == MaxDt)
	return _make_intensity_clipper4<S,Df,MaxDt> (axis, nt_chunk, sigma, niter, iter_sigma);
    else
	return _make_intensity_clipper3<S,Df,(MaxDt/2)> (Dt, axis, nt_chunk, sigma, niter, iter_sigma);
}


template<unsigned int S, unsigned int MaxDf, unsigned int MaxDt, typename std::enable_if<(MaxDf==0),int>::type = 0>
inline shared_ptr<clipper_transform_base> _make_intensity_clipper2(int Df, int Dt, axis_type axis, int nt_chunk, double sigma, int niter, double iter_sigma)
{
    throw runtime_error("rf_pipelines internal error: Df=" + to_string(Df) + " not found in intensity_clipper template chain");
}

template<unsigned int S, unsigned int MaxDf, unsigned int MaxDt, typename std::enable_if<(MaxDf > 0),int>::type = 0>
inline shared_ptr<clipper_transform_base> _make_intensity_clipper2(int Df, int Dt, axis_type axis, int nt_chunk, double sigma, int niter, double iter_sigma)
{
    if (Df == MaxDf)
	return _make_intensity_clipper3<S,MaxDf,MaxDt> (Dt, axis, nt_chunk, sigma, niter, iter_sigma);
    else
	return _make_intensity_clipper2<S,(MaxDf/2),MaxDt> (Df, Dt, axis, nt_chunk, sigma, niter, iter_sigma);
}


static void check_params(int Df, int Dt, axis_type axis, int nfreq, int nt, int stride, double sigma, int niter, double iter_sigma)
{
    static constexpr int S = constants::single_precision_simd_length;
    static constexpr int MaxDf = constants::intensity_clipper_max_frequency_downsampling;
    static constexpr int MaxDt = constants::intensity_clipper_max_time_downsampling;

    if (_unlikely((Df <= 0) || !is_power_of_two(Df)))
	throw runtime_error("rf_pipelines intensity clipper: Df=" + to_string(Df) + " must be a power of two");

    if (_unlikely((Dt <= 0) || !is_power_of_two(Dt)))
	throw runtime_error("rf_pipelines intensity clipper: Dt=" + to_string(Dt) + " must be a power of two");

    if (_unlikely((axis != AXIS_FREQ) && (axis != AXIS_TIME) && (axis != AXIS_NONE)))
	throw runtime_error("rf_pipelines intensity clipper: axis=" + stringify(axis) + " is not defined for this transform");

    if (_unlikely(nfreq <= 0))
	throw runtime_error("rf_pipelines intensity clipper: nfreq=" + to_string(nfreq) + ", positive value was expected");

    if (_unlikely(nt <= 0))
	throw runtime_error("rf_pipelines intensity clipper: nt=" + to_string(nt) + ", positive value was expected");

    if (_unlikely(abs(stride) < nt))
	throw runtime_error("rf_pipelines intensity clipper: stride=" + to_string(stride) + " must be >= nt");

    if (_unlikely(sigma < 1.0))
	throw runtime_error("rf_pipelines intensity clipper: sigma=" + to_string(sigma) + " must be >= 1.0");

    if (_unlikely(niter < 1))
	throw runtime_error("rf_pipelines intensity clipper: niter=" + to_string(niter) + " must be >= 1");

    if (_unlikely((nfreq % Df) != 0))
	throw runtime_error("rf_pipelines intensity clipper: nfreq=" + to_string(nfreq)
			    + " must be a multiple of the downsampling factor Df=" + to_string(Df));
    
    if (_unlikely((nt % (Dt*S)) != 0))
	throw runtime_error("rf_pipelines intensity clipper: nt=" + to_string(nt)
			    + " must be a multiple of the downsampling factor Dt=" + to_string(Dt)
			    + " multiplied by constants::single_precision_simd_length=" + to_string(S));

    if (_unlikely((Df > MaxDf) || (Dt > MaxDt)))
	throw runtime_error("rf_pipelines intensity clipper: (Df,Dt)=(" + to_string(Df) + "," + to_string(Dt) + ")"
			    + " exceeds compile time limits; to fix this see 'constants' in rf_pipelines.hpp");
}


// Non-static
// Caller must call check_params() first.
shared_ptr<clipper_transform_base> _make_intensity_clipper(int Df, int Dt, axis_type axis, int nt_chunk, double sigma, int niter, double iter_sigma)
{
    static constexpr int S = constants::single_precision_simd_length;
    static constexpr int MaxDf = constants::intensity_clipper_max_frequency_downsampling;
    static constexpr int MaxDt = constants::intensity_clipper_max_time_downsampling;

    shared_ptr<clipper_transform_base> ret = _make_intensity_clipper2<S,MaxDf,MaxDt> (Df, Dt, axis, nt_chunk, sigma, niter, iter_sigma);
    
    // Sanity check the template instantiation
    assert(ret->nds_f == Df);
    assert(ret->nds_t == Dt);
    assert(ret->axis == axis);
    assert(ret->niter == niter);
    assert(ret->sigma == sigma);
    assert(ret->iter_sigma == iter_sigma);
    
    return ret;
}


// Externally callable "bottom line" routine which returns clipper transform with specified downsampling and axis.
shared_ptr<wi_transform> make_intensity_clipper(int Df, int Dt, axis_type axis, int nt_chunk, double sigma, int niter, double iter_sigma)
{
    int dummy_nfreq = Df;         // arbitrary
    int dummy_stride = nt_chunk;  // arbitrary

    check_params(Df, Dt, axis, dummy_nfreq, nt_chunk, dummy_stride, sigma, niter, iter_sigma);
    return _make_intensity_clipper(Dt, Dt, axis, nt_chunk, sigma, niter, iter_sigma);
}


// -------------------------------------------------------------------------------------------------


// Compile-time integer-valued log_2()
template<unsigned int D, typename std::enable_if<(D==1),int>::type = 0>
inline constexpr int IntegerLog2() { return 0; }

template<unsigned int D, typename std::enable_if<(D>1 && (D%2)==0),int>::type = 0>
inline constexpr int IntegerLog2() { return IntegerLog2<(D/2)>() + 1; }


struct clipper_table {
    static constexpr int MaxDf = constants::intensity_clipper_max_frequency_downsampling;
    static constexpr int MaxDt = constants::intensity_clipper_max_time_downsampling;

    static constexpr int NDf = IntegerLog2<MaxDf>() + 1;
    static constexpr int NDt = IntegerLog2<MaxDt>() + 1;

    struct kernels {
	kernel_nds_t f_nds;
	kernel_clip_t f_clip;
    };

    using ktab4_t = std::array<kernels, 2>;
    using ktab3_t = std::array<ktab4_t, NDt>;
    using ktab2_t = std::array<ktab3_t, NDf>;
    using ktab_t = std::array<ktab2_t, 3>;

    ktab_t entries;

    // Lookup table for integer-valued log_2
    std::vector<int> ilog2_lookup;

    clipper_table()
    {
	static constexpr int S = constants::single_precision_simd_length;

	for (axis_type axis: { AXIS_FREQ, AXIS_TIME, AXIS_NONE }) {
	    for (int idf = 0; idf < NDf; idf++) {
		for (int idt = 0; idt < NDt; idt++) {
		    for (int niter = 1; niter < 3; niter++) {
			int Df = 1 << idf;
			int Dt = 1 << idt;
			int iflag = (niter > 1) ? 1 : 0;

			auto t = _make_intensity_clipper(Df, Dt, axis, Dt*S, 2.0, niter, 2.0);   // nt_chunk, sigma arbitrary
			entries.at(axis).at(idf).at(idt).at(iflag) = { t->f_nds, t->f_clip };
		    }
		}
	    }
	}

	// Weird: for some reason using std::max() here gives a clang linker (not compiler) error.
	int n = (MaxDf > MaxDt) ? (MaxDf+1) : (MaxDt+1);

	ilog2_lookup = vector<int> (n, -1);

	for (int i = 0; (1<<i) < n; i++)
	    ilog2_lookup[1<<i] = i;
    }

    // Caller must call check_params()!
    inline kernels get(axis_type axis, int Df, int Dt, int niter)
    {
	int idf = ilog2_lookup[Df];
	int idt = ilog2_lookup[Dt];
	int iflag = (niter > 1) ? 1 : 0;
	
	return entries.at(axis).at(idf).at(idt).at(iflag);
    }
};


void apply_intensity_clipper(const float *intensity, float *weights, int nfreq, int nt, int stride, axis_type axis, double sigma, int niter, double iter_sigma, int Df, int Dt)
{
    static clipper_table ktab;

    check_params(Df, Dt, axis, nfreq, nt, stride, sigma, niter, iter_sigma);

    clipper_table::kernels k = ktab.get(axis, Df, Dt, niter);

    int nds_int, nds_wt;
    k.f_nds(nds_int, nds_wt, nfreq, nt);

    float *ds_int = aligned_alloc<float> (nds_int);
    float *ds_wt = aligned_alloc<float> (nds_wt);

    k.f_clip(const_cast<float *> (intensity), weights, nfreq, nt, stride, niter, sigma, iter_sigma, ds_int, ds_wt);

    free(ds_int);
    free(ds_wt);
}


}  // namespace rf_pipelines

