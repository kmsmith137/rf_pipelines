// FIXME (low-priority) a nuisance issue when working with this code is that functions
// which are very similar have different argument orderings, e.g.
//
//          make_intensity_clipper(nt_chunk, axis, sigma, niter, iter_sigma, Df, Dt)
//   calls _make_intensity_clipper(Df, Dt, axis, nt_chunk, sigma, niter, iter_sigma)

// FIXME: currently we need to compile a new kernel for every (Df,Dt) pair, where
// Df,Dt are the frequency/time downsampling factors.  Eventually I'd like to 
// improve this by having special kernels to handle the large-Df and large-Dt cases.

// #include <simd_helpers/simd_debug.hpp>

#include <array>
#include "rf_pipelines_internals.hpp"
#include "kernels/intensity_clippers.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif


// -------------------------------------------------------------------------------------------------
//
// intensity_clipper_kernel_table


struct intensity_clipper_kernels {
    // f_nds(nds_int, nds_wt, nfreq, nt)
    // f_clip(intensity, weights, nfreq, nt, stride, niter, sigma, iter_sigma, ds_int, ds_wt)

    void (*f_nds)(int &, int &, int, int);
    void (*f_clip)(float *, float *, int, int, int, int, double, double, float *, float *);
};


// Fills length-3 array indexed by axis
template<unsigned int S, unsigned int Df, unsigned int Dt, bool IterFlag>
inline void fill_1d_intensity_clipper_kernel_table(intensity_clipper_kernels *out)
{
    out[AXIS_FREQ] = { _kernel_nds_1d_f<S,Df,Dt,IterFlag>, _kernel_clip_1d_f<S,Df,Dt,IterFlag> };
    out[AXIS_TIME] = { _kernel_nds_1d_t<S,Df,Dt,IterFlag>, _kernel_clip_1d_t<S,Df,Dt,IterFlag> };
    out[AXIS_NONE] = { _kernel_nds_2d<S,Df,Dt,IterFlag>, _kernel_clip_2d<float,S,Df,Dt,IterFlag> };
}

// Fills shape-(2,3) array indexed by (iterflag, axis)
template<unsigned int S, unsigned int Df, unsigned int Dt>
inline void fill_2d_intensity_clipper_table(intensity_clipper_kernels *out)
{
    fill_1d_intensity_clipper_kernel_table<S,Df,Dt,false> (out);
    fill_1d_intensity_clipper_kernel_table<S,Df,Dt,true> (out+3);
}

// Fills shape-(NDt,2,3) array indexed by (Dt, iterflag, axis)
template<unsigned int S, unsigned int Df, unsigned int NDt, typename enable_if<(NDt==0),int>::type = 0>
inline void fill_3d_intensity_clipper_table(intensity_clipper_kernels *out) { }

template<unsigned int S, unsigned int Df, unsigned int NDt, typename enable_if<(NDt>0),int>::type = 0>
inline void fill_3d_intensity_clipper_table(intensity_clipper_kernels *out) 
{ 
    fill_3d_intensity_clipper_table<S,Df,(NDt-1)> (out);
    fill_2d_intensity_clipper_table<S,Df,(1<<(NDt-1))> (out + 6*(NDt-1));
}

// Fills shape-(NDf,NDt,2,3) array indexed by (Df, Dt, iterflag, axis)
template<unsigned int S, unsigned int NDf, unsigned int NDt, typename enable_if<(NDf==0),int>::type = 0>
inline void fill_4d_intensity_clipper_table(intensity_clipper_kernels *out) { }

template<unsigned int S, unsigned int NDf, unsigned int NDt, typename enable_if<(NDf>0),int>::type = 0>
inline void fill_4d_intensity_clipper_table(intensity_clipper_kernels *out)
{
    fill_4d_intensity_clipper_table<S,(NDf-1),NDt> (out);
    fill_3d_intensity_clipper_table<S,(1<<(NDf-1)),NDt> (out + 6*(NDf-1)*NDt);
}


struct intensity_clipper_kernel_table {
    static constexpr int S = constants::single_precision_simd_length;
    static constexpr int MaxDf = constants::max_frequency_downsampling;
    static constexpr int MaxDt = constants::max_time_downsampling;
    static constexpr int NDf = IntegerLog2<MaxDf>() + 1;
    static constexpr int NDt = IntegerLog2<MaxDt>() + 1;

    // Weird: for some reason using std::max() here gives a clang linker (not compiler) error.
    static constexpr int MaxD = (MaxDf > MaxDt) ? MaxDf : MaxDt;

    vector<intensity_clipper_kernels> kernels;

    integer_log2_lookup_table ilog2_lookup;

    intensity_clipper_kernel_table() :
	kernels(6*NDf*NDt), ilog2_lookup(MaxD)
    {
	fill_4d_intensity_clipper_table<S,NDf,NDt> (&kernels[0]);
    }

    // Caller must call check_params()!
    inline intensity_clipper_kernels get_kernel(axis_type axis, int Df, int Dt, int niter)
    {
	int idf = ilog2_lookup(Df);
	int idt = ilog2_lookup(Dt);
	int iflag = (niter > 1) ? 1 : 0;

	return kernels[6*(idf*NDt+idt) + 3*iflag + axis];
    }
};


static intensity_clipper_kernel_table global_intensity_clipper_kernel_table;


// -------------------------------------------------------------------------------------------------
//
// clipper_transform


struct clipper_transform : public wi_transform 
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
    intensity_clipper_kernels kernels;

    // Noncopyable
    clipper_transform(const clipper_transform &) = delete;
    clipper_transform &operator=(const clipper_transform &) = delete;

    clipper_transform(int nds_f_, int nds_t_, axis_type axis_, int nt_chunk_, double sigma_, int niter_, double iter_sigma_, intensity_clipper_kernels kernels_)
	: nds_f(nds_f_), nds_t(nds_t_), axis(axis_), niter(niter_), sigma(sigma_), iter_sigma(iter_sigma_ ? iter_sigma_ : sigma_), kernels(kernels_)
    {
	stringstream ss;
        ss << "intensity_clipper_cpp(nt_chunk=" << nt_chunk_ << ", axis=" << axis << ", sigma=" << sigma
           << ", niter=" << niter << ", iter_sigma=" << iter_sigma << ", Df=" << nds_f << ", Dt=" << nds_t << ")";

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

    virtual ~clipper_transform()
    {
	free(ds_intensity);
	free(ds_weights);
	ds_intensity = ds_weights = nullptr;
    }

    virtual void set_stream(const wi_stream &stream) override
    {
	rf_assert(stream.nfreq % nds_f == 0);

	int nds_int, nds_wt;
	this->kernels.f_nds(nds_int, nds_wt, stream.nfreq, nt_chunk);

	this->nfreq = stream.nfreq;
	this->ds_intensity = aligned_alloc<float> (nds_int);
	this->ds_weights = aligned_alloc<float> (nds_wt);
    }

    virtual void process_chunk(double t0, double t1, float *intensity, float *weights, ssize_t stride, float *pp_intensity, float *pp_weights, ssize_t pp_stride) override
    {
	this->kernels.f_clip(intensity, weights, nfreq, nt_chunk, stride, niter, sigma, iter_sigma, ds_intensity, ds_weights);
    }

    virtual void start_substream(int isubstream, double t0) override { }
    virtual void end_substream() override { }
};


// -------------------------------------------------------------------------------------------------


static void check_params(int Df, int Dt, axis_type axis, int nfreq, int nt, int stride, double sigma, int niter, double iter_sigma)
{
    static constexpr int S = constants::single_precision_simd_length;
    static constexpr int MaxDf = constants::max_frequency_downsampling;
    static constexpr int MaxDt = constants::max_time_downsampling;

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


// externally visible
shared_ptr<wi_transform> make_intensity_clipper(int nt_chunk, axis_type axis, double sigma, int niter, double iter_sigma, int Df, int Dt)
{
    int dummy_nfreq = Df;         // arbitrary
    int dummy_stride = nt_chunk;  // arbitrary

    check_params(Df, Dt, axis, dummy_nfreq, nt_chunk, dummy_stride, sigma, niter, iter_sigma);
    
    intensity_clipper_kernels kernels = global_intensity_clipper_kernel_table.get_kernel(axis, Df, Dt, niter);
    return make_shared<clipper_transform> (Df, Dt, axis, nt_chunk, sigma, niter, iter_sigma, kernels);
}


// externally visible
void apply_intensity_clipper(const float *intensity, float *weights, int nfreq, int nt, int stride, axis_type axis, double sigma, int niter, double iter_sigma, int Df, int Dt)
{
    check_params(Df, Dt, axis, nfreq, nt, stride, sigma, niter, iter_sigma);

    intensity_clipper_kernels kernels = global_intensity_clipper_kernel_table.get_kernel(axis, Df, Dt, niter);

    int nds_int, nds_wt;
    kernels.f_nds(nds_int, nds_wt, nfreq, nt);

    float *ds_int = aligned_alloc<float> (nds_int);
    float *ds_wt = aligned_alloc<float> (nds_wt);

    kernels.f_clip(const_cast<float *> (intensity), weights, nfreq, nt, stride, niter, sigma, iter_sigma, ds_int, ds_wt);

    free(ds_int);
    free(ds_wt);
}


// Externally visible
void weighted_mean_and_rms(float &mean, float &rms, const float *intensity, const float *weights, int nfreq, int nt, int stride, int niter, double sigma)
{
    static constexpr int S = constants::single_precision_simd_length;

    if (_unlikely((nfreq <= 0) || (nt <= 0)))
	throw runtime_error("weighted_mean_and_rms(): (nfreq,nt)=(" + to_string(nfreq) + "," + to_string(nt) + ") is invalid");
    
    if (_unlikely(stride < nt))
	throw runtime_error("weighted_mean_and_rms(): stride=" + to_string(stride) + " is < nt=" + to_string(nt));	

    if (_unlikely(sigma < 1.0))
	throw runtime_error("weighted_mean_and_rms(): sigma=" + to_string(sigma) + ", expected >= 1.0");

    if (_unlikely(niter < 1))
	throw runtime_error("weighted_mean_and_rms(): niter=" + to_string(niter) + " is invalid");

    if (_unlikely(nt % S))
	throw runtime_error("weighted_mean_and_rms(): nt=" + to_string(nt) + " must be divisible by S=" + to_string(S) + ", the single-precision simd length on this machine");

    if (_unlikely(!intensity))
	throw runtime_error("weighted_mean_and_rms(): 'intensity' argument is a NULL pointer");

    if (_unlikely(!weights))
	throw runtime_error("weighted_mean_and_rms(): 'weights' argument is a NULL pointer");

    simd_t<float,S> mean_x, rms_x;

    // Note: in the case (Df,Dt)=(1,1), the value of the compile-time argument IterFlag doesn't matter,
    // and it's always OK to take ds_int = ds_wt = NULL.
    _kernel_iterative_wrms_2d<float,S,1,1,false> (mean_x, rms_x, intensity, weights, nfreq, nt, stride, niter, sigma, NULL, NULL);

    mean = mean_x.template extract<0> ();
    rms = rms_x.template extract<0> ();
}


}  // namespace rf_pipelines
