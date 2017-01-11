// FIXME: currently we need to compile a new kernel for every (Df,Dt) pair, where
// Df,Dt are the frequency/time downsampling factors.  Eventually I'd like to 
// improve this by having special kernels to handle the large-Df and large-Dt cases.

#include <array>
#include <cassert>

#include "rf_pipelines_internals.hpp"
#include "kernels/std_dev_clippers.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif


using mask_t = simd_helpers::smask_t<float>;

// f_ntmp(nfreq, nt)
// f_clip(intensity, weights, nfreq, nt, stride, sigma, tmp_sd, tmp_valid)

using kernel_ntmp_t = int (*) (int, int);
using kernel_clip_t = void (*) (float *, float *, int, int, int, double, float *, mask_t *);


// -------------------------------------------------------------------------------------------------
//
// sd_clipper_transform_base


struct sd_clipper_transform_base : public wi_transform 
{
    // (Frequency, time) downsampling factors and axis.
    const int nds_f;
    const int nds_t;
    const axis_type axis;
    
    // Clipping threshold.
    const double sigma;

    kernel_ntmp_t f_ntmp;
    kernel_clip_t f_clip;

    // Allocated in set_stream()
    float *tmp_sd = nullptr;
    mask_t *tmp_valid = nullptr;

    // Noncopyable
    sd_clipper_transform_base(const sd_clipper_transform_base &) = delete;
    sd_clipper_transform_base &operator=(const sd_clipper_transform_base &) = delete;

    sd_clipper_transform_base(int nds_f_, int nds_t_, axis_type axis_, int nt_chunk_, double sigma_, kernel_ntmp_t f_ntmp_, kernel_clip_t f_clip_)
	: nds_f(nds_f_), nds_t(nds_t_), axis(axis_), sigma(sigma_), f_ntmp(f_ntmp_), f_clip(f_clip_)
    {
	stringstream ss;
	ss << "std_dev_clipper_transform(Df=" << nds_f << ",Dt=" << nds_t << ",axis=" << axis
	   << ",nt_chunk=" << nt_chunk_ << ",sigma=" << sigma << ")";

	this->name = ss.str();
	this->nt_chunk = nt_chunk_;
	this->nt_prepad = 0;
	this->nt_postpad = 0;

	// No need to make these asserts "verbose", since they should have been checked in make_intensity_clipper2d().
	rf_assert(sigma >= 1.0);
	rf_assert(nt_chunk > 0);
	rf_assert(nt_chunk % nds_t == 0);
    }

    virtual ~sd_clipper_transform_base()
    {
	free(tmp_sd);
	free(tmp_valid);

	tmp_sd = nullptr;
	tmp_valid = nullptr;
    }

    virtual void set_stream(const wi_stream &stream) override
    {
	rf_assert(stream.nfreq % nds_f == 0);

	int ntmp = f_ntmp(nfreq, nt_chunk);

	this->nfreq = stream.nfreq;
	this->tmp_sd = aligned_alloc<float> (ntmp);
	this->tmp_valid = aligned_alloc<mask_t> (ntmp);
    }

    virtual void process_chunk(double t0, double t1, float *intensity, float *weights, ssize_t stride, float *pp_intensity, float *pp_weights, ssize_t pp_stride) override
    {
	f_clip(intensity, weights, nfreq, nt_chunk, stride, sigma, tmp_sd, tmp_valid);
    }

    virtual void start_substream(int isubstream, double t0) override { }
    virtual void end_substream() override { }
};


// -------------------------------------------------------------------------------------------------


static void clip_1d(int n, float *tmp_sd, mask_t *tmp_valid, double sigma)
{
    float acc0 = 0.0;
    float acc1 = 0.0;
    
    for (int i = 0; i < n; i++) {
	if (tmp_valid[i]) {
	    acc0 += 1.0;
	    acc1 += tmp_sd[i];
	}
    }
    
    if (acc0 < 1.5) {
	memset(tmp_valid, 0, n * sizeof(mask_t));
	return;
    }
    
    float mean = acc1 / acc0;
    float acc2 = 0.0;
    
    for (int i = 0; i < n; i++)
	if (tmp_valid[i])
	    acc2 += square(tmp_sd[i] - mean);
    
    float thresh = sigma * sqrtf(acc2/(acc0-1.0));

    for (int i = 0; i < n; i++) {
	if (fabs(tmp_sd[i] - mean) >= thresh)
	    tmp_valid[i] = 0;
    }
}


template<unsigned int Df> static int _kernel_ntmp_time_axis(int nfreq, int nt) { return nfreq/Df; }
template<unsigned int Dt> static int _kernel_ntmp_freq_axis(int nfreq, int nt) { return nt/Dt; }


template<unsigned int S, unsigned int Df, unsigned int Dt>
static void _kernel_clip_time_axis(float *intensity, float *weights, int nfreq, int nt, int stride, double sigma, float *tmp_sd, mask_t *tmp_valid)
{
    _kernel_std_dev_t<float,S,Df,Dt> (tmp_sd, tmp_valid, intensity, weights, nfreq, nt, stride);
    
    clip_1d(nfreq/Df, tmp_sd, tmp_valid, sigma);

    for (int ifreq = 0; ifreq < nfreq; ifreq++) {
	if (!tmp_valid[ifreq])
	    memset(weights + ifreq*stride, 0, nt * sizeof(float));
    }
}


template<unsigned int S, unsigned int Df, unsigned int Dt>
struct sd_clipper_transform_time_axis : public sd_clipper_transform_base
{
    sd_clipper_transform_time_axis(int nt_chunk_, double sigma_)
	: sd_clipper_transform_base(Df, Dt, AXIS_TIME, nt_chunk_, sigma_, _kernel_ntmp_time_axis<Df>, _kernel_clip_time_axis<S,Df,Dt>)
    { }
};


// -------------------------------------------------------------------------------------------------
//
// Boilerplate needed to instantiate templates and export factory functions.


template<unsigned int S, unsigned int Df, unsigned int Dt>
inline shared_ptr<sd_clipper_transform_base> _make_std_dev_clipper4(axis_type axis, int nt_chunk, double sigma)
{
    if (axis == AXIS_FREQ)
	throw runtime_error("sd_clipper_transform_freq_axis not written yet");
    if (axis == AXIS_TIME)
	return make_shared<sd_clipper_transform_time_axis<S,Df,Dt>> (nt_chunk, sigma);
    if (axis == AXIS_NONE)
	throw runtime_error("rf_pipelines::make_std_dev_clipper(): axis=None is not supported for this transform");

    throw runtime_error("rf_pipelines::make_std_dev_clipper(): axis='" + to_string(axis) + "' is not a valid value");
}


template<unsigned int S, unsigned int Df, unsigned int MaxDt, typename std::enable_if<(MaxDt==0),int>::type = 0>
inline shared_ptr<sd_clipper_transform_base> _make_std_dev_clipper3(int Dt, axis_type axis, int nt_chunk, double sigma)
{
    throw runtime_error("rf_pipelines internal error: Dt=" + to_string(Dt) + " not found in std_dev_clipper template chain");
}

template<unsigned int S, unsigned int Df, unsigned int MaxDt, typename std::enable_if<(MaxDt > 0),int>::type = 0>
inline shared_ptr<sd_clipper_transform_base> _make_std_dev_clipper3(int Dt, axis_type axis, int nt_chunk, double sigma)
{
    if (Dt == MaxDt)
	return _make_std_dev_clipper4<S,Df,MaxDt> (axis, nt_chunk, sigma);
    else
	return _make_std_dev_clipper3<S,Df,(MaxDt/2)> (Dt, axis, nt_chunk, sigma);
}


template<unsigned int S, unsigned int MaxDf, unsigned int MaxDt, typename std::enable_if<(MaxDf==0),int>::type = 0>
inline shared_ptr<sd_clipper_transform_base> _make_std_dev_clipper2(int Df, int Dt, axis_type axis, int nt_chunk, double sigma)
{
    throw runtime_error("rf_pipelines internal error: Df=" + to_string(Df) + " not found in std_dev_clipper template chain");
}

template<unsigned int S, unsigned int MaxDf, unsigned int MaxDt, typename std::enable_if<(MaxDf > 0),int>::type = 0>
inline shared_ptr<sd_clipper_transform_base> _make_std_dev_clipper2(int Df, int Dt, axis_type axis, int nt_chunk, double sigma)
{
    if (Df == MaxDf)
	return _make_std_dev_clipper3<S,MaxDf,MaxDt> (Dt, axis, nt_chunk, sigma);
    else
	return _make_std_dev_clipper2<S,(MaxDf/2),MaxDt> (Df, Dt, axis, nt_chunk, sigma);
}


// non-static
shared_ptr<sd_clipper_transform_base> _make_std_dev_clipper(int Df, int Dt, axis_type axis, int nt_chunk, double sigma)
{
    static constexpr int S = constants::single_precision_simd_length;
    static constexpr int MaxDf = constants::std_dev_clipper_max_frequency_downsampling;
    static constexpr int MaxDt = constants::std_dev_clipper_max_time_downsampling;
    
    auto ret = _make_std_dev_clipper2<S,MaxDf,MaxDt> (Df, Dt, axis, nt_chunk, sigma);

    // Sanity check on the template instantiation
    assert(ret->nds_f == Df);
    assert(ret->nds_t == Dt);
    assert(ret->axis == axis);
    
    return ret;
}


static void check_params(int Df, int Dt, axis_type axis, int nfreq, int nt, int stride, double sigma)
{
    static constexpr int S = constants::single_precision_simd_length;
    static constexpr int MaxDf = constants::std_dev_clipper_max_frequency_downsampling;
    static constexpr int MaxDt = constants::std_dev_clipper_max_time_downsampling;

    if (_unlikely((Df <= 0) || !is_power_of_two(Df)))
	throw runtime_error("rf_pipelines std_dev clipper: Df=" + to_string(Df) + " must be a power of two");

    if (_unlikely((Dt <= 0) || !is_power_of_two(Dt)))
	throw runtime_error("rf_pipelines std_dev clipper: Dt=" + to_string(Dt) + " must be a power of two");

    if (_unlikely((axis != AXIS_FREQ) && (axis != AXIS_TIME)))
	throw runtime_error("rf_pipelines std_dev clipper: axis=" + stringify(axis) + " is not defined for this transform");

    if (_unlikely(nfreq <= 0))
	throw runtime_error("rf_pipelines std_dev clipper: nfreq=" + to_string(nfreq) + ", positive value was expected");

    if (_unlikely(nt <= 0))
	throw runtime_error("rf_pipelines std_dev clipper: nt=" + to_string(nt) + ", positive value was expected");

    if (_unlikely(abs(stride) < nt))
	throw runtime_error("rf_pipelines std_dev clipper: stride=" + to_string(stride) + " must be >= nt");

    if (_unlikely(sigma < 2.0))
	throw runtime_error("rf_pipelines std_dev clipper: sigma=" + to_string(sigma) + " must be >= 2.0");
    
    if (_unlikely((nt % (Dt*S)) != 0))
	throw runtime_error("rf_pipelines std_dev clipper: nt=" + to_string(nt)
			    + " must be a multiple of the downsampling factor Dt=" + to_string(Dt)
			    + " multiplied by constants::single_precision_simd_length=" + to_string(S));

    if (_unlikely((Df > MaxDf) || (Dt > MaxDt)))
	throw runtime_error("rf_pipelines std_dev clipper: (Df,Dt)=(" + to_string(Df) + "," + to_string(Dt) + ")"
			    + " exceeds compile time limits; to fix this see 'constants' in rf_pipelines.hpp");
}


// Externally callable "bottom line" routine which returns clipper transform with specified downsampling and axis.

shared_ptr<wi_transform> make_std_dev_clipper(int Df, int Dt, axis_type axis, int nt_chunk, double sigma)
{
    int dummy_nfreq = 16;         // arbitrary
    int dummy_stride = nt_chunk;  // arbitrary

    check_params(Df, Dt, axis, dummy_nfreq, nt_chunk, dummy_stride, sigma);
    return _make_std_dev_clipper(Df, Dt, axis, nt_chunk, sigma);
}


// -------------------------------------------------------------------------------------------------


// Compile-time integer-valued log_2()
template<unsigned int D, typename std::enable_if<(D==1),int>::type = 0>
inline constexpr int IntegerLog2() { return 0; }

template<unsigned int D, typename std::enable_if<(D>1 && (D%2)==0),int>::type = 0>
inline constexpr int IntegerLog2() { return IntegerLog2<(D/2)>() + 1; }


struct sd_clipper_table {
    static constexpr int MaxDf = constants::std_dev_clipper_max_frequency_downsampling;
    static constexpr int MaxDt = constants::std_dev_clipper_max_time_downsampling;

    static constexpr int NDf = IntegerLog2<MaxDf>() + 1;
    static constexpr int NDt = IntegerLog2<MaxDt>() + 1;

    struct kernels {
	kernel_ntmp_t f_ntmp;
	kernel_clip_t f_clip;
    };

    using ktab3_t = std::array<kernels, (NDt)>;
    using ktab2_t = std::array<ktab3_t, (NDf)>;
    using ktab_t = std::array<ktab2_t, 2>;

    ktab_t entries;

    // Lookup table for integer-valued log_2
    std::vector<int> ilog2_lookup;

    sd_clipper_table()
    {
	for (axis_type axis: { AXIS_FREQ, AXIS_TIME }) {
	    for (int idf = 0; idf < NDf; idf++) {
		for (int idt = 0; idt < NDt; idt++) {
		    auto t = _make_std_dev_clipper(1 << idf, 1 << idt, axis, 16, 2.0);   // nt_chunk, sigma arbitrary
		    entries.at(axis).at(idf).at(idt) = { t->f_ntmp, t->f_clip };
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
    inline kernels get(axis_type axis, int Df, int Dt)
    {
	int idf = ilog2_lookup[Df];
	int idt = ilog2_lookup[Dt];
	
	return entries.at(axis).at(idf).at(idt);
    }
};


void apply_std_dev_clipper(const float *intensity, float *weights, int Df, int Dt, axis_type axis, int nfreq, int nt, int stride, double sigma)
{
    static sd_clipper_table ktab;

    check_params(Df, Dt, axis, nfreq, nt, stride, sigma);

    sd_clipper_table::kernels k = ktab.get(axis, Df, Dt);

    int ntmp = k.f_ntmp(nfreq, nt);
    float *tmp_sd = aligned_alloc<float> (ntmp);
    mask_t *tmp_valid = aligned_alloc<mask_t> (ntmp);

    k.f_clip(const_cast<float *> (intensity), weights, nfreq, nt, stride, sigma, tmp_sd, tmp_valid);

    free(tmp_sd);
    free(tmp_valid);
}


}  // namespace rf_pipelines
