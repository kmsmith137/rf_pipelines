#include "rf_pipelines_internals.hpp"
#include "kernels/intensity_clippers.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif


// -------------------------------------------------------------------------------------------------
//
// sd_clipper_transform_base
//
// The compile-time arguments Df,Dt are the frequency/time downsampling factors.
//
// FIXME: currently we need to compile a new kernel for every (Df,Dt) pair.  Eventually I'd
// like to improve this by having special kernels to handle the large-Df and large-Dt cases.


struct sd_clipper_transform_base : public wi_transform 
{
    using mask_t = simd_helpers::smask_t<float>;

    // (Frequency, time) downsampling factors and axis.
    const int nds_f;
    const int nds_t;
    const axis_type axis;
    
    // Clipping threshold.
    const double sigma;

    // Allocated in set_stream()
    float *tmp_sd = nullptr;
    mask_t *tmp_mask = nullptr;

    // Noncopyable
    sd_clipper_transform_base(const sd_clipper_transform_base &) = delete;
    sd_clipper_transform_base &operator=(const sd_clipper_transform_base &) = delete;

    sd_clipper_transform_base(int nds_f_, int nds_t_, axis_type axis_, int nt_chunk_, double sigma_)
	: nds_f(nds_f_), nds_t(nds_t_), axis(axis_), sigma(sigma_)
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
	free(tmp_mask);

	tmp_sd = nullptr;
	tmp_mask = nullptr;
    }

    virtual void set_stream(const wi_stream &stream) override
    {
	rf_assert(stream.nfreq % nds_f == 0);

	this->nfreq = stream.nfreq;
	this->allocate_tmp_arrays();
    }

    virtual void start_substream(int isubstream, double t0) override { }
    virtual void end_substream() override { }

    // Subclass must define allocate_tmp_arrays() and process_chunk()
    virtual void allocate_tmp_arrays() = 0;
};


// -------------------------------------------------------------------------------------------------
//
// sd_clipper_transform_time_axis


template<unsigned int S, unsigned int Df, unsigned int Dt>
struct sd_clipper_transform_time_axis : public sd_clipper_transform_base
{
    sd_clipper_transform_time_axis(int nt_chunk_, double sigma_)
	: sd_clipper_transform_base(Df, Dt, AXIS_TIME, nt_chunk_, sigma_)
    { }

    
    virtual void allocate_tmp_arrays() override
    {
	this->tmp_sd = aligned_alloc<float> (nfreq/Df);
	this->tmp_mask = aligned_alloc<mask_t> (nfreq/Df);
    }


    virtual void process_chunk(double t0, double t1, float *intensity, float *weights, ssize_t stride, float *pp_intensity, float *pp_weights, ssize_t pp_stride) override
    {
	throw runtime_error("sd_clipper_transform_time_axis::process_chunk() not written yet");
    }
};


// -------------------------------------------------------------------------------------------------
//
// Boilerplate needed to instantiate templates and export factory functions.


template<unsigned int S, unsigned int Df, unsigned int Dt>
inline shared_ptr<sd_clipper_transform_base> _make_std_dev_clipper3(axis_type axis, int nt_chunk, double sigma)
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
inline shared_ptr<sd_clipper_transform_base> _make_std_dev_clipper2(int Dt, axis_type axis, int nt_chunk, double sigma)
{
    throw runtime_error("rf_pipelines internal error: Dt=" + to_string(Dt) + " not found in std_dev_clipper template chain");
}

template<unsigned int S, unsigned int Df, unsigned int MaxDt, typename std::enable_if<(MaxDt > 0),int>::type = 0>
inline shared_ptr<sd_clipper_transform_base> _make_std_dev_clipper2(int Dt, axis_type axis, int nt_chunk, double sigma)
{
    if (Dt == MaxDt)
	return _make_std_dev_clipper3<S,Df,MaxDt> (axis, nt_chunk, sigma);
    else
	return _make_std_dev_clipper2<S,Df,(MaxDt/2)> (Dt, axis, nt_chunk, sigma);
}


template<unsigned int S, unsigned int MaxDf, unsigned int MaxDt, typename std::enable_if<(MaxDf==0),int>::type = 0>
inline shared_ptr<sd_clipper_transform_base> _make_std_dev_clipper(int Df, int Dt, axis_type axis, int nt_chunk, double sigma)
{
    throw runtime_error("rf_pipelines internal error: Df=" + to_string(Df) + " not found in std_dev_clipper template chain");
}

template<unsigned int S, unsigned int MaxDf, unsigned int MaxDt, typename std::enable_if<(MaxDf > 0),int>::type = 0>
inline shared_ptr<sd_clipper_transform_base> _make_std_dev_clipper(int Df, int Dt, axis_type axis, int nt_chunk, double sigma)
{
    if (Df == MaxDf)
	return _make_std_dev_clipper2<S,MaxDf,MaxDt> (Dt, axis, nt_chunk, sigma);
    else
	return _make_std_dev_clipper<S,(MaxDf/2),MaxDt> (Df, Dt, axis, nt_chunk, sigma);
}


// Externally callable "bottom line" routine which returns clipper transform with specified downsampling and axis.

shared_ptr<wi_transform> make_std_dev_clipper(int Df, int Dt, axis_type axis, int nt_chunk, double sigma)
{
    // SIMD word length on this machine (AVX instruction set assumed)
    static constexpr int S = 8;

    // MaxDf, MaxDt are the max downsampling factors allowed in the frequency and time directions.
    // If you get the error message "no precompiled kernel available..." then you'll need to change these (see below).
    // IMPORTANT: Df,Dt should be powers of two, or unpredictable problems will result!

    static constexpr int MaxDf = 32;
    static constexpr int MaxDt = 32;

    if ((Df <= 0) || !is_power_of_two(Df))
	throw runtime_error("rf_pipelines::make_std_dev_clipper(): Df must be a power of two (value received = " + to_string(Df) + ")");
    if ((Dt <= 0) || !is_power_of_two(Dt))
	throw runtime_error("rf_pipelines::make_std_dev_clipper(): Dt must be a power of two (value received = " + to_string(Dt) + ")");
    if (nt_chunk <= 0)
	throw runtime_error("rf_pipelines::make_std_dev_clipper(): nt_chunk must be > 0 (value received = " + to_string(nt_chunk) + ")");
    if (sigma < 2.0)
	throw runtime_error("rf_pipelines::make_std_dev_clipper(): sigma must be >= 2.0 (value received = " + to_string(sigma) + ")");
    
    if (nt_chunk % (Dt*S)) {
	stringstream ss;
	ss << "rf_pipelines::make_std_dev_clipper(): nt_chunk=" << nt_chunk << " is not a multiple of S*Dt\n"
	   << "Here, S=" << S << " is the simd length on this machine, and Dt=" << Dt << " is the requested time downsampling factor.\n"
	   << "Note that the value of S may be machine-dependent!\n";

	throw runtime_error(ss.str());
    }

    if ((Df > MaxDf) || (Dt > MaxDt)) {
	stringstream ss;
	ss << "rf_pipelines::make_std_dev_clipper(): no precompiled kernel is available for (Df,Dt)=(" << Df << "," << Dt << ")\n"
	   << "Eventually, this will be fixed in a systematic way.  As a temporary workaround, you can increase the values of\n"
	   << "MaxDf and MaxDt in rf_pipelines/clipper_transforms.cpp::make_std_dev_clipper() and recompile.\n";
	
	throw runtime_error(ss.str());
    }

    shared_ptr<sd_clipper_transform_base> ret = _make_std_dev_clipper<S,MaxDf,MaxDt> (Df, Dt, axis, nt_chunk, sigma);
    
    // Sanity check on the template instantiation
    assert(ret->nds_f == Df);
    assert(ret->nds_t == Dt);
    assert(ret->axis == axis);
    
    return ret;
}


}  // namespace rf_pipelines
