#include "rf_pipelines_internals.hpp"
#include "kernels/intensity_clippers.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif


// -------------------------------------------------------------------------------------------------
//
// clipper2d_transform
//
// The compile-time arguments Df,Dt are the frequency/time downsampling factors, and the
// compile-time boolean argument IterFlag should be set to 'true' if and only if niter > 1.
//
// FIXME: currently we need to compile a new kernel for every (Df,Dt) pair.  Eventually I'd
// like to improve this by having special kernels to handle the large-Df and large-Dt cases.


template<unsigned int S, unsigned int Df, unsigned int Dt, bool IterFlag>
struct clipper2d_transform : public wi_transform 
{
    // These compile-time flags determine whether downsampled intensity/weights
    // arrays are written in _kernel_clip2d_wrms().
    static constexpr bool DsiFlag = (Df > 1) || (Dt > 1);
    static constexpr bool DswFlag = IterFlag && ((Df > 1) || (Dt > 1));

    const int niter;
    const double sigma;
    const double iter_sigma;

    float *ds_intensity = nullptr;
    float *ds_weights = nullptr;

    // Noncopyable
    clipper2d_transform(const clipper2d_transform<S,Df,Dt,IterFlag> &) = delete;
    clipper2d_transform<S,Df,Dt,IterFlag> operator&(const clipper2d_transform<S,Df,Dt,IterFlag> &) = delete;


    clipper2d_transform(int nt_chunk_, double sigma_, int niter_, double iter_sigma_)
	: niter(niter_), sigma(sigma_), iter_sigma(iter_sigma_ ? iter_sigma_ : sigma_)
    {
	stringstream ss;
	ss << "intensity_clipper_transform(Df=" << Df << ",Dt=" << Dt << ",axis=AXIS_NONE" 
	   << ", nt_chunk=" << nt_chunk_ << ",sigma=" << sigma << ",niter=" << niter 
	   << ",iter_sigma=" << iter_sigma << ")";

	this->name = ss.str();
	this->nt_chunk = nt_chunk_;
	this->nt_prepad = 0;
	this->nt_postpad = 0;

	// No need to make these asserts "verbose", since they should have been checked in make_intensity_clipper2d().
	rf_assert(sigma >= 1.0);
	rf_assert(iter_sigma >= 1.0);
	rf_assert(nt_chunk > 0);
	rf_assert(nt_chunk % Dt == 0);

	if (IterFlag) rf_assert(niter > 1);
	if (!IterFlag) rf_assert(niter == 1);
    }


    ~clipper2d_transform()
    {
	free(ds_intensity);
	free(ds_weights);
	ds_intensity = ds_weights = nullptr;
    }


    virtual void set_stream(const wi_stream &stream) override
    {
	rf_assert(stream.nfreq % Df == 0);

	this->nfreq = stream.nfreq;
	
	if (DsiFlag)
	    this->ds_intensity = aligned_alloc<float> ((nfreq/Df) * (nt_chunk/Dt));
	if (DswFlag)
	    this->ds_weights = aligned_alloc<float> ((nfreq/Df) * (nt_chunk/Dt));
    }


    virtual void process_chunk(double t0, double t1, float *intensity, float *weights, ssize_t stride, float *pp_intensity, float *pp_weights, ssize_t pp_stride) override
    {
	simd_t<float,S> mean, rms;
	_kernel_clip2d_wrms<float,S,Df,Dt,DsiFlag,DswFlag,float,S> (mean, rms, intensity, weights, nfreq, nt_chunk, stride, ds_intensity, ds_weights);

	const float *s_intensity = DsiFlag ? ds_intensity : intensity;
	const float *s_weights = DswFlag ? ds_weights : weights;
	int s_stride = DsiFlag ? (nt_chunk/Dt) : stride;   // must use DsiFlag here, not DswFlag
	
	for (int iter = 1; iter < niter; iter++) {
	    // (s_intensity, s_weights, iter_sigma) here
	    simd_t<float,S> thresh = simd_t<float,S>(iter_sigma) * rms;
	    _kernel_clip2d_iterate<float,S> (mean, rms, s_intensity, s_weights, mean, thresh, nfreq/Df, nt_chunk/Dt, s_stride);
	}

	// (s_intensity, weights, sigma) here
	simd_t<float,S> thresh = simd_t<float,S>(sigma) * rms;
	_kernel_clip2d_mask<float,S,Df,Dt> (weights, s_intensity, mean, thresh, nfreq, nt_chunk, stride, s_stride);
    }

    virtual void start_substream(int isubstream, double t0) override { }
    virtual void end_substream() override { }
};


// -------------------------------------------------------------------------------------------------
//
// clipper1d_t_transform
//
// Currently implemented by calling the 2d kernels many times with nfreq=1.
//
// Note: If we want to optimize the clipper transforms further, it might be worth exploring
// a generalization in which there is a template parameter R controlling the number of rows
// of the array (i.e. frequency channels) read in each pass.  In this case, we would need
// separate 1d_t and 2d kernels.


template<unsigned int S, unsigned int Df, unsigned int Dt, bool IterFlag>
struct clipper1d_t_transform : public wi_transform 
{
    // These compile-time flags determine whether downsampled intensity/weights
    // arrays are written in _kernel_clip2d_wrms().
    static constexpr bool DsiFlag = (Df > 1) || (Dt > 1);
    static constexpr bool DswFlag = IterFlag && ((Df > 1) || (Dt > 1));

    const int niter;
    const double sigma;
    const double iter_sigma;

    float *ds_intensity = nullptr;
    float *ds_weights = nullptr;

    // Noncopyable
    clipper1d_t_transform(const clipper1d_t_transform<S,Df,Dt,IterFlag> &) = delete;
    clipper1d_t_transform<S,Df,Dt,IterFlag> operator&(const clipper1d_t_transform<S,Df,Dt,IterFlag> &) = delete;


    clipper1d_t_transform(int nt_chunk_, double sigma_, int niter_, double iter_sigma_)
	: niter(niter_), sigma(sigma_), iter_sigma(iter_sigma_ ? iter_sigma_ : sigma_)
    {
	stringstream ss;
	ss << "intensity_clipper_transform(Df=" << Df << ",Dt=" << Dt << ",axis=AXIS_TIME" 
	   << ", nt_chunk=" << nt_chunk_ << ",sigma=" << sigma << ",niter=" << niter 
	   << ",iter_sigma=" << iter_sigma << ")";

	this->name = ss.str();
	this->nt_chunk = nt_chunk_;
	this->nt_prepad = 0;
	this->nt_postpad = 0;

	// No need to make these asserts "verbose", since they should have been checked in make_intensity_clipper2d().
	rf_assert(sigma >= 1.0);
	rf_assert(iter_sigma >= 1.0);
	rf_assert(nt_chunk > 0);
	rf_assert(nt_chunk % Dt == 0);

	if (IterFlag) rf_assert(niter > 1);
	if (!IterFlag) rf_assert(niter == 1);
    }


    ~clipper1d_t_transform()
    {
	free(ds_intensity);
	free(ds_weights);
	ds_intensity = ds_weights = nullptr;
    }


    virtual void set_stream(const wi_stream &stream) override
    {
	rf_assert(stream.nfreq % Df == 0);

	this->nfreq = stream.nfreq;

	// These are now 1d arrays.
	if (DsiFlag)
	    this->ds_intensity = aligned_alloc<float> (nt_chunk/Dt);
	if (DswFlag)
	    this->ds_weights = aligned_alloc<float> (nt_chunk/Dt);
    }


    virtual void process_chunk(double t0, double t1, float *intensity, float *weights, ssize_t stride, float *pp_intensity, float *pp_weights, ssize_t pp_stride) override
    {
	simd_t<float,S> mean, rms;

	for (int ifreq = 0; ifreq < nfreq; ifreq += Df) {
	    float *irow = intensity + ifreq * stride;
	    float *wrow = intensity + ifreq * stride;

	    // We use nfreq=1 and stride=0 here and throughout this routine.
	    _kernel_clip2d_wrms<float,S,Df,Dt,DsiFlag,DswFlag,float,S> (mean, rms, irow, wrow, 1, nt_chunk, 0, ds_intensity, ds_weights);
									
	    const float *irow2 = DsiFlag ? ds_intensity : irow;
	    const float *wrow2 = DswFlag ? ds_weights : wrow;
	
	    for (int iter = 1; iter < niter; iter++) {
		// (irow2, wrow2, iter_sigma)
		simd_t<float,S> thresh = simd_t<float,S>(iter_sigma) * rms;
		_kernel_clip2d_iterate<float,S> (mean, rms, irow2, wrow2, mean, thresh, 1, nt_chunk/Dt, 0);   // nfreq=1, stride=0
	    }

	    // (irow2, wrow, sigma)
	    simd_t<float,S> thresh = simd_t<float,S>(sigma) * rms;
	    _kernel_clip2d_mask<float,S,Df,Dt> (wrow, irow2, mean, thresh, 1, nt_chunk, 0, 0);
	}
    }

    virtual void start_substream(int isubstream, double t0) override { }
    virtual void end_substream() override { }
};


// -------------------------------------------------------------------------------------------------
//
// Boilerplate needed to instantiate templates and export factory functions.


template<unsigned int S, unsigned int Df, unsigned int Dt, unsigned int IterFlag>
inline shared_ptr<wi_transform> _make_intensity_clipper4(axis_type axis, int nt_chunk, double sigma, int niter, double iter_sigma)
{
    if (axis == AXIS_FREQ)
	throw runtime_error("rf_pipelines::make_intensity_clipper(): AXIS_FREQ is not implemented yet");
    if (axis == AXIS_TIME)
	return make_shared<clipper1d_t_transform<S,Df,Dt,IterFlag>> (nt_chunk, sigma, niter, iter_sigma);
    if (axis == AXIS_NONE)
	return make_shared<clipper2d_transform<S,Df,Dt,IterFlag>> (nt_chunk, sigma, niter, iter_sigma);

    throw runtime_error("rf_pipelines::make_intensity_clipper(): axis='" + to_string(axis) + "' is not a valid value");
}


template<unsigned int S, unsigned int Df, unsigned int Dt>
inline shared_ptr<wi_transform> _make_intensity_clipper3(axis_type axis, int nt_chunk, double sigma, int niter, double iter_sigma)
{
    if (niter > 1)
	return _make_intensity_clipper4<S,Df,Dt,true> (axis, nt_chunk, sigma, niter, iter_sigma);
    else
	return _make_intensity_clipper4<S,Df,Dt,false> (axis, nt_chunk, sigma, niter, iter_sigma);
}


template<unsigned int S, unsigned int Df, unsigned int MaxDt, typename std::enable_if<(MaxDt==0),int>::type = 0>
inline shared_ptr<wi_transform> _make_intensity_clipper2(int Dt, axis_type axis, int nt_chunk, double sigma, int niter, double iter_sigma)
{
    throw runtime_error("rf_pipelines internal error: Dt=" + to_string(Dt) + " not found in template chain");
}

template<unsigned int S, unsigned int Df, unsigned int MaxDt, typename std::enable_if<(MaxDt > 0),int>::type = 0>
inline shared_ptr<wi_transform> _make_intensity_clipper2(int Dt, axis_type axis, int nt_chunk, double sigma, int niter, double iter_sigma)
{
    if (Dt == MaxDt)
	throw _make_intensity_clipper3<S,Df,MaxDt> (axis, nt_chunk, sigma, niter, iter_sigma);
    else
	throw _make_intensity_clipper2<S,Df,(MaxDt/2)> (Dt, axis, nt_chunk, sigma, niter, iter_sigma);
}


template<unsigned int S, unsigned int MaxDf, unsigned int MaxDt, typename std::enable_if<(MaxDf==0),int>::type = 0>
inline shared_ptr<wi_transform> _make_intensity_clipper(int Df, int Dt, axis_type axis, int nt_chunk, double sigma, int niter, double iter_sigma)
{
    throw runtime_error("rf_pipelines internal error: Df=" + to_string(Df) + " not found in template chain");
}

template<unsigned int S, unsigned int MaxDf, unsigned int MaxDt, typename std::enable_if<(MaxDf > 0),int>::type = 0>
inline shared_ptr<wi_transform> _make_intensity_clipper(int Df, int Dt, axis_type axis, int nt_chunk, double sigma, int niter, double iter_sigma)
{
    if (Df == MaxDf)
	return _make_intensity_clipper2<S,MaxDf,MaxDt> (Dt, axis, nt_chunk, sigma, niter, iter_sigma);
    else
	return _make_intensity_clipper<S,(MaxDf/2),MaxDt> (Df, Dt, axis, nt_chunk, sigma, niter, iter_sigma);
}


// externally callable
shared_ptr<wi_transform> make_intensity_clipper(int Df, int Dt, axis_type axis, int nt_chunk, double sigma, int niter, double iter_sigma)
{
    // SIMD length on this machine
    static constexpr int S = 8;

    // MaxDf, MaxDt are the max downsampling factors allowed in the frequency and time directions.
    // If you get the error message "no precompiled kernel available..." then you'll need to change these (see below).
    // IMPORTANT: Df,Dt should be powers of two, or unpredictable problems will result!

    static constexpr int MaxDf = 32;
    static constexpr int MaxDt = 32;

    if ((Df <= 0) || !is_power_of_two(Df))
	throw runtime_error("rf_pipelines::make_intensity_clipper(): Df must be a power of two (value received = " + to_string(Df) + ")");
    if ((Dt <= 0) || !is_power_of_two(Dt))
	throw runtime_error("rf_pipelines::make_intensity_clipper(): Dt must be a power of two (value received = " + to_string(Dt) + ")");
    if (nt_chunk <= 0)
	throw runtime_error("rf_pipelines::make_intensity_clipper(): nt_chunk must be > 0 (value received = " + to_string(nt_chunk) + ")");
    if (sigma < 2.0)
	throw runtime_error("rf_pipelines::make_intensity_clipper(): sigma must be >= 2.0 (value received = " + to_string(sigma) + ")");
    if (niter < 1)
	throw runtime_error("rf_pipelines::make_intensity_clipper(): niter must be >= 1 (value received = " + to_string(niter) + ")");
    if (iter_sigma < 2.0)
	throw runtime_error("rf_pipelines::make_intensity_clipper(): iter_sigma must be >= 2.0 (value received = " + to_string(iter_sigma) + ")");
    
    if (nt_chunk % (Dt*S)) {
	stringstream ss;
	ss << "rf_pipelines::make_intensity_clipper(): nt_chunk=" << nt_chunk << " is not a multiple of S*Dt\n"
	   << "Here, S=" << S << " is the simd length on this machine, and Dt=" << Dt << " is the requested time downsampling factor.\n"
	   << "Note that the value of S may be machine-dependent!\n";

	throw runtime_error(ss.str());
    }

    if ((Df > MaxDf) || (Dt > MaxDt)) {
	stringstream ss;
	ss << "rf_pipelines::make_intensity_clipper(): no precompiled kernel is available for (Df,Dt)=(" << Df << "," << Dt << ")\n"
	   << "Eventually, this will be fixed in a systematic way.  As a temporary workaround, you can increase the values of\n"
	   << "MaxDf and MaxDt in rf_pipelines/clipper_transforms.cpp::make_intensity_clipper() and recompile.\n";
	
	throw runtime_error(ss.str());
    }

    return _make_intensity_clipper<S,MaxDf,MaxDt> (Df, Dt, axis, nt_chunk, sigma, niter, iter_sigma);
}


}  // namespace rf_pipelines
