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
	ss << "intensity_clipper2d_transform(Df=" << Df << ",Dt=" << Dt << ",nt_chunk=" << nt_chunk_
	   << ",sigma=" << sigma << ",niter=" << niter << ",iter_sigma=" << iter_sigma << ")";

	this->name = ss.str();
	this->nt_chunk = nt_chunk_;
	this->nt_prepad = 0;
	this->nt_postpad = 0;

	// No need to make these asserts "verbose", since they should have been checked in make_clipper2d().
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
	int s_stride = DsiFlag ? stride : (nt_chunk/Dt);   // must use DsiFlag here, not DswFlag

	for (int iter = 1; iter < niter; iter++) {
	    simd_t<float,S> thresh = simd_t<float,S>(iter_sigma) * rms;   // iter_sigma here
	    _kernel_clip2d_iterate<float,S> (mean, rms, s_intensity, s_weights, mean, thresh, nfreq/Df, nt_chunk/Dt, s_stride);
	}

	simd_t<float,S> thresh = simd_t<float,S>(sigma) * rms;           // sigma here
	_kernel_clip2d_mask<float,S,Df,Dt> (weights, s_intensity, mean, thresh, nfreq, nt_chunk, stride, s_stride);
    }

    virtual void start_substream(int isubstream, double t0) override { }
    virtual void end_substream() override { }
};


// -------------------------------------------------------------------------------------------------
//
// Boilerplate needed to instantiate templates and export factory function.


template<unsigned int S, unsigned int Df, unsigned int MaxDt, typename std::enable_if<(MaxDt==0),int>::type = 0>
static inline shared_ptr<wi_transform> _make_clipper2d_df(int Dt, int nt_chunk, double sigma, int niter, double iter_sigma)
{
    throw runtime_error("rf_pipelines: internal error in _make_clipper2d_df()");
}

template<unsigned int S, unsigned int Df, unsigned int MaxDt, typename std::enable_if<(MaxDt>0),int>::type = 0>
static inline shared_ptr<wi_transform> _make_clipper2d_df(int Dt, int nt_chunk, double sigma, int niter, double iter_sigma)
{
    if (Dt != MaxDt) 
	return _make_clipper2d_df<S,Df,(MaxDt/2)> (Dt, nt_chunk, sigma, niter, iter_sigma);
    if (niter > 1)
	return make_shared< clipper2d_transform<S,Df,MaxDt,true> >(nt_chunk, sigma, niter, iter_sigma);
    else
	return make_shared< clipper2d_transform<S,Df,MaxDt,false> >(nt_chunk, sigma, niter, iter_sigma);
}

template<unsigned int S, unsigned int MaxDf, unsigned int MaxDt, typename std::enable_if<(MaxDf==0),int>::type = 0>
static inline shared_ptr<wi_transform> _make_clipper2d(int Df, int Dt, int nt_chunk, double sigma, int niter, double iter_sigma)
{
    throw runtime_error("rf_pipelines: internal error in _make_clipper2d()");
}

template<unsigned int S, unsigned int MaxDf, unsigned int MaxDt, typename std::enable_if<(MaxDf>0),int>::type = 0>
static inline shared_ptr<wi_transform> _make_clipper2d(int Df, int Dt, int nt_chunk, double sigma, int niter, double iter_sigma)
{
    if (Df == MaxDf) 
	return _make_clipper2d_df<S,MaxDf,MaxDt> (Dt, nt_chunk, sigma, niter, iter_sigma);
    return _make_clipper2d<S,(MaxDf/2),MaxDt> (Df, Dt, nt_chunk, sigma, niter, iter_sigma);
};


shared_ptr<wi_transform> make_intensity_clipper2d(int Df, int Dt, int nt_chunk, double sigma, int niter, double iter_sigma)
{
    // SIMD length on this machine
    static constexpr int S = 8;

    // MaxDf, MaxDt are the max downsampling factors allowed in the frequency and time directions.
    // If you get the error message "no precompiled kernel available..." then you'll need to change these (see below).
    // Increasing Df,Dt will increase compile time and object file size, but this shouldn't be a big deal.
    // IMPORTANT: Df,Dt should be powers of two, or unpredictable problems will result!
    static constexpr int MaxDf = 32;
    static constexpr int MaxDt = 32;

    if ((Df <= 0) || !is_power_of_two(Df))
	throw runtime_error("rf_pipelines::make_intensity_clipper2d(): Df must be a power of two (value received = " + to_string(Df) + ")");
    if ((Dt <= 0) || !is_power_of_two(Dt))
	throw runtime_error("rf_pipelines::make_intensity_clipper2d(): Dt must be a power of two (value received = " + to_string(Dt) + ")");
    if (nt_chunk <= 0)
	throw runtime_error("rf_pipelines::make_intensity_clipper2d(): nt_chunk must be > 0 (value received = " + to_string(nt_chunk) + ")");
    if (sigma < 2.0)
	throw runtime_error("rf_pipelines::make_intensity_clipper2d(): sigma must be >= 2.0 (value received = " + to_string(sigma) + ")");
    if (niter < 1)
	throw runtime_error("rf_pipelines::make_intensity_clipper2d(): niter must be >= 1 (value received = " + to_string(niter) + ")");
    if (iter_sigma < 2.0)
	throw runtime_error("rf_pipelines::make_intensity_clipper2d(): iter_sigma must be >= 2.0 (value received = " + to_string(iter_sigma) + ")");
    
    if (nt_chunk % (Dt*S)) {
	stringstream ss;
	ss << "rf_pipelines::make_intensity_clipper2d(): nt_chunk=" << nt_chunk << " is not a multiple of S*Dt\n"
	   << "Here, S=" << S << " is the simd length on this machine, and Dt=" << Dt << " is the requested time downsampling factor.\n"
	   << "Note that the value of S may be machine-dependent!\n";

	throw runtime_error(ss.str());
    }

    if ((Df > MaxDf) || (Dt > MaxDt)) {
	stringstream ss;
	ss << "rf_pipelines::make_intensity_clipper2d(): no precompiled kernel is available for (Df,Dt)=(" << Df << "," << Dt << ")\n"
	   << "Eventually, this will be fixed in a systematic way.  As a temporary workaround, you can change the values of\n"
	   << "MaxDf and MaxDt in rf_pipelines/clipper_transforms.cpp::make_intensity_clipper2d() and recompile.\n";
	
	throw runtime_error(ss.str());
    }

    return _make_clipper2d<S,MaxDf,MaxDt> (Df, Dt, nt_chunk, sigma, niter, iter_sigma);
}


}  // namespace rf_pipelines
