#include "rf_pipelines_internals.hpp"
#include "kernels/polyfit.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif


// -------------------------------------------------------------------------------------------------
//
// polynomial_detrender_t<T,S,N>, polynomial_detrender_f<T,S,N>:
//
// These transform classes (subclasses of wi_transform) perform detrending along either the time
// or frequency axis.  Note that the template parameter N is (polydeg+1), not (polydeg).


struct polynomial_detrender_base : public wi_transform
{
    double epsilon;

    polynomial_detrender_base(const char *transform_name, int S, int nt_chunk_, int polydeg, double epsilon_)
	: epsilon(epsilon_)
    {
	stringstream ss;
	ss << transform_name << "(nt_chunk=" << nt_chunk_ << ",polydeg=" << polydeg << ",epsilon=" << epsilon_ << ")";

	this->name = ss.str();
	this->nt_chunk = nt_chunk_;
	this->nt_prepad = 0;
	this->nt_postpad = 0;

	// No need to make these asserts "verbose", since they should have been checked in make_clipper2d().
	rf_assert(nt_chunk > 0 && (nt_chunk % S == 0));
	rf_assert(polydeg >= 0);
	rf_assert(epsilon > 0.0);
    }
    
    virtual void set_stream(const wi_stream &stream) override
    {
	this->nfreq = stream.nfreq;
    }

    virtual void start_substream(int isubstream, double t0) override { }
    virtual void end_substream() override { }
};


template<unsigned int S, unsigned int N>
struct polynomial_detrender_t : public polynomial_detrender_base
{
    polynomial_detrender_t(int nt_chunk_, double epsilon_)
	: polynomial_detrender_base("polynomial_detrender_time_axis", S, nt_chunk_, N-1, epsilon_)
    { }

    virtual void process_chunk(double t0, double t1, float *intensity, float *weights, ssize_t stride, float *pp_intensity, float *pp_weights, ssize_t pp_stride) override
    {
	_kernel_detrend_t<float,S,N> (nfreq, nt_chunk, intensity, weights, stride, epsilon);
    }
};


template<unsigned int S, unsigned int N>
struct polynomial_detrender_f : public polynomial_detrender_base
{
    polynomial_detrender_f(int nt_chunk_, double epsilon_)
	: polynomial_detrender_base("polynomial_detrender_freq_axis", S, nt_chunk_, N-1, epsilon_)
    { }

    virtual void process_chunk(double t0, double t1, float *intensity, float *weights, ssize_t stride, float *pp_intensity, float *pp_weights, ssize_t pp_stride) override
    {
	_kernel_detrend_f<float,S,N> (nfreq, nt_chunk, intensity, weights, stride, epsilon);
    }
};


// -------------------------------------------------------------------------------------------------
//
// Boilerplate needed to instantiate templates and export factory functions.


static void check_args(int nt_chunk, int polydeg, double epsilon, int S, int MaxDeg, const char *fn_name)
{
    if ((nt_chunk <= 0) || (nt_chunk % S) != 0) {
	stringstream ss;
	ss << "rf_pipelines::" << fn_name << "(): nt_chunk=" << nt_chunk << " must be a multiple of S=" << S << "\n"
	   << "Here, S is the simd vector size, and may vary between machines.";
	throw ss.str();
    }

    if (polydeg > MaxDeg) {
	stringstream ss;
	ss << "rf_pipelines::" << fn_name << "(): polydeg=" << polydeg << " must be <= " << MaxDeg << "\n"
	   << "The upper limit can be increased by editing rf_pipelines/polynomial_detrenders.cpp, changing\n"
	   << "MaxDeg in " << fn_name << "() and recompling rf_pipelines";
	throw ss.str();
    }

    if (epsilon <= 0.0)
	throw runtime_error(string("rf_pipelines::") + fn_name + ": epsilon must be > 0");
}


template<unsigned int S, unsigned int N, typename std::enable_if<(N==0),int>::type = 0>
static shared_ptr<wi_transform> _make_polynomial_detrender_t(int nt_chunk, int polydeg, double epsilon)
{
    throw runtime_error("make_polynomial_detrender_t: internal error");
}

template<unsigned int S, unsigned int N, typename std::enable_if<(N>0),int>::type = 0>
static shared_ptr<wi_transform> _make_polynomial_detrender_t(int nt_chunk, int polydeg, double epsilon)
{
    if (N == polydeg + 1)
	return make_shared< polynomial_detrender_t<S,N> > (nt_chunk, epsilon);
    return _make_polynomial_detrender_t<S,N-1> (nt_chunk, polydeg, epsilon);
}

shared_ptr<wi_transform> make_polynomial_detrender_time_axis(int nt_chunk, int polydeg, double epsilon)
{
    static constexpr int S = 8;
    static constexpr int MaxDeg = 20;

    check_args(nt_chunk, polydeg, epsilon, S, MaxDeg, "make_polynomial_detrender_time_axis");
    return _make_polynomial_detrender_t<S,MaxDeg+1> (nt_chunk, polydeg, epsilon);
}


template<unsigned int S, unsigned int N, typename std::enable_if<(N==0),int>::type = 0>
static shared_ptr<wi_transform> _make_polynomial_detrender_f(int nt_chunk, int polydeg, double epsilon)
{
    throw runtime_error("make_polynomial_detrender_f: internal error");
}

template<unsigned int S, unsigned int N, typename std::enable_if<(N>0),int>::type = 0>
static shared_ptr<wi_transform> _make_polynomial_detrender_f(int nt_chunk, int polydeg, double epsilon)
{
    if (N == polydeg + 1)
	return make_shared< polynomial_detrender_f<S,N> > (nt_chunk, epsilon);
    return _make_polynomial_detrender_f<S,N-1> (nt_chunk, polydeg, epsilon);
}

shared_ptr<wi_transform> make_polynomial_detrender_freq_axis(int nt_chunk, int polydeg, double epsilon)
{
    static constexpr int S = 8;
    static constexpr int MaxDeg = 20;

    check_args(nt_chunk, polydeg, epsilon, S, MaxDeg, "make_polynomial_detrender_freq_axis");
    return _make_polynomial_detrender_f<S,MaxDeg+1> (nt_chunk, polydeg, epsilon);
}


}  // namespace rf_pipelines
