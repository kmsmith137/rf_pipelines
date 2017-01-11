#include <cassert>

#include "rf_pipelines_internals.hpp"
#include "kernels/polyfit.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif


// Usage: kernel(nfreq, nt, intensity, weights, stride, epsilon)
using detrending_kernel_t = void (*)(int, int, float *, float *, int, double);


// -------------------------------------------------------------------------------------------------
//
// polynomial_detrender_t<T,S,N>, polynomial_detrender_f<T,S,N>:
//
// These transform classes (subclasses of wi_transform) perform detrending along either the time
// or frequency axis.  Note that the template parameter N is (polydeg+1), not (polydeg).


struct polynomial_detrender_base : public wi_transform
{
    const int axis;
    const int polydeg;
    const double epsilon;
    const detrending_kernel_t kernel;

    polynomial_detrender_base(int axis_, int nt_chunk_, int polydeg_, double epsilon_, detrending_kernel_t kernel_, int S) :
	axis(axis_), polydeg(polydeg_), epsilon(epsilon_), kernel(kernel_)
    {
	stringstream ss;
	ss << "polynomial_detrender(axis=" << axis << ",nt_chunk=" << nt_chunk_ << ",polydeg=" << polydeg << ",epsilon=" << epsilon_ << ")";

	this->name = ss.str();
	this->nt_chunk = nt_chunk_;
	this->nt_prepad = 0;
	this->nt_postpad = 0;

	// No need to make these asserts "verbose", since they should have been checked in make_polynomial_detrender().
	rf_assert(nt_chunk > 0 && (nt_chunk % S == 0));
	rf_assert(polydeg >= 0);
	rf_assert(epsilon > 0.0);
    }
    
    virtual void set_stream(const wi_stream &stream) override
    {
	this->nfreq = stream.nfreq;
    }

    virtual void process_chunk(double t0, double t1, float *intensity, float *weights, ssize_t stride, float *pp_intensity, float *pp_weights, ssize_t pp_stride) override
    {
	this->kernel(nfreq, nt_chunk, intensity, weights, stride, epsilon);
    }

    virtual void start_substream(int isubstream, double t0) override { }
    virtual void end_substream() override { }
};


template<unsigned int S, unsigned int N>
struct polynomial_detrender_time_axis : public polynomial_detrender_base
{
    polynomial_detrender_time_axis(int nt_chunk_, double epsilon_) :
	polynomial_detrender_base(AXIS_TIME, nt_chunk_, N-1, epsilon_, _kernel_detrend_t<float,S,N>, S)
    { }
};


template<unsigned int S, unsigned int N>
struct polynomial_detrender_freq_axis : public polynomial_detrender_base
{
    polynomial_detrender_freq_axis(int nt_chunk_, double epsilon_) :
	polynomial_detrender_base(AXIS_FREQ, nt_chunk_, N-1, epsilon_, _kernel_detrend_f<float,S,N>, S)
    { }
};


// -------------------------------------------------------------------------------------------------
//
// Boilerplate needed to instantiate templates and export factory functions.


template<unsigned int S, unsigned int N, typename std::enable_if<(N==0),int>::type = 0>
static shared_ptr<polynomial_detrender_base> _make_polynomial_detrender(axis_type axis, int nt_chunk, int polydeg, double epsilon)
{
    throw runtime_error("rf_pipelines::make_polynomial_detrender: internal error");
}

template<unsigned int S, unsigned int N, typename std::enable_if<(N>0),int>::type = 0>
static shared_ptr<polynomial_detrender_base> _make_polynomial_detrender(axis_type axis, int nt_chunk, int polydeg, double epsilon)
{
    if (N != polydeg + 1)
	return _make_polynomial_detrender<S,N-1> (axis, nt_chunk, polydeg, epsilon);

    if (axis == AXIS_FREQ)
	return make_shared<polynomial_detrender_freq_axis<S,N> > (nt_chunk, epsilon);
    if (axis == AXIS_TIME)
	return make_shared<polynomial_detrender_time_axis<S,N> > (nt_chunk, epsilon);
    if (axis == AXIS_NONE)
	throw runtime_error("rf_pipelines::make_polynomial_detrender(): axis=None is not supported for this transform");

    throw runtime_error("rf_pipelines::make_polynomial_detrender(): axis='" + to_string(axis) + "' is not a valid value");
}


// Externally callable factory function
shared_ptr<wi_transform> make_polynomial_detrender(axis_type axis, int nt_chunk, int polydeg, double epsilon)
{
    // S = single-precision simd vector length on this machine (AVX instruction set assumed)
    static constexpr int S = 8;
    static constexpr int MaxDeg = 20;

    if ((nt_chunk <= 0) || (nt_chunk % S) != 0) {
	stringstream ss;
	ss << "rf_pipelines::make_polynomial_detrender(): nt_chunk=" << nt_chunk << " must be a multiple of S=" << S << "\n"
	   << "Here, S is the simd vector size, and may vary between machines.";
	throw ss.str();
    }

    if (polydeg > MaxDeg) {
	stringstream ss;
	ss << "rf_pipelines::make_polynomial_detrender(): polydeg=" << polydeg << " must be <= " << MaxDeg << "\n"
	   << "The upper limit can be increased by editing rf_pipelines/polynomial_detrenders.cpp, changing\n"
	   << "MaxDeg in make_polynomial_detrender() and recompling rf_pipelines";
	throw ss.str();
    }

    if (epsilon <= 0.0)
	throw runtime_error(string("rf_pipelines::make_polynomial_detrender(): epsilon must be > 0"));

    shared_ptr<polynomial_detrender_base> ret = _make_polynomial_detrender<S,MaxDeg+1> (axis, nt_chunk, polydeg, epsilon);
    
    assert(ret->axis == axis);
    assert(ret->nt_chunk == nt_chunk);
    assert(ret->polydeg == polydeg);
    assert(ret->epsilon == epsilon);

    return ret;
}


}  // namespace rf_pipelines
