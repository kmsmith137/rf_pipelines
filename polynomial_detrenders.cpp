// FIXME (low-priority) a nuisance issue when working with this code is that functions
// which are very similar have different argument orderings, e.g.
//
//          make_polynomial_detrender(nt_chunk, axis, polydeg, epsilon)
//   calls _make_polynomial_detrender(axis, nt_chunk, polydeg, epsilon)

#include <cassert>
#include <array>

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


// Helper function to check parameters passed to a detrending call
static void check_params(axis_type axis, int nfreq, int nt, int stride, int polydeg, double epsilon)
{
    static constexpr int MaxDeg = constants::polynomial_detrender_max_degree;
    static constexpr int S = constants::single_precision_simd_length;

    if (_unlikely((axis != AXIS_FREQ) && (axis != AXIS_TIME)))
	throw runtime_error("rf_pipelines polynomial detrender: axis=" + stringify(axis) + " is not defined for this transform");

    if (_unlikely(nfreq <= 0))
	throw runtime_error("rf_pipelines polynomial detrender: nfreq=" + to_string(nfreq) + ", positive number expected");

    if (_unlikely(nt <= 0))
	throw runtime_error("rf_pipelines polynomial detrender: nt_chunk=" + to_string(nt) + ", positive number expected");
    
    if (_unlikely((nt % S) != 0))
	throw runtime_error("rf_pipelines polynomial detrender: nt_chunk=" + to_string(nt)
			    + " must be a multiple of constants::single_precision_simd_length=" + to_string(S));
    
    if (_unlikely(abs(stride) < nt))
	throw runtime_error("rf_pipelines polynomial detrender: stride=" + to_string(stride) + " must be >= nt");

    if (_unlikely(polydeg < 0))
	throw runtime_error("rf_pipelines polynomial detrender: polydeg=" + to_string(polydeg) + ", positive number expected");

    if (_unlikely(polydeg > MaxDeg))
	throw runtime_error("rf_pipelines polynomial detrender: polydeg=" + to_string(polydeg)
			    + " exceeds compile time limits; to fix this see 'constants' in rf_pipelines.hpp");

    if (_unlikely(epsilon <= 0.0))
	throw runtime_error("rf_pipelines polynomial detrender: epsilon=" + to_string(epsilon) + ", positive number expected");
}


template<unsigned int S, unsigned int N, typename std::enable_if<(N==0),int>::type = 0>
static shared_ptr<polynomial_detrender_base> _make_polynomial_detrender2(axis_type axis, int nt_chunk, int polydeg, double epsilon)
{
    throw runtime_error("rf_pipelines::make_polynomial_detrender: internal error");
}

template<unsigned int S, unsigned int N, typename std::enable_if<(N>0),int>::type = 0>
static shared_ptr<polynomial_detrender_base> _make_polynomial_detrender2(axis_type axis, int nt_chunk, int polydeg, double epsilon)
{
    if (N != polydeg + 1)
	return _make_polynomial_detrender2<S,N-1> (axis, nt_chunk, polydeg, epsilon);

    if (axis == AXIS_FREQ)
	return make_shared<polynomial_detrender_freq_axis<S,N> > (nt_chunk, epsilon);
    if (axis == AXIS_TIME)
	return make_shared<polynomial_detrender_time_axis<S,N> > (nt_chunk, epsilon);
    if (axis == AXIS_NONE)
	throw runtime_error("rf_pipelines::make_polynomial_detrender(): axis=None is not supported for this transform");

    throw runtime_error("rf_pipelines::make_polynomial_detrender(): axis='" + to_string(axis) + "' is not a valid value");
}


// Shouldn't be declared static, in order to avoid potentially instantiating every template twice.
shared_ptr<polynomial_detrender_base> _make_polynomial_detrender(axis_type axis, int nt_chunk, int polydeg, double epsilon)
{
    static constexpr int MaxDeg = constants::polynomial_detrender_max_degree;
    static constexpr int S = constants::single_precision_simd_length;

    shared_ptr<polynomial_detrender_base> ret = _make_polynomial_detrender2<S,MaxDeg+1> (axis, nt_chunk, polydeg, epsilon);

    // Sanity check the template instantiation
    assert(ret->axis == axis);
    assert(ret->nt_chunk == nt_chunk);
    assert(ret->polydeg == polydeg);
    assert(ret->epsilon == epsilon);

    return ret;
}


// Externally callable factory function
shared_ptr<wi_transform> make_polynomial_detrender(int nt_chunk, axis_type axis, int polydeg, double epsilon)
{
    int dummy_nfreq = 16;         // arbitrary
    int dummy_stride = nt_chunk;  // arbitrary

    check_params(axis, dummy_nfreq, nt_chunk, dummy_stride, polydeg, epsilon);
    return _make_polynomial_detrender(axis, nt_chunk, polydeg, epsilon);
}


// -------------------------------------------------------------------------------------------------


struct detrending_kernel_table {
    static constexpr int MaxDeg = constants::polynomial_detrender_max_degree;

    using ktab2_t = std::array<detrending_kernel_t, (MaxDeg+1)>;
    using ktab_t = std::array<ktab2_t, 2>;

    ktab_t entries;

    detrending_kernel_table()
    {
	for (axis_type axis: { AXIS_FREQ, AXIS_TIME }) {
	    for (int polydeg = 0; polydeg <= MaxDeg; polydeg++) {
		auto t = _make_polynomial_detrender (axis, 16, polydeg, 1.0e-2);   // (nt_chunk, epsilon) arbitrary here
		entries.at(axis).at(polydeg) = t->kernel;
	    }
	}
    }

    inline detrending_kernel_t get(axis_type axis, int polydeg)
    {
	// Bounds-checked
	return entries.at(axis).at(polydeg);
    }
};


void apply_polynomial_detrender(float *intensity, const float *weights, int nfreq, int nt, int stride, axis_type axis, int polydeg, double epsilon)
{
    static detrending_kernel_table ktab;

    check_params(axis, nfreq, nt, stride, polydeg, epsilon);

    detrending_kernel_t k = ktab.get(axis, polydeg);
    k(nfreq, nt, intensity, const_cast<float *> (weights), stride, epsilon);
}


}  // namespace rf_pipelines
