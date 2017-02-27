// FIXME (low-priority) a nuisance issue when working with this code is that functions
// which are very similar have different argument orderings, e.g.
//
//          make_polynomial_detrender(nt_chunk, axis, polydeg, epsilon)
//   calls _make_polynomial_detrender(axis, nt_chunk, polydeg, epsilon)

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


struct polynomial_detrender : public wi_transform
{
    const int axis;
    const int polydeg;
    const double epsilon;
    const detrending_kernel_t kernel;

    polynomial_detrender(int axis_, int nt_chunk_, int polydeg_, double epsilon_, detrending_kernel_t kernel_) :
	axis(axis_), polydeg(polydeg_), epsilon(epsilon_), kernel(kernel_)
    {
	stringstream ss;
        ss << "polynomial_detrender_cpp(nt_chunk=" << nt_chunk_ << ", axis=" << axis << ", polydeg=" << polydeg << ", epsilon=" << epsilon_ << ")";

	this->name = ss.str();
	this->nt_chunk = nt_chunk_;
	this->nt_prepad = 0;
	this->nt_postpad = 0;
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


// -------------------------------------------------------------------------------------------------
//
// _fill_detrending_kernel_table<S,N>(): fills shape (N,2) array with kernels.
// The outer index is a polynomial degree 0 <= polydeg < N, and the inner index is the axis.


template<unsigned int S, unsigned int N, typename std::enable_if<(N==0),int>::type = 0>
inline void fill_detrending_kernel_table(detrending_kernel_t *out) { }

template<unsigned int S, unsigned int N, typename std::enable_if<(N>0),int>::type = 0>
inline void fill_detrending_kernel_table(detrending_kernel_t *out)
{
    static_assert(AXIS_FREQ == 0, "polynomial_detrenders: current implementation assumes AXIS_FREQ==0");
    static_assert(AXIS_TIME == 1, "polynomial_detrenders: current implementation assumes AXIS_TIME==1");

    fill_detrending_kernel_table<S,N-1> (out);
    out[2*(N-1) + AXIS_FREQ] = _kernel_detrend_f<float,S,N>;
    out[2*(N-1) + AXIS_TIME] = _kernel_detrend_t<float,S,N>;
}


struct detrending_kernel_table {
    static constexpr int S = constants::single_precision_simd_length;
    static constexpr int MaxDeg = constants::polynomial_detrender_max_degree;

    std::vector<detrending_kernel_t> entries;

    detrending_kernel_table() : entries(2*MaxDeg+2)
    {
	fill_detrending_kernel_table<S,MaxDeg+1> (&entries[0]);
    }

    // Caller must argument-check by calling check_params()!
    inline detrending_kernel_t get_kernel(int axis, int polydeg)
    {
	return entries[2*polydeg + axis];
    }
};


static detrending_kernel_table global_detrending_kernel_table;


// -------------------------------------------------------------------------------------------------


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


// Externally callable factory function
shared_ptr<wi_transform> make_polynomial_detrender(int nt_chunk, axis_type axis, int polydeg, double epsilon)
{
    int dummy_nfreq = 16;         // arbitrary
    int dummy_stride = nt_chunk;  // arbitrary

    check_params(axis, dummy_nfreq, nt_chunk, dummy_stride, polydeg, epsilon);

    detrending_kernel_t kernel = global_detrending_kernel_table.get_kernel(axis, polydeg);
    return make_shared<polynomial_detrender> (axis, nt_chunk, polydeg, epsilon, kernel);
}


void apply_polynomial_detrender(float *intensity, float *weights, int nfreq, int nt, int stride, axis_type axis, int polydeg, double epsilon)
{
    check_params(axis, nfreq, nt, stride, polydeg, epsilon);

    if (_unlikely(!intensity))
	throw runtime_error("rf_pipelines: apply_polynomial_detrender(): NULL intensity pointer");
    if (_unlikely(!weights))
	throw runtime_error("rf_pipelines: apply_polynomial_detrender(): NULL weights pointer");

    detrending_kernel_t kernel = global_detrending_kernel_table.get_kernel(axis, polydeg);
    kernel(nfreq, nt, intensity, weights, stride, epsilon);
}


}  // namespace rf_pipelines
