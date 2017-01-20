// FIXME: currently we need to compile a new kernel for every (Df,Dt) pair, where
// Df,Dt are the frequency/time downsampling factors.  Eventually I'd like to 
// improve this by having special kernels to handle the large-Df and large-Dt cases.

#include <array>
#include <cassert>

#include "rf_pipelines_internals.hpp"

#include "kernels/mask.hpp"
#include "kernels/std_dev_clippers.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif


using mask_t = simd_helpers::smask_t<float>;


// Externally-linkable helper function, declared "extern" in kernels/std_dev_clippers.hpp.
// If the arguments change here, then declaration should be changed there as well!
void clip_1d(int n, float *tmp_sd, mask_t *tmp_valid, double sigma)
{
#if 0
    cerr << "clip_1d: [";
    for (int i = 0; i < n; i++) {
	if (tmp_valid[i])
	    cerr << " " << tmp_sd[i];
	else
	    cerr << " -";
    }
    cerr << " ]\n";
#endif

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
    
    float stdv = sqrtf(acc2/acc0);
    float thresh = sigma * stdv;

    for (int i = 0; i < n; i++) {
	if (fabs(tmp_sd[i] - mean) >= thresh)
	    tmp_valid[i] = 0;
    }
}


template<typename T>
inline void allocate_buffers(std_dev_clipper_buffers<T> &buf, int nfreq, int nt, axis_type axis, int Df, int Dt, bool two_pass)
{
    int sd_nalloc = (axis == AXIS_FREQ) ? (nt/Dt) : (nfreq/Df);

    buf.sd = aligned_alloc<T> (sd_nalloc);
    buf.sd_valid = aligned_alloc<smask_t<T,1>> (sd_nalloc);

    if (two_pass && ((Df > 1) || (Dt > 1))) {
	buf.ds_int = aligned_alloc<T> ((nfreq*nt) / (Df*Dt));
	buf.ds_wt = aligned_alloc<T> ((nfreq*nt) / (Df*Dt));
    }
}


// -------------------------------------------------------------------------------------------------
//
// std_dev_clipper_kernel_table


// kernel(intensity, weights, nfreq, nt, stride, sigma, tmp_sd, tmp_valid)
using std_dev_clipper_kernel_t = void (*)(const std_dev_clipper_buffers<float> &, const float *, float *, int, int, int, double);


// Fills shape-(NDt,2,2) array indexed by (Dt,axis,two_pass)
template<unsigned int S, unsigned int Df, unsigned int NDt, typename enable_if<(NDt==0),int>::type = 0>
inline void fill_3d_std_dev_clipper_kernel_table(std_dev_clipper_kernel_t *out) { }

template<unsigned int S, unsigned int Df, unsigned int NDt, typename enable_if<(NDt>0),int>::type = 0>
inline void fill_3d_std_dev_clipper_kernel_table(std_dev_clipper_kernel_t *out) 
{ 
    static_assert(AXIS_FREQ == 0, "expected AXIS_FREQ==0");
    static_assert(AXIS_TIME == 1, "expected AXIS_TIME==1");

    fill_3d_std_dev_clipper_kernel_table<S,Df,NDt-1> (out);

    constexpr unsigned int Dt = 1 << (NDt-1);
    out[4*(NDt-1) + 2*AXIS_FREQ] = _kernel_std_dev_clip_freq_axis<float,S,Df,Dt,false>;
    out[4*(NDt-1) + 2*AXIS_FREQ+1] = _kernel_std_dev_clip_freq_axis<float,S,Df,Dt,true>;
    out[4*(NDt-1) + 2*AXIS_TIME] = _kernel_std_dev_clip_time_axis<float,S,Df,Dt,false>;
    out[4*(NDt-1) + 2*AXIS_TIME+1] = _kernel_std_dev_clip_time_axis<float,S,Df,Dt,true>;
}

// Fills shape-(NDf,NDt,2,2) array indexed by (Df,Dt,axis,two_pass)
template<unsigned int S, unsigned int NDf, unsigned int NDt, typename enable_if<(NDf==0),int>::type = 0>
inline void fill_4d_std_dev_clipper_kernel_table(std_dev_clipper_kernel_t *out) { }

template<unsigned int S, unsigned int NDf, unsigned int NDt, typename enable_if<(NDf>0),int>::type = 0>
inline void fill_4d_std_dev_clipper_kernel_table(std_dev_clipper_kernel_t *out) 
{ 
    fill_4d_std_dev_clipper_kernel_table<S,NDf-1,NDt> (out);
    fill_3d_std_dev_clipper_kernel_table<S,(1<<(NDf-1)),NDt> (out + 4*(NDf-1)*NDt);
}


struct std_dev_clipper_kernel_table {
    static constexpr int S = constants::single_precision_simd_length;
    static constexpr int MaxDf = constants::max_frequency_downsampling;
    static constexpr int MaxDt = constants::max_time_downsampling;
    static constexpr int NDf = IntegerLog2<MaxDf>() + 1;
    static constexpr int NDt = IntegerLog2<MaxDt>() + 1;

    // Weird: for some reason using std::max() here gives a clang linker (not compiler) error.
    static constexpr int MaxD = (MaxDf > MaxDt) ? MaxDf : MaxDt;

    vector<std_dev_clipper_kernel_t> kernels;

    integer_log2_lookup_table ilog2_lookup;

    std_dev_clipper_kernel_table() :
	kernels(4*NDf*NDt), ilog2_lookup(MaxD)
    {
	fill_4d_std_dev_clipper_kernel_table<S,NDf,NDt> (&kernels[0]);
    }

    // Caller must call check_params()!
    inline std_dev_clipper_kernel_t get_kernel(axis_type axis, int Df, int Dt, bool two_pass)
    {
	int idf = ilog2_lookup(Df);
	int idt = ilog2_lookup(Dt);

	return kernels[4*(idf*NDt+idt) + 2*axis + (two_pass ? 1 : 0)];
    }
};


static std_dev_clipper_kernel_table global_std_dev_clipper_kernel_table;


// -------------------------------------------------------------------------------------------------
//
// sd_clipper_transform


struct std_dev_clipper_transform : public wi_transform 
{
    // (Frequency, time) downsampling factors and axis.
    const int nds_f;
    const int nds_t;
    const axis_type axis;
    const bool two_pass;
    
    // Clipping threshold.
    const double sigma;

    std_dev_clipper_kernel_t kernel;
    std_dev_clipper_buffers<float> buf;

    // Noncopyable
    std_dev_clipper_transform(const std_dev_clipper_transform &) = delete;
    std_dev_clipper_transform &operator=(const std_dev_clipper_transform &) = delete;

    std_dev_clipper_transform(int nds_f_, int nds_t_, axis_type axis_, int nt_chunk_, double sigma_, bool two_pass_, std_dev_clipper_kernel_t(kernel_))
	: nds_f(nds_f_), nds_t(nds_t_), axis(axis_), two_pass(two_pass_), sigma(sigma_), kernel(kernel_)
    {
	stringstream ss;
        ss << "std_dev_clipper_transform_cpp(nt_chunk=" << nt_chunk_ << ", axis=" << axis
           << ", sigma=" << sigma << ", Df=" << nds_f << ", Dt=" << nds_t << ", two_pass=" << two_pass << ")";
	
        this->name = ss.str();
	this->nt_chunk = nt_chunk_;
	this->nt_prepad = 0;
	this->nt_postpad = 0;
    }

    virtual void set_stream(const wi_stream &stream) override
    {
	if (stream.nfreq % nds_f)
	    throw runtime_error("rf_pipelines std_dev_clipper: stream nfreq (=" + to_string(stream.nfreq) 
				+ ") is not divisible by frequency downsampling factor Df=" + to_string(nds_f));

	this->nfreq = stream.nfreq;
	
	allocate_buffers(this->buf, nfreq, nt_chunk, axis, nds_f, nds_t, two_pass);
    }

    virtual void process_chunk(double t0, double t1, float *intensity, float *weights, ssize_t stride, float *pp_intensity, float *pp_weights, ssize_t pp_stride) override
    {
	this->kernel(buf, intensity, weights, nfreq, nt_chunk, stride, sigma);
    }

    virtual void start_substream(int isubstream, double t0) override { }
    virtual void end_substream() override { }
};


// -------------------------------------------------------------------------------------------------


static void check_params(int Df, int Dt, axis_type axis, int nfreq, int nt, int stride, double sigma)
{
    static constexpr int S = constants::single_precision_simd_length;
    static constexpr int MaxDf = constants::max_frequency_downsampling;
    static constexpr int MaxDt = constants::max_time_downsampling;

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

    if (_unlikely(sigma < 1.0))
	throw runtime_error("rf_pipelines std_dev clipper: sigma=" + to_string(sigma) + " must be >= 1.0");

    if (_unlikely((nfreq % Df) != 0))
	throw runtime_error("rf_pipelines std_dev clipper: nfreq=" + to_string(nfreq)
			    + " must be a multiple of the downsampling factor Df=" + to_string(Df));

    if (_unlikely((nt % (Dt*S)) != 0))
	throw runtime_error("rf_pipelines std_dev clipper: nt=" + to_string(nt)
			    + " must be a multiple of the downsampling factor Dt=" + to_string(Dt)
			    + " multiplied by constants::single_precision_simd_length=" + to_string(S));

    if (_unlikely((Df > MaxDf) || (Dt > MaxDt)))
	throw runtime_error("rf_pipelines std_dev clipper: (Df,Dt)=(" + to_string(Df) + "," + to_string(Dt) + ")"
			    + " exceeds compile time limits; to fix this see 'constants' in rf_pipelines.hpp");
}


// Externally callable
shared_ptr<wi_transform> make_std_dev_clipper(int nt_chunk, axis_type axis, double sigma, int Df, int Dt, bool two_pass)
{
    int dummy_nfreq = Df;         // arbitrary
    int dummy_stride = nt_chunk;  // arbitrary
    check_params(Df, Dt, axis, dummy_nfreq, nt_chunk, dummy_stride, sigma);

    auto kernel = global_std_dev_clipper_kernel_table.get_kernel(axis, Df, Dt, two_pass);
    return make_shared<std_dev_clipper_transform> (Df, Dt, axis, nt_chunk, sigma, two_pass, kernel);
}


// Externally callable
void apply_std_dev_clipper(const float *intensity, float *weights, int nfreq, int nt, int stride, axis_type axis, double sigma, int Df, int Dt, bool two_pass)
{
    check_params(Df, Dt, axis, nfreq, nt, stride, sigma);

    std_dev_clipper_buffers<float> buf;
    allocate_buffers(buf, nfreq, nt, axis, Df, Dt, two_pass);

    auto kernel = global_std_dev_clipper_kernel_table.get_kernel(axis, Df, Dt, two_pass);
    kernel(buf, intensity, weights, nfreq, nt, stride, sigma);
}


}  // namespace rf_pipelines
