#include "rf_pipelines_internals.hpp"

#include <simd_helpers/simd_float32.hpp>
#include <simd_helpers/downsample_bitwise_or.hpp>

using namespace std;
using namespace simd_helpers;

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif


// -------------------------------------------------------------------------------------------------
//
// make_bitmask() kernel: slow reference version.
//
// Input: float32 array with shape (nfreq,nt)
// Output: uint8 array with shape (nfreq,nt/8)
//
// Notes: nt must be a multiple of 8, and input array can have an arbitrary stride.


void make_bitmask_reference(uint8_t *out_bitmask, int nfreq, int nt, const float *in_weights, int in_stride)
{
    rf_assert(out_bitmask != nullptr);
    rf_assert(in_weights != nullptr);
    rf_assert(nfreq > 0);
    rf_assert(nt > 0);
    rf_assert(nt % 8 == 0);
    rf_assert(in_stride >= nt);

    for (int ifreq = 0; ifreq < nfreq; ifreq++) {
	for (int it = 0; it < nt; it += 8) {
	    uint8_t out = 0;

	    for (int j = 0; j < 8; j++) {
		float w = in_weights[ifreq*in_stride + it + j];
		uint8_t b = (w > 0.0) ? 1 : 0;
		out |= (b << j);
	    }

	    out_bitmask[(ifreq*nt+it)/8] = out;
	}
    }
}


// -------------------------------------------------------------------------------------------------
//
// make_bitmask() kernel: fast version.
//
// Input: float32 array with shape (nfreq,nt)
// Output: uint8 array with shape (nfreq,nt/8)
//
// Note: In the fast kernel, nt must be a multiple of 256(!).  


// Helper class for make_bitmask().
struct make_bitmask_helper {
    const simd_t<int,8> c0;
    const simd_t<int,8> c1;
    const simd_t<int,8> c2;
    const simd_t<int,8> c3;

    make_bitmask_helper() :
	c0(_mm256_set_epi32(1U, 1U<<1, 1U<<2, 1U<<3, 1U<<4, 1U<<5, 1U<<6, 1U<<7)),
	c1(_mm256_set_epi32(1U<<8, 1U<<9, 1U<<10, 1U<<11, 1U<<12, 1U<<13, 1U<<14, 1U<<15)),
	c2(_mm256_set_epi32(1U<<16, 1U<<17, 1U<<18, 1U<<19, 1U<<20, 1U<<21, 1U<<22, 1U<<23)),
	c3(_mm256_set_epi32(1U<<24, 1U<<25, 1U<<26, 1U<<27, 1U<<28, 1U<<29, 1U<<30, 1U<<31))
    { }

    inline simd_t<int,8> process8(const float *p, simd_t<int,8> c)
    {
	simd_t<float,8> w = simd_load<float,8> (p);
	simd_t<int,8> b = w.compare_gt(simd_t<float,8>::zero());
	return b.bitwise_and(c);
    }

    inline simd_t<int,8> process32(const float *p)
    {
	simd_t<int,8> b0 = process8(p, c0);
	simd_t<int,8> b1 = process8(p+8, c1);
	simd_t<int,8> b2 = process8(p+16, c2);
	simd_t<int,8> b3 = process8(p+24, c3);

	b0 = b0.bitwise_or(b1);
	b2 = b2.bitwise_or(b3);
	return b0.bitwise_or(b2);
    }

    inline void process256(int *out, const float *p)
    {
	simd_ntuple<int,8,8> t;
	t.template extract<0>() = process32(p);
	t.template extract<1>() = process32(p+32);
	t.template extract<2>() = process32(p+64);
	t.template extract<3>() = process32(p+96);
	t.template extract<4>() = process32(p+128);
	t.template extract<5>() = process32(p+160);
	t.template extract<6>() = process32(p+192);
	t.template extract<7>() = process32(p+224);

	simd_t<int,8> u = downsample_bitwise_or(t);
	simd_store(out, u);
    }
};


void make_bitmask(uint8_t *out_bitmask, int nfreq, int nt, const float *in_weights, int in_stride)
{
    rf_assert(out_bitmask != nullptr);
    rf_assert(in_weights != nullptr);
    rf_assert(nfreq > 0);
    rf_assert(nt > 0);
    rf_assert(in_stride >= nt);

    if ((nt % 256) != 0)
	throw runtime_error("make_bitmask: assembly-language-kernelized implementation assumes nt(=" + to_string(nt) + ") is a multiple of 256");

    make_bitmask_helper h;

    int *out = reinterpret_cast<int *> (out_bitmask);
    int n = nt / 256;

    for (int ifreq = 0; ifreq < nfreq; ifreq++) {
	for (int i = 0; i < n; i++)
	    h.process256(out + 8*(ifreq*n+i), in_weights + ifreq*in_stride + 256*i);
    }
}


// -------------------------------------------------------------------------------------------------
//
// chunk_guard: helper class used by the bitmask_saver, to ensure that calls to
// bitmask_chunk_manager::get_chunk() and bitmask_chunk_manager::put_chunk()
// are correctly "paired".


struct chunk_guard {
    bitmask_chunk_manager *mp;
    uint8_t *p;

    chunk_guard(bitmask_chunk_manager *mp, double t0, ssize_t nfreq, ssize_t nt_chunk);
    ~chunk_guard();

    // Noncopyable
    chunk_guard(const chunk_guard &) = delete;
    chunk_guard &operator=(const chunk_guard &) = delete;
};


chunk_guard::chunk_guard(bitmask_chunk_manager *mp_, double t0, ssize_t nfreq, ssize_t nt_chunk)
{
    this->mp = mp_;
    this->p = mp->get_chunk(t0, nfreq, nt_chunk);

    if (!p)
	throw runtime_error("rf_pipelines: bitmask_chunk_manager::get_chunk() returned null pointer");
}


chunk_guard::~chunk_guard()
{
    mp->put_chunk();
    mp = nullptr;
    p = nullptr;
}


// -------------------------------------------------------------------------------------------------
//
// bitmask_saver transform (just a wrapper around the make_bitmask() kernel above)
    

struct bitmask_saver : public wi_transform
{
    shared_ptr<bitmask_chunk_manager> mp;

    bitmask_saver(const shared_ptr<bitmask_chunk_manager> &mp, ssize_t nt_chunk);

    virtual void set_stream(const wi_stream &stream) override;
    virtual void start_substream(int isubstream, double t0) override { }
    virtual void process_chunk(double t0, double t1, float *intensity, float *weights, ssize_t stride, float *pp_intensity, float *pp_weights, ssize_t pp_stride) override;
    virtual void end_substream() override { }
};


bitmask_saver::bitmask_saver(const shared_ptr<bitmask_chunk_manager> &mp_, ssize_t nt_chunk_)
{
    this->mp = mp_;
    this->nt_chunk = nt_chunk_;

    rf_assert(mp);
    rf_assert(nt_chunk > 0);
    rf_assert(nt_chunk % 256 == 0);
}


void bitmask_saver::set_stream(const wi_stream &stream)
{
    this->nfreq = stream.nfreq;
}


void bitmask_saver::process_chunk(double t0, double t1, float *intensity, float *weights, ssize_t stride, float *pp_intensity, float *pp_weights, ssize_t pp_stride)
{
    // Use of chunk_guard here ensures that bitmask_chunk_manager::put_chunk()
    // always gets called, even if make_bitmask() throws an exception.
    chunk_guard chunk(this->mp.get(), t0, this->nfreq, this->nt_chunk);
    make_bitmask(chunk.p, this->nfreq, this->nt_chunk, weights, stride);
}


// -------------------------------------------------------------------------------------------------
//
// Externally visible factory function


shared_ptr<wi_transform> make_bitmask_saver(const shared_ptr<bitmask_chunk_manager> &mp, ssize_t nt_chunk)
{
    return make_shared<bitmask_saver> (mp, nt_chunk);
}


}  // namespace rf_pipelines
