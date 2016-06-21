#ifndef _RF_PIPELINES_HPP
#define _RF_PIPELINES_HPP

#if (__cplusplus < 201103) && !defined(__GXX_EXPERIMENTAL_CXX0X__)
#error "This source file needs to be compiled with C++0x support (g++ -std=c++0x)"
#endif

#include <vector>

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif


struct noncopyable
{
    noncopyable() { }
    noncopyable(const noncopyable &) = delete;
    noncopyable& operator=(const noncopyable &) = delete;
};


// A helper class
struct wraparound_buf : noncopyable {
    // specified at construction
    int nfreq;
    int nt_contig;
    int nt_logical_size;

    // 2d arrays of shape (nfreq, nt_tot)
    std::vector<float> intensity;
    std::vector<float> weights;
    int nt_tot;

    int ipos;

    // Main constructor syntax
    wraparound_buf(int nfreq, int nt_contig, int nt_logical_size);

    // Alternate syntax: use default constuctor, then call construct()
    wraparound_buf();
    void construct(int nfreq, int nt_contig, int nt_logical_size);

    void setup_read(int it0, int nt, const float* &intensityp, const float* &weightp, int &stride) const;
    void setup_write(int it0, int nt, float* &intensityp, float* &weightp, int &stride, bool zero_flag);
    void finalize_write(int it0, int nt);

    void _copy(int it_dst, int it_src, int nt);
    void _check_integrity();

    static void run_unit_tests();
};


}  // namespace rf_pipelines

#endif // _RF_PIPELINES_HPP
