// Note: I haven't systematically documented the C++ interface to rf_pipelines,
// so the level of documentation will be hit-or-miss.  Also please note that the
// python-wrapping in rf_pipelines_c.cpp is kind of a mess which I hope to improve
// soon.  In the meantime if you want to python-wrap a C++ class, just email me
// and I'll help navigate the mess!

#include "rf_pipelines_internals.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif


wraparound_buf::wraparound_buf() :
    nfreq(0), nt_contig(0), nt_ring(0), ipos(0)
{ }


wraparound_buf::wraparound_buf(ssize_t nfreq_, ssize_t nt_contig_, ssize_t nt_ring_) :
    wraparound_buf()
{
    this->construct(nfreq_, nt_contig_, nt_ring_);
}


void wraparound_buf::construct(ssize_t nfreq_, ssize_t nt_contig_, ssize_t nt_ring_)
{
    if (this->nfreq != 0)
	throw runtime_error("double call to wraparound_buf::construct()");

    if (nfreq_ <= 0)
	throw runtime_error("wraparound_buf::construct(): invalid nfreq");
    if (nt_contig_ <= 0)
	throw runtime_error("wraparound_buf::construct(): invalid nt_contig");
    if (nt_ring_ <= 0)
	throw runtime_error("wraparound_buf::construct(): invalid nt_ring");

    this->nfreq = nfreq_;
    this->nt_contig = nt_contig_;

    // The property nt_ring >= 2*nt_contig is assumed in a few places
    this->nt_ring = max(nt_ring_, 2*nt_contig_);
    
    this->nt_tot = nt_ring + nt_contig;
    this->intensity.resize(nfreq * nt_tot, 0.0);
    this->weights.resize(nfreq * nt_tot, 0.0);
    this->ipos = 0;
}


void wraparound_buf::reset()
{
    this->nfreq = 0;
    this->nt_contig = 0;
    this->nt_ring = 0;
    this->nt_tot = 0;
    this->ipos = 0;

    deallocate(this->intensity);
    deallocate(this->weights);
}


void wraparound_buf::setup_write(ssize_t it0, ssize_t nt, float* &intensityp, float* &weightp, ssize_t &stride)
{
    if ((nt <= 0) || (nt > nt_contig))
	throw runtime_error("wraparound_buf::setup_write(): invalid value of nt");
    if ((it0 < 0) || (it0 < ipos-nt_ring) || (it0 + nt > ipos))
	throw runtime_error("wraparound_buf::setup_write(): invalid value of it0");

    intensityp = &intensity[it0 % nt_ring];
    weightp = &weights[it0 % nt_ring];
    stride = nt_tot;
}


void wraparound_buf::finalize_write(ssize_t it0, ssize_t nt)
{
    if ((nt <= 0) || (nt > nt_contig))
	throw runtime_error("wraparound_buf::setup_write(): invalid value of nt");
    if ((it0 < 0) || (it0 < ipos-nt_ring) || (it0 + nt > ipos))
	throw runtime_error("wraparound_buf::setup_write(): invalid value of it0");

    it0 %= nt_ring;
    ssize_t it1 = it0 + nt;

    if (it0 < nt_contig)
	this->_copy(it0 + nt_ring, it0, min(nt_contig,it1) - it0);
    else if (it1 > nt_ring)
	this->_copy(0, nt_ring, it1-nt_ring);
}


void wraparound_buf::setup_append(ssize_t nt, float* &intensityp, float* &weightp, ssize_t &stride, bool zero_flag)
{
    if ((nt <= 0) || (nt > nt_contig))
	throw runtime_error("wraparound_buf::setup_append(): invalid value of nt");

    this->ipos += nt;
    this->setup_write(ipos-nt, nt, intensityp, weightp, stride);

    if (!zero_flag)
	return;

    for (ssize_t ifreq = 0; ifreq < nfreq; ifreq++) {
	memset(intensityp + ifreq*nt_tot, 0, nt * sizeof(float));
	memset(weightp + ifreq*nt_tot, 0, nt * sizeof(float));
    }
}


void wraparound_buf::finalize_append(ssize_t nt)
{
    if ((nt <= 0) || (nt > nt_contig) || (nt > ipos))
	throw runtime_error("wraparound_buf::finalize_append(): invalid value of nt");

    this->finalize_write(ipos-nt, nt);
}


void wraparound_buf::append_zeros(ssize_t nt)
{
    if ((nt <= 0) || (nt > nt_contig))
	throw runtime_error("wraparound_buf::append_zeros(): invalid value of nt");
    
    float *dummy_intensity;
    float *dummy_weights;
    ssize_t dummy_stride;
    bool zero_flag = true;

    this->setup_append(nt, dummy_intensity, dummy_weights, dummy_stride, zero_flag);
    this->finalize_append(nt);
}


void wraparound_buf::_check_integrity()
{
    for (ssize_t ifreq = 0; ifreq < nfreq; ifreq++) {
	ssize_t s1 = ifreq*nt_tot;
	ssize_t s2 = ifreq*nt_tot + nt_ring;

	rf_assert(!memcmp(&intensity[s1], &intensity[s2], nt_contig * sizeof(float)));
	rf_assert(!memcmp(&weights[s1], &weights[s2], nt_contig * sizeof(float)));
    }
}


void wraparound_buf::_copy(ssize_t it0_dst, ssize_t it0_src, ssize_t nt)
{
    ssize_t it1_dst = it0_dst + nt;
    ssize_t it1_src = it0_src + nt;

    // Argument checking
    rf_assert(nt > 0);
    rf_assert(it0_dst >= 0 && it1_dst <= nt_tot);
    rf_assert(it0_src >= 0 && it1_src <= nt_tot);
    rf_assert(max(it0_src,it0_dst) >= min(it1_src,it1_dst));
    
    for (ssize_t ifreq = 0; ifreq < nfreq; ifreq++) {
	memcpy(&intensity[ifreq*nt_tot + it0_dst], &intensity[ifreq*nt_tot + it0_src], nt * sizeof(float));
	memcpy(&weights[ifreq*nt_tot + it0_dst], &weights[ifreq*nt_tot + it0_src], nt * sizeof(float));
    }
}

// static member function
void wraparound_buf::run_unit_tests(std::mt19937 &rng)
{
    cerr << "wraparound_buf::run_unit_tests()";

    for (int i = 0; i < 1000; i++) {
	if (i % 10 == 0)
	    cerr << ".";

	ssize_t nfreq = randint(rng, 1, 10);
	ssize_t nt_contig = randint(rng, 1, 10);
	ssize_t nt_ring = randint(rng, 1, 20);
	ssize_t nt_linear = randint(rng, 1000, 2000);

	vector<float> linear_ibuf(nfreq * nt_linear, 0.0);
	vector<float> linear_wbuf(nfreq * nt_linear, 0.0);

	wraparound_buf wrap_buf(nfreq, nt_contig, nt_ring);
	ssize_t ipos = 0;

	for (;;) {
	    float *intensity;
	    float *weights;
	    ssize_t stride;
	    
	    // Appending write

	    ssize_t nt_append = randint(rng, 1, nt_contig+1);
	    if (ipos + nt_append > nt_linear)
		break;

	    bool zero_flag = true;
	    wrap_buf.setup_append(nt_append, intensity, weights, stride, zero_flag);

	    for (ssize_t ifreq = 0; ifreq < nfreq; ifreq++) {
		for (ssize_t it = 0; it < nt_append; it++) {
		    ssize_t swrap = ifreq*stride + it;
		    ssize_t slin = ifreq*nt_linear + ipos + it;

		    intensity[swrap] = linear_ibuf[slin] = uniform_rand(rng);
		    weights[swrap] = linear_wbuf[slin] = uniform_rand(rng);
		}
	    }
	
	    wrap_buf.finalize_append(nt_append);
	    ipos += nt_append;

	    // Non-appending writes

	    for (int j = 0; j < 5; j++) {
		ssize_t it0_min = max(ipos - nt_contig, (ssize_t)0);
		ssize_t it0_write = randint(rng, it0_min, ipos);
		ssize_t it1_write = randint(rng, it0_write+1, ipos+1);
		ssize_t nt_write = it1_write - it0_write;
	    
		wrap_buf.setup_write(it0_write, nt_write, intensity, weights, stride);
	    
		for (ssize_t ifreq = 0; ifreq < nfreq; ifreq++) {
		    for (ssize_t it = 0; it < nt_write; it++) {
			ssize_t swrap = ifreq*stride + it;
			ssize_t slin = ifreq*nt_linear + it0_write + it;

			rf_assert(intensity[swrap] == linear_ibuf[slin]);
			rf_assert(weights[swrap] == linear_wbuf[slin]);

			intensity[swrap] = linear_ibuf[slin] = uniform_rand(rng);
			weights[swrap] = linear_wbuf[slin] = uniform_rand(rng);
		    }
		}

		wrap_buf.finalize_write(it0_write, nt_write);
	    }

	    wrap_buf._check_integrity();
	}
    }

    cerr << "done\n";
}


}  // namespace rf_pipelines
