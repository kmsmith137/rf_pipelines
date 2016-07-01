#include "rf_pipelines_internals.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif


wraparound_buf::wraparound_buf() :
    nfreq(0), nt_contig(0), nt_ring(0), ipos(0)
{ }


wraparound_buf::wraparound_buf(int nfreq_, int nt_contig_, int nt_ring_) :
    wraparound_buf()
{
    this->construct(nfreq_, nt_contig_, nt_ring_);
}


void wraparound_buf::construct(int nfreq_, int nt_contig_, int nt_ring_)
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


void wraparound_buf::setup_write(int it0, int nt, float* &intensityp, float* &weightp, int &stride)
{
    if ((nt <= 0) || (nt > nt_contig))
	throw runtime_error("wraparound_buf::setup_write(): invalid value of nt");
    if ((it0 < 0) || (it0 < ipos-nt_ring) || (it0 + nt > ipos))
	throw runtime_error("wraparound_buf::setup_write(): invalid value of it0");

    intensityp = &intensity[it0 % nt_ring];
    weightp = &weights[it0 % nt_ring];
    stride = nt_tot;
}


void wraparound_buf::finalize_write(int it0, int nt)
{
    if ((nt <= 0) || (nt > nt_contig))
	throw runtime_error("wraparound_buf::setup_write(): invalid value of nt");
    if ((it0 < 0) || (it0 < ipos-nt_ring) || (it0 + nt > ipos))
	throw runtime_error("wraparound_buf::setup_write(): invalid value of it0");

    it0 %= nt_ring;
    int it1 = it0 + nt;

    if (it0 < nt_contig)
	this->_copy(it0 + nt_ring, it0, min(nt_contig,it1) - it0);
    else if (it1 > nt_ring)
	this->_copy(0, nt_ring, it1-nt_ring);
}


void wraparound_buf::setup_append(int nt, float* &intensityp, float* &weightp, int &stride, bool zero_flag)
{
    if ((nt <= 0) || (nt > nt_contig))
	throw runtime_error("wraparound_buf::setup_append(): invalid value of nt");

    this->ipos += nt;
    this->setup_write(ipos-nt, nt, intensityp, weightp, stride);

    if (!zero_flag)
	return;

    for (int ifreq = 0; ifreq < nfreq; ifreq++) {
	memset(intensityp + ifreq*nt_tot, 0, nt * sizeof(float));
	memset(weightp + ifreq*nt_tot, 0, nt * sizeof(float));
    }
}


void wraparound_buf::finalize_append(int nt)
{
    if ((nt <= 0) || (nt > nt_contig) || (nt > ipos))
	throw runtime_error("wraparound_buf::finalize_append(): invalid value of nt");

    this->finalize_write(ipos-nt, nt);
}


void wraparound_buf::append_zeros(int nt)
{
    if ((nt <= 0) || (nt > nt_contig))
	throw runtime_error("wraparound_buf::append_zeros(): invalid value of nt");
    
    float *dummy_intensity;
    float *dummy_weights;
    int dummy_stride;
    bool zero_flag = true;

    this->setup_append(nt, dummy_intensity, dummy_weights, dummy_stride, zero_flag);
    this->finalize_append(nt);
}


void wraparound_buf::_check_integrity()
{
    for (int ifreq = 0; ifreq < nfreq; ifreq++) {
	int s1 = ifreq*nt_tot;
	int s2 = ifreq*nt_tot + nt_ring;

	rf_assert(!memcmp(&intensity[s1], &intensity[s2], nt_contig * sizeof(float)));
	rf_assert(!memcmp(&weights[s1], &weights[s2], nt_contig * sizeof(float)));
    }
}


void wraparound_buf::_copy(int it0_dst, int it0_src, int nt)
{
    int it1_dst = it0_dst + nt;
    int it1_src = it0_src + nt;

    // Argument checking
    rf_assert(nt > 0);
    rf_assert(it0_dst >= 0 && it1_dst <= nt_tot);
    rf_assert(it0_src >= 0 && it1_src <= nt_tot);
    rf_assert(max(it0_src,it0_dst) >= min(it1_src,it1_dst));
    
    for (int ifreq = 0; ifreq < nfreq; ifreq++) {
	memcpy(&intensity[ifreq*nt_tot + it0_dst], &intensity[ifreq*nt_tot + it0_src], nt * sizeof(float));
	memcpy(&weights[ifreq*nt_tot + it0_dst], &weights[ifreq*nt_tot + it0_src], nt * sizeof(float));
    }
}

// static member function
void wraparound_buf::run_unit_tests()
{
    cerr << "wraparound_buf::run_unit_tests()";

    for (int i = 0; i < 1000; i++) {
	if (i % 10 == 0)
	    cerr << ".";

	int nfreq = randint(1, 10);
	int nt_contig = randint(1, 10);
	int nt_ring = randint(1, 20);
	int nt_linear = randint(1000, 2000);

	vector<float> linear_ibuf(nfreq * nt_linear, 0.0);
	vector<float> linear_wbuf(nfreq * nt_linear, 0.0);

	wraparound_buf wrap_buf(nfreq, nt_contig, nt_ring);
	int ipos = 0;

	for (;;) {
	    float *intensity;
	    float *weights;
	    int stride;
	    
	    // Appending write

	    int nt_append = randint(1, nt_contig+1);
	    if (ipos + nt_append > nt_linear)
		break;

	    bool zero_flag = true;
	    wrap_buf.setup_append(nt_append, intensity, weights, stride, zero_flag);

	    for (int ifreq = 0; ifreq < nfreq; ifreq++) {
		for (int it = 0; it < nt_append; it++) {
		    int swrap = ifreq*stride + it;
		    int slin = ifreq*nt_linear + ipos + it;

		    intensity[swrap] = linear_ibuf[slin] = uniform_rand();
		    weights[swrap] = linear_wbuf[slin] = uniform_rand();
		}
	    }
	
	    wrap_buf.finalize_append(nt_append);
	    ipos += nt_append;

	    // Non-appending writes

	    for (int j = 0; j < 5; j++) {
		int it0_min = max(ipos - nt_contig, 0);
		int it0_write = randint(it0_min, ipos);
		int it1_write = randint(it0_write+1, ipos+1);
		int nt_write = it1_write - it0_write;
	    
		wrap_buf.setup_write(it0_write, nt_write, intensity, weights, stride);
	    
		for (int ifreq = 0; ifreq < nfreq; ifreq++) {
		    for (int it = 0; it < nt_write; it++) {
			int swrap = ifreq*stride + it;
			int slin = ifreq*nt_linear + it0_write + it;

			rf_assert(intensity[swrap] == linear_ibuf[slin]);
			rf_assert(weights[swrap] == linear_wbuf[slin]);

			intensity[swrap] = linear_ibuf[slin] = uniform_rand();
			weights[swrap] = linear_wbuf[slin] = uniform_rand();
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
