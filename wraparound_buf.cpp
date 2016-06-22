#include "rf_pipelines_internals.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif


wraparound_buf::wraparound_buf() :
    nfreq(0), nt_contig(0), nt_logical_size(0), ipos(0)
{ }


wraparound_buf::wraparound_buf(int nfreq_, int nt_contig_, int nt_logical_size_) :
    wraparound_buf()
{
    this->construct(nfreq_, nt_contig_, nt_logical_size_);
}


void wraparound_buf::construct(int nfreq_, int nt_contig_, int nt_logical_size_)
{
    if (this->nfreq != 0)
	throw runtime_error("double call to wraparound_buf::construct()");

    if (nfreq_ <= 0)
	throw runtime_error("wraparound_buf::construct(): invalid nfreq");
    if (nt_contig_ <= 0)
	throw runtime_error("wraparound_buf::construct(): invalid nt_contig");
    if (nt_logical_size_ <= 0)
	throw runtime_error("wraparound_buf::construct(): invalid nt_logical_size");

    this->nfreq = nfreq_;
    this->nt_contig = nt_contig_;

    // The property nt_logical_size >= 2*nt_contig is assumed in a few places
    this->nt_logical_size = max(nt_logical_size_, 2*nt_contig_);
    
    this->nt_tot = nt_logical_size + nt_contig;
    this->intensity.resize(nfreq * nt_tot, 0.0);
    this->weights.resize(nfreq * nt_tot, 0.0);
    this->ipos = 0;
}


void wraparound_buf::reset()
{
    this->nfreq = 0;
    this->nt_contig = 0;
    this->nt_logical_size = 0;
    this->nt_tot = 0;
    this->ipos = 0;

    deallocate(this->intensity);
    deallocate(this->weights);
}


void wraparound_buf::setup_read(int it0, int nt, float* &intensityp, float* &weightp, int &stride)
{
    if ((nt <= 0) || (nt > this->nt_contig))
	throw runtime_error("wraparound_buf::setup_read(): invalid value of nt");
    if ((it0 < 0) || (it0 < ipos-nt_logical_size) || (it0+nt > ipos))
	throw runtime_error("wraparound_buf::setup_read(): invalid value of it0");

    intensityp = &intensity[it0 % nt_logical_size];
    weightp = &weights[it0 % nt_logical_size];
    stride = nt_tot;
}


void wraparound_buf::setup_write(int it0, int nt, float* &intensityp, float* &weightp, int &stride, bool zero_flag)
{
    if ((nt <= 0) || (nt > nt_contig))
	throw runtime_error("wraparound_buf::setup_write(): invalid value of nt");
    if (it0 != ipos)
	throw runtime_error("wraparound_buf::setup_write(): invalid value of it0");

    intensityp = &intensity[it0 % nt_logical_size];
    weightp = &weights[it0 % nt_logical_size];
    stride = nt_tot;

    if (!zero_flag)
	return;
	
    for (int ifreq = 0; ifreq < nfreq; ifreq++) {
	memset(intensityp + ifreq*nt_tot, 0, nt * sizeof(float));
	memset(weightp + ifreq*nt_tot, 0, nt * sizeof(float));
    }
}


void wraparound_buf::finalize_write(int it0, int nt)
{
    if ((nt <= 0) || (nt > nt_contig))
	throw runtime_error("wraparound_buf::setup_write(): invalid value of nt");
    if (it0 != ipos)
	throw runtime_error("wraparound_buf::setup_write(): invalid value of it0");

    it0 %= nt_logical_size;
    int it1 = it0 + nt;

    if (it0 < nt_contig)
	this->_copy(it0 + nt_logical_size, it0, min(nt_contig,it1) - it0);
    else if (it1 > nt_logical_size)
	this->_copy(0, nt_logical_size, it1-nt_logical_size);

    this->ipos += nt;
}


void wraparound_buf::_check_integrity()
{
    for (int ifreq = 0; ifreq < nfreq; ifreq++) {
	int s1 = ifreq*nt_tot;
	int s2 = ifreq*nt_tot + nt_logical_size;

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

    for (int i = 0; i < 100; i++) {
	cerr << ".";

	int nfreq = randint(1, 10);
	int nt_contig = randint(1, 10);
	int nt_logical_size = randint(1, 20);
	float a = uniform_rand(100, 200);
	float b = uniform_rand(10, 20);
	float c = uniform_rand(1, 2);

	wraparound_buf wbuf(nfreq, nt_contig, nt_logical_size);
	int ipos = 0;

	for (int j = 0; j < 1000; j++) {

	    // Write

	    int nt_write = randint(1, nt_contig+1);
	    float *intensity = nullptr;
	    float *weight = nullptr;
	    int stride = 0;

	    bool zero_flag = true;
	    wbuf.setup_write(ipos, nt_write, intensity, weight, stride, zero_flag);

	    for (int ifreq = 0; ifreq < nfreq; ifreq++) {
		for (int it = 0; it < nt_write; it++) {
		    intensity[ifreq*stride + it] = b*ifreq + c*(ipos+it);
		    weight[ifreq*stride + it] = a + b*ifreq + c*(ipos+it);
		}
	    }
	
	    wbuf.finalize_write(ipos, nt_write);
	    ipos += nt_write;

	    wbuf._check_integrity();

	    // Read

	    int it0_min = max(ipos - nt_contig, 0);
	    int it0_read = randint(it0_min, ipos);
	    int nt_read = randint(1, ipos-it0_read+1);
	    
	    float *intensity2 = nullptr;
	    float *weight2 = nullptr;
	    stride = 0;
	    
	    wbuf.setup_read(it0_read, nt_read, intensity2, weight2, stride);
	    
	    for (int ifreq = 0; ifreq < nfreq; ifreq++) {
		for (int it = 0; it < nt_read; it++) {
		    rf_assert_close(intensity2[ifreq*stride + it], b*ifreq + c*(it0_read+it), 1.0e-2);
		    rf_assert_close(weight2[ifreq*stride + it], a + b*ifreq + c*(it0_read+it), 1.0e-2);
		}
	    }
	}
    }

    cerr << "done\n";
}


}  // namespace rf_pipelines
