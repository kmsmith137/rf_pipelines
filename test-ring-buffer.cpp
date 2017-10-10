#include "rf_pipelines_internals.hpp"

using namespace std;
using namespace rf_pipelines;


static vector<ssize_t> make_random_cdims(std::mt19937 &rng)
{
    vector<ssize_t> cdims;
    
    do {
	cdims.push_back(randint(rng,1,8));
    } while ((prod(cdims) < 20) && (cdims.size() <= 4) && (uniform_rand(rng) < 0.5));

    return cdims;
}


// Note: for testing a single ring buffer, outside the context of a larger pipeline, it's
// convenient to use downsampled indices throughout, and multiply by 'nds' when calling 
// ring_buffer methods which expect non-downsampled indices (e.g. ring_buffer::get()).

static void test_ring_buffer(std::mt19937 &rng, const vector<ssize_t> &cdims, ssize_t nds, ssize_t nt_contig, ssize_t nt_maxlag)
{
    shared_ptr<ring_buffer> rb = make_shared<ring_buffer> (cdims, nds);
    rb->update_params(nt_contig * nds, nt_maxlag * nds);
    rb->allocate();

    ssize_t buf_pos0 = 0;
    ssize_t buf_pos1 = 0;
    ssize_t csize = rb->csize;
    ssize_t stride = rb->get_stride();
    
    // Reference ring buffer
    vector<float> rb_ref(nt_maxlag * csize, 0.0);

    for (int iter = 0; iter < 10000; iter++) {

	// Randomly generate pos0, pos1, mode

	ssize_t pos0 = -1;
	ssize_t pos1 = -1;
	int mode = ring_buffer::ACCESS_NONE;

	double append_prob = nt_contig / double(nt_contig + 2*buf_pos1 - 2*buf_pos0);
	
	if (uniform_rand(rng) < append_prob) {
	    pos0 = buf_pos1;
	    pos1 = pos0 + randint(rng, 1, nt_contig+1);
	    mode = ring_buffer::ACCESS_APPEND;
	}
	else {
	    ssize_t nmax = min(buf_pos1 - buf_pos0, nt_contig);
	    ssize_t n = randint(rng, 1, nmax+1);
	    pos0 = randint(rng, buf_pos0, buf_pos1-n+1);
	    pos1 = pos0 + n;
	    mode = randint(rng, ring_buffer::ACCESS_READ, ring_buffer::ACCESS_RW+1);
	}
	
	float *p = rb->get(pos0 * nds, pos1 * nds, mode);

	if (mode & ring_buffer::ACCESS_READ) {
	    for (ssize_t it = 0; it < (pos1-pos0); it++) {
		float *r = &rb_ref[((it+pos0) % nt_maxlag) * csize];
		for (ssize_t ic = 0; ic < csize; ic++)
		    rf_assert(p[ic*stride+it] == r[ic]);
	    }
	}

	if (mode & ring_buffer::ACCESS_WRITE) {
	    for (ssize_t it = 0; it < (pos1-pos0); it++) {
		float *r = &rb_ref[((it+pos0) % nt_maxlag) * csize];
		for (ssize_t ic = 0; ic < csize; ic++)
		    p[ic*stride+it] = r[ic] = uniform_rand(rng);
	    }
	}

	rb->put(p, pos0 * nds, pos1 * nds, mode);

	buf_pos1 = max(buf_pos1, pos1);
	buf_pos0 = max(buf_pos0, buf_pos1 - nt_maxlag);
    }
}


int main(int argc, char **argv)
{
    const int nouter = 1000;

    std::random_device rd;
    std::mt19937 rng(rd());

    for (int iouter = 0; iouter < nouter; iouter++) {
	if (iouter % 50 == 0)
	    cout << "test-ring-buffer: iteration " << iouter << "/" << nouter << endl;

	vector<ssize_t> cdims = make_random_cdims(rng);
	ssize_t nt_contig = randint(rng, 1, 50);
	ssize_t nt_maxlag = randint(rng, nt_contig, 10 * nt_contig);
	ssize_t nds = randint(rng, 1, 5);

	test_ring_buffer(rng, cdims, nds, nt_contig, nt_maxlag);
    }

    cout << "test-ring-buffer: pass" << endl;
    return 0;
}
