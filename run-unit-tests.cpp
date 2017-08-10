// Note: I haven't systematically documented the C++ interface to rf_pipelines,
// so the level of documentation will be hit-or-miss.  Also please note that the
// python-wrapping in rf_pipelines_c.cpp is kind of a mess which I hope to improve
// soon.  In the meantime if you want to python-wrap a C++ class, just email me
// and I'll help navigate the mess!

#include "rf_pipelines_internals.hpp"

using namespace std;
using namespace rf_pipelines;


struct affine_map2 {
    double af;
    double at;
    double b;

    inline double apply(int ifreq, int it)
    {
	return af*ifreq + at*it + b;
    }

    static affine_map2 make_random(std::mt19937 &rng)
    {
	affine_map2 ret;
	ret.af = uniform_rand(rng, 0.1, 10.);
	ret.at = uniform_rand(rng, 0.01, 1.0);
	ret.b = uniform_rand(rng, 1., 10.);
	return ret;
    }
};


struct affine_map1 {
    double a;
    double b;

    // default constructor gives the identity
    affine_map1() : a(1.0), b(0.0) { }
    
    inline double apply(double x)
    {
	return a*x + b;
    }

    affine_map1 compose(const affine_map1 &m)
    {
	affine_map1 ret;
	ret.a = a*m.a;
	ret.b = a*m.b + b;
	return ret;
    }

    static affine_map1 make_random(std::mt19937 &rng)
    {
	affine_map1 ret;
	ret.a = uniform_rand(rng, 0.5, 1.0);
	ret.b = uniform_rand(rng, 0.0, 10.0);
	return ret;
    }
};


//
// The test_wi_stream simulates intensities and weights of the form
//   (si_f*ifreq + si_t*it + si0)    [ intensities ]
//   (sw_f*ifreq + sw_t*it + sw0)    [ weights ]
//
struct test_wi_stream : public wi_stream {
    std::mt19937 &rng;
    ssize_t nt_stream;
    affine_map2 intensity_map;
    affine_map2 weight_map;

    test_wi_stream(std::mt19937 &rng_) :
	rng(rng_)
    {
	this->nfreq = randint(rng, 1, 9);
	this->nt_maxwrite = randint(rng, 10, 21);
	this->nt_stream = randint(rng, 200, 401);
	this->intensity_map = affine_map2::make_random(rng);
	this->weight_map = affine_map2::make_random(rng);

	// arbitrary
	this->freq_lo_MHz = 400.;
	this->freq_hi_MHz = 800.;
	this->dt_sample = 1.0e-3;
    }

    virtual void stream_body(wi_run_state &run_state)
    {
	// arbitrary
	double t0 = 0.0;
	run_state.start_substream(t0);

	ssize_t ipos = 0;

	while (ipos < nt_stream) {
	    ssize_t nt = randint(rng, 1, nt_maxwrite+1);
	    nt = min(nt, nt_stream - ipos);

	    float *intensity;
	    float *weights;
	    ssize_t stride;
	    
	    bool zero_flag = true;
	    run_state.setup_write(nt, intensity, weights, stride, zero_flag);
	    
	    for (ssize_t ifreq = 0; ifreq < nfreq; ifreq++) {
		for (ssize_t it = 0; it < nt; it++) {
		    intensity[ifreq*stride + it] = intensity_map.apply(ifreq, ipos+it);
		    weights[ifreq*stride + it] = weight_map.apply(ifreq, ipos+it);
		}
	    }

	    run_state.finalize_write(nt);
	    ipos += nt;
	}

	run_state.end_substream();
    }
};


struct test_wi_transform : public wi_transform {
    affine_map1 my_imap, my_wmap;
    affine_map1 in_imap, in_wmap;
    affine_map1 out_imap, out_wmap;
    affine_map2 stream_imap, stream_wmap;
    
    double t0_substream;
    double dt_sample;

    ssize_t nt_stream;
    ssize_t curr_it;

 
    test_wi_transform(std::mt19937 &rng, const test_wi_stream &stream, const std::shared_ptr<test_wi_transform> &prev_transform)
    {
	// initialize fields in base class
	this->name = "test_transform";
	this->nfreq = stream.nfreq;
	this->nt_chunk = randint(rng, 1, 21);
	this->nt_prepad = max(randint(rng,-15,21), (ssize_t)0);    // order-one probability of zero
	this->nt_postpad = max(randint(rng,-15,21), (ssize_t)0);   // order-one probability of zero

	this->my_imap = affine_map1::make_random(rng);
	this->my_wmap = affine_map1::make_random(rng);

	if (prev_transform) {
	    this->in_imap = prev_transform->out_imap;
	    this->in_wmap = prev_transform->out_wmap;
	}

	this->out_imap = my_imap.compose(in_imap);
	this->out_wmap = my_wmap.compose(in_wmap);

	this->stream_imap = stream.intensity_map;
	this->stream_wmap = stream.weight_map;

	this->t0_substream = 0.0;
	this->dt_sample = stream.dt_sample;
	this->nt_stream = stream.nt_stream;
	this->curr_it = 0;
    }

    virtual void set_stream(const wi_stream &stream) override { return; }
    virtual void start_substream(int isubstream, double t0) override { this->t0_substream = t0; }
    virtual void end_substream() override { return; }

    virtual void process_chunk(double t0, double t1, float *intensity, float *weights, ssize_t stride, float *pp_intensity, float *pp_weights, ssize_t pp_stride) override
    {
	double t0_expected = t0_substream + curr_it * dt_sample;
	double t1_expected = t0_substream + (curr_it + nt_chunk) * dt_sample;
	rf_assert(fabs(t0 - t0_expected) < 1.0e-3 * dt_sample);
	rf_assert(fabs(t1 - t1_expected) < 1.0e-3 * dt_sample);	

	//
	// Check chunk + postpadded region
	//
	for (ssize_t ifreq = 0; ifreq < nfreq; ifreq++) {
	    for (ssize_t it = 0; it < nt_chunk+nt_postpad; it++) {
		ssize_t it2 = curr_it + it;
		double s_int = (it2 < nt_stream) ? stream_imap.apply(ifreq,it2) : 0.0;
		double s_wt = (it2 < nt_stream) ? stream_wmap.apply(ifreq,it2) : 0.0;

		rf_assert(reldist(intensity[ifreq*stride+it], in_imap.apply(s_int)) < 1.0e-5);
		rf_assert(reldist(weights[ifreq*stride+it], in_wmap.apply(s_wt)) < 1.0e-5);
	    }
	}

	//
	// Check prepadded region
	//
	for (ssize_t ifreq = 0; ifreq < nfreq; ifreq++) {
	    for (ssize_t it = 0; it < nt_prepad; it++) {
		ssize_t it2 = curr_it - nt_prepad + it;
		double expected_intensity = 0.0;
		double expected_weight = 0.0;

		if (it2 >= 0) {
		    double s_int = (it2 < nt_stream) ? stream_imap.apply(ifreq,it2) : 0.0;
		    double s_wt = (it2 < nt_stream) ? stream_wmap.apply(ifreq,it2) : 0.0;
		    
		    expected_intensity = in_imap.apply(s_int);
		    expected_weight = in_wmap.apply(s_wt);
		}

		rf_assert(reldist(pp_intensity[ifreq*pp_stride+it], expected_intensity) < 1.0e-5);
		rf_assert(reldist(pp_weights[ifreq*pp_stride+it], expected_weight) < 1.0e-5);
	    }
	}
	
	// apply transform
	for (ssize_t ifreq = 0; ifreq < nfreq; ifreq++) {
	    for (ssize_t it = 0; it < nt_chunk; it++) {
		intensity[ifreq*stride+it] = my_imap.apply(intensity[ifreq*stride+it]);
		weights[ifreq*stride+it] = my_wmap.apply(weights[ifreq*stride+it]);
	    }
	}

	this->curr_it += nt_chunk;
    }
};


static void run_pipeline_unit_tests(std::mt19937 &rng)
{
    cerr << "run_pipeline_unit_tests()";

    for (int iouter = 0; iouter < 1000; iouter++) {
	if (iouter % 10 == 0)
	    cerr << ".";
	
	// make random test stream
	test_wi_stream stream(rng);
	
	// make random transforms
	int ntransforms = randint(rng, 1, 10);

	vector<shared_ptr<wi_transform> > transforms(ntransforms);
	shared_ptr<test_wi_transform> prev;

	for (int itr = 0; itr < ntransforms; itr++)
	    transforms[itr] = prev = make_shared<test_wi_transform> (rng, stream, prev);

	int verbosity=0;
	stream.run(transforms, ".", nullptr, verbosity);
    }

    cerr << "done\n";
}


// -------------------------------------------------------------------------------------------------
//
// test_make_bitmask(): verifies that make_bitmask_reference() and make_bitmask() are equivalent


static void test_make_bitmask(std::mt19937 &rng, int nfreq, int nt, int in_stride)
{
    vector<float> in_weights(nfreq*in_stride, 0.0);
    vector<uint8_t> out_reference(nfreq*(nt/8), 0);
    vector<uint8_t> out_fast(nfreq*(nt/8), 0);

    for (unsigned int i = 0; i < in_weights.size(); i++)
	in_weights[i] = (uniform_rand(rng) < 0.25) ? 0.0 : uniform_rand(rng, -1.0, 2.0);

    make_bitmask_reference(&out_reference[0], nfreq, nt, &in_weights[0], in_stride);
    make_bitmask(&out_fast[0], nfreq, nt, &in_weights[0], in_stride);

    for (unsigned int i = 0; i < out_reference.size(); i++) {
	if (out_reference[i] != out_fast[i])
	    throw runtime_error("test_make_bitmask failed: nfreq=" + to_string(nfreq) + ", nt=" + to_string(nt) + ", in_stride=" + to_string(in_stride));
    }
}


static void test_make_bitmask(std::mt19937 &rng)
{
#ifdef __AVX__
    cerr << "test_make_bitmask()";

    for (int iouter = 0; iouter < 1000; iouter++) {
	if (iouter % 10 == 0)
	    cerr << ".";

	// Note: make_bitmask() requires nt to be a multiple of 256.
	int nfreq = randint(rng,1,100);
	int nt = randint(rng,1,6) * 256;
	int stride = randint(rng,nt,2*nt);

	test_make_bitmask(rng, nfreq, nt, stride);
    }

    cout << "done\n";
#else
    cout << "test_make_bitmask(): skipped on this machine (requires AVX)\n";
#endif
}


// -------------------------------------------------------------------------------------------------


// namespace rf_pipelines { extern void run_online_mask_filler_unit_tests(); }


int main(int argc, char **argv)
{
    std::random_device rd;
    std::mt19937 rng(rd());

    // run_online_mask_filler_unit_tests();
    wraparound_buf::run_unit_tests(rng);
    run_pipeline_unit_tests(rng);
    test_make_bitmask(rng);

    return 0;
}
