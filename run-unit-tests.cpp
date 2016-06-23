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

    static affine_map2 make_random()
    {
	affine_map2 ret;
	ret.af = uniform_rand(0.1, 10.);
	ret.at = uniform_rand(0.01, 1.0);
	ret.b = uniform_rand(1., 10.);
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

    static affine_map1 make_random()
    {
	affine_map1 ret;
	ret.a = uniform_rand(0.5, 1.0);
	ret.b = uniform_rand(0.0, 10.0);
	return ret;
    }
};


//
// The test_wi_stream simulates intensities and weights of the form
//   (si_f*ifreq + si_t*it + si0)    [ intensities ]
//   (sw_f*ifreq + sw_t*it + sw0)    [ weights ]
//
struct test_wi_stream : public wi_stream {
    int nt_stream;
    affine_map2 intensity_map;
    affine_map2 weight_map;

    test_wi_stream()
    {
	this->nfreq = randint(1, 9);
	this->nt_maxwrite = randint(10, 21);
	this->nt_stream = randint(200, 401);
	this->intensity_map = affine_map2::make_random();
	this->weight_map = affine_map2::make_random();

	// arbitrary
	this->freq_lo_MHz = 400.;
	this->freq_hi_MHz = 800.;
	this->dt_sample = 1.0e-3;
    }

    virtual void run_stream(wi_run_state &rstate)
    {
	rstate.start_stream();
	int ipos = 0;

	while (ipos < nt_stream) {
	    int nt = randint(1, nt_maxwrite+1);
	    nt = min(nt, nt_stream - ipos);

	    float *intensity;
	    float *weights;
	    int stride;
	    
	    bool zero_flag = true;
	    rstate.setup_write(nt, intensity, weights, stride, zero_flag);
	    
	    for (int ifreq = 0; ifreq < nfreq; ifreq++) {
		for (int it = 0; it < nt; it++) {
		    intensity[ifreq*stride + it] = intensity_map.apply(ifreq, ipos+it);
		    weights[ifreq*stride + it] = weight_map.apply(ifreq, ipos+it);
		}
	    }

	    rstate.finalize_write(nt);
	    ipos += nt;
	}

	rstate.end_stream();
    }
};


struct test_wi_transform : public wi_transform {
    affine_map1 my_imap, my_wmap;
    affine_map1 in_imap, in_wmap;
    affine_map1 out_imap, out_wmap;
    affine_map2 stream_imap, stream_wmap;

    int nt_stream;
    int curr_it;

 
    test_wi_transform(const test_wi_stream &stream, const std::shared_ptr<test_wi_transform> &prev_transform)
    {
	this->nfreq = stream.nfreq;
	this->nt_stream = stream.nt_stream;
	this->nt_chunk = randint(1, 21);
	this->nt_prepad = max(randint(-15,21), 0);    // order-one probability of zero
	this->nt_postpad = max(randint(-15,21), 0);   // order-one probability of zero

	this->my_imap = affine_map1::make_random();
	this->my_wmap = affine_map1::make_random();

	if (prev_transform) {
	    this->in_imap = prev_transform->out_imap;
	    this->in_wmap = prev_transform->out_wmap;
	}

	this->out_imap = my_imap.compose(in_imap);
	this->out_wmap = my_wmap.compose(in_wmap);

	this->stream_imap = stream.intensity_map;
	this->stream_wmap = stream.weight_map;
    }

    
    virtual void start_stream(int nfreq_, double freq_lo_MHz, double freq_hi_MHz, double dt_sample)
    {
	rf_assert(this->nfreq == nfreq_);
	this->curr_it = 0;
    }

    virtual void process_chunk(float *intensity, float *weight, int stride, float *pp_intensity, float *pp_weight, int pp_stride)
    {
	//
	// Check chunk + postpadded region
	//
	for (int ifreq = 0; ifreq < nfreq; ifreq++) {
	    for (int it = 0; it < nt_chunk+nt_postpad; it++) {
		int it2 = curr_it + it;
		double s_int = (it2 < nt_stream) ? stream_imap.apply(ifreq,it2) : 0.0;
		double s_wt = (it2 < nt_stream) ? stream_wmap.apply(ifreq,it2) : 0.0;

		rf_assert(reldist(intensity[ifreq*stride+it], in_imap.apply(s_int)) < 1.0e-5);
		rf_assert(reldist(weight[ifreq*stride+it], in_wmap.apply(s_wt)) < 1.0e-5);
	    }
	}

	//
	// Check prepadded region
	//
	for (int ifreq = 0; ifreq < nfreq; ifreq++) {
	    for (int it = 0; it < nt_prepad; it++) {
		int it2 = curr_it - nt_prepad + it;
		double expected_intensity = 0.0;
		double expected_weight = 0.0;

		if (it2 >= 0) {
		    double s_int = (it2 < nt_stream) ? stream_imap.apply(ifreq,it2) : 0.0;
		    double s_wt = (it2 < nt_stream) ? stream_wmap.apply(ifreq,it2) : 0.0;
		    
		    expected_intensity = in_imap.apply(s_int);
		    expected_weight = in_wmap.apply(s_wt);
		}

		rf_assert(reldist(pp_intensity[ifreq*pp_stride+it], expected_intensity) < 1.0e-5);
		rf_assert(reldist(pp_weight[ifreq*pp_stride+it], expected_weight) < 1.0e-5);
	    }
	}
	
	// apply transform
	for (int ifreq = 0; ifreq < nfreq; ifreq++) {
	    for (int it = 0; it < nt_chunk; it++) {
		intensity[ifreq*stride+it] = my_imap.apply(intensity[ifreq*stride+it]);
		weight[ifreq*stride+it] = my_wmap.apply(weight[ifreq*stride+it]);
	    }
	}

	this->curr_it += nt_chunk;
    }

    virtual void end_stream() { }
};


static void run_pipeline_unit_tests()
{
    cerr << "run_pipeline_unit_tests()";

    for (int iouter = 0; iouter < 1000; iouter++) {
	if (iouter % 10 == 0)
	    cerr << ".";
	
	// make random test stream
	test_wi_stream stream;
	
	// make random transforms
	int ntransforms = randint(1, 10);

	vector<shared_ptr<wi_transform> > transforms(ntransforms);
	shared_ptr<test_wi_transform> prev;

	for (int itr = 0; itr < ntransforms; itr++)
	    transforms[itr] = prev = make_shared<test_wi_transform> (stream, prev);

	stream.run_transforms(transforms);
    }

    cerr << "done\n";
}


// -------------------------------------------------------------------------------------------------


int main(int argc, char **argv)
{
    wraparound_buf::run_unit_tests();
    run_pipeline_unit_tests();

    return 0;
}
