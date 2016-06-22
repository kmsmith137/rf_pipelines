#include "rf_pipelines_internals.hpp"

using namespace std;
using namespace rf_pipelines;


struct affine_map {
    double cf;
    double ct;
    double c0;

    inline double eval(int ifreq, int it)
    {
	return cf*ifreq + ct*it + c0;
    }

    static affine_map make_random()
    {
	affine_map ret;
	ret.cf = uniform_rand(1., 100.);
	ret.ct = uniform_rand(0.01, 1.0);
	ret.c0 = uniform_rand(1., 200.);
	return ret;
    }
};


struct linear_map {
    double a;
    double b;
    
    inline double apply(double x)
    {
	return a*x + b;
    }
    
    affine_map apply(const affine_map &m)
    {
	affine_map ret;
	ret.cf = a * m.cf;
	ret.ct = a * m.ct;
	ret.c0 = a * m.c0 + b;
	return ret;
    }

    static linear_map make_random()
    {
	linear_map ret;
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
    int stream_len;
    affine_map intensity_map;
    affine_map weight_map;

    test_wi_stream()
    {
	this->nfreq = randint(1, 9);
	this->nt_maxwrite = randint(10, 21);
	this->stream_len = randint(200, 401);
	this->intensity_map = affine_map::make_random();
	this->weight_map = affine_map::make_random();

	// arbitrary
	this->freq_lo_MHz = 400.;
	this->freq_hi_MHz = 800.;
	this->dt_sample = 1.0e-3;
    }

    virtual void run_stream(wi_run_state &rstate)
    {
	rstate.start_stream();
	int ipos = 0;

	while (ipos < stream_len) {
	    int nt = randint(1, nt_maxwrite+1);
	    nt = min(nt, stream_len - ipos);

	    float *intensity;
	    float *weights;
	    int stride;
	    
	    bool zero_flag = true;
	    rstate.setup_write(nt, intensity, weights, stride, zero_flag);
	    
	    for (int ifreq = 0; ifreq < nfreq; ifreq++) {
		for (int it = 0; it < nt; it++) {
		    intensity[ifreq*stride + it] = intensity_map.eval(ifreq, ipos+it);
		    weights[ifreq*stride + it] = weight_map.eval(ifreq, ipos+it);
		}
	    }

	    rstate.finalize_write(nt);
	    ipos += nt;
	}
	    
	rstate.end_stream();
    }
};


struct test_wi_transform : public wi_transform {
    int nfreq;
    int curr_it;

    linear_map intensity_transform, weight_transform;
    affine_map input_intensity_map, input_weight_map;
    affine_map output_intensity_map, output_weight_map;

 
    test_wi_transform(const test_wi_stream &stream, const std::shared_ptr<test_wi_transform> &prev_transform)
    {
	this->nt_chunk = randint(1, 21);
	this->nt_prepad = max(randint(-15,21), 0);    // order-one probability of zero
	this->nt_postpad = max(randint(-15,21), 0);   // order-one probability of zero
	this->nfreq = 0;  // to be initialized later in start_stream()

	this->intensity_transform = linear_map::make_random();
	this->weight_transform = linear_map::make_random();
	
	if (prev_transform) {
	    this->input_intensity_map = prev_transform->output_intensity_map;
	    this->input_weight_map = prev_transform->output_weight_map;
	}
	else {
	    this->input_intensity_map = stream.intensity_map;
	    this->input_weight_map = stream.weight_map;
	}

	this->output_intensity_map = intensity_transform.apply(input_intensity_map);
	this->output_weight_map = weight_transform.apply(input_weight_map);
    }

    
    virtual void start_stream(int nfreq_, double freq_lo_MHz, double freq_hi_MHz, double dt_sample)
    {
	this->nfreq = nfreq_;
	this->curr_it = 0;
    }

    virtual void process_chunk(float *intensity, float *weight, int stride, float *pp_intensity, float *pp_weight, int pp_stride)
    {
	// check that data is what we expect
	for (int ifreq = 0; ifreq < nfreq; ifreq++) {
	    for (int it = 0; it < nt_chunk+nt_postpad; it++) {
		rf_assert_close(intensity[ifreq*stride+it], input_intensity_map.eval(ifreq,curr_it+it), 1.0e-5);
		rf_assert_close(weight[ifreq*stride+it], input_weight_map.eval(ifreq,curr_it+it), 1.0e-5);
	    }
	}

	// check that prepradded data is what we expect
	for (int ifreq = 0; ifreq < nfreq; ifreq++) {
	    for (int it = max(nt_prepad-curr_it,0); it < nt_prepad; it++) {
		rf_assert_close(pp_intensity[ifreq*pp_stride+it], input_intensity_map.eval(ifreq,curr_it-nt_prepad+it), 1.0e-5);
		rf_assert_close(pp_weight[ifreq*pp_stride+it], input_weight_map.eval(ifreq,curr_it-nt_prepad+it), 1.0e-5);
	    }
	}
	
	// apply transform
	for (int ifreq = 0; ifreq < nfreq; ifreq++) {
	    for (int it = 0; it < nt_chunk; it++) {
		intensity[ifreq*stride+it] = intensity_transform.apply(intensity[ifreq*stride+it]);
		weight[ifreq*stride+it] = weight_transform.apply(weight[ifreq*stride+it]);
	    }
	}

	this->curr_it += nt_chunk;
    }

    virtual void end_stream() { }
};


static void run_pipeline_unit_tests()
{
    cerr << "run_pipeline_unit_tests().";

    for (int iouter = 0; iouter < 100; iouter++) {
	cerr << ".";
	
	// make random test stream
	test_wi_stream stream;
	
	// make random transforms
	int ntransforms = randint(1, 5);

	vector<shared_ptr<wi_transform> > transforms(ntransforms);
	shared_ptr<test_wi_transform> prev;

	for (int itr = 0; itr < ntransforms; itr++)
	    transforms[itr] = prev = make_shared<test_wi_transform> (stream, prev);

	stream.run_transforms(transforms);

	// final check?
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
