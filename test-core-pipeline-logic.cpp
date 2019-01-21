#include "rf_pipelines_internals.hpp"

using namespace std;
using namespace rf_pipelines;


static string bufname_from_index(int ix)
{
    rf_assert(ix >= 0 && ix <= 100);
    return "BUFFER" + to_string(ix);
}


// -------------------------------------------------------------------------------------------------


struct reference_pipeline_object {
    ssize_t csize = 0;
    ssize_t nds = 0;

    virtual void apply_reference(vector<vector<float>> &v, ssize_t nt_end) = 0;
    virtual shared_ptr<pipeline_object> make_real_pipeline_object(ssize_t nt_end) = 0;

    void run_test(ssize_t nt_end);
};


struct reference_pipeline : reference_pipeline_object 
{
    vector<shared_ptr<reference_pipeline_object>> stages;

    void add(const shared_ptr<reference_pipeline_object> &p)
    {
	if (stages.size() == 0) {
	    this->csize = p->csize;
	    this->nds = p->nds;
	}

	rf_assert(p->csize == this->csize);
	rf_assert(p->nds == this->nds);
	stages.push_back(p);
    }

    virtual void apply_reference(vector<vector<float>> &v, ssize_t nt_end) override
    {
	for (const auto &p: stages)
	    p->apply_reference(v, nt_end);
    }

    virtual shared_ptr<pipeline_object> make_real_pipeline_object(ssize_t nt_end) override
    {
	shared_ptr<pipeline> ret = make_shared<pipeline> ();

	for (const auto &p: this->stages) {
	    auto rp = p->make_real_pipeline_object(nt_end);
	    ret->add(rp);
	}

	return ret;
    }
};


void reference_pipeline_object::run_test(ssize_t nt_end)
{
    // Not currently using 'nds' in the code (this is a FIXME, not a deliberate decision)
    rf_assert(this->nds == 1); 

    // Run reference_pipeline.
    vector<vector<float>> ref_buf;
    this->apply_reference(ref_buf, nt_end);

    // The output of the reference_pipeline is a list of buffers.  Some buffers are
    // absent (indicated by zero-length vectors), and the rest have shape (nt_end, csize).
    for (unsigned int ix = 0; ix < ref_buf.size(); ix++)
	rf_assert((ref_buf[ix].size() == nt_end*csize) || (ref_buf[ix].size() == 0));

    vector<string> spool_bufnames;
    for (unsigned int ix = 0; ix < ref_buf.size(); ix++)
	if (ref_buf[ix].size() > 0)
	    spool_bufnames.push_back(bufname_from_index(ix));

    auto p1 = this->make_real_pipeline_object(nt_end);
    auto p2 = make_shared<pipeline_spool> (spool_bufnames);

    auto p = make_shared<pipeline> ();
    p->add(p1);
    p->add(p2);

    run_params params;
    params.outdir = "";
    params.verbosity = 0;
    params.debug = true;
    p->run(params);

    for (unsigned int ix = 0; ix < ref_buf.size(); ix++) {
	if (ref_buf[ix].size() == 0)
	    continue;

	auto sp = p2->get_spooled_buffer(bufname_from_index(ix));
	rf_assert(sp->csize == this->csize);
	rf_assert(sp->nds == this->nds);
	rf_assert(sp->nt >= nt_end);   // note: can be padded past nt_end
	rf_assert(sp->data.size() == csize * sp->nt);

	for (ssize_t ic = 0; ic < csize; ic++) {
	    for (ssize_t it = 0; it < nt_end; it++) {
		float x = sp->data[ic*sp->nt + it];
		float y = ref_buf[ix][it*csize + ic];
		rf_assert(abs(x-y) < 1.0e-5);
	    }
	}
    }
}


// -------------------------------------------------------------------------------------------------


// Called in non-reference pipeline
static void _randomize(ring_buffer_subarray &arr, std::mt19937 &rng)
{
    ssize_t nc = arr.buf->csize;
    ssize_t nt = arr.pos1 - arr.pos0;
    ssize_t stride = arr.stride;
    float *data = arr.data;

    // Note loop ordering!
    for (ssize_t it = 0; it < nt; it++)
	for (ssize_t ic = 0; ic < nc; ic++)
	    data[ic*stride + it] = uniform_rand(rng, -1.0, 1.0);
}


static void _rotate(float *dst1, float *dst2, const float *src1, const float *src2, ssize_t n, float cos_theta, float sin_theta)
{
    for (ssize_t i = 0; i < n; i++) {
	float x = src1[i];
	float y = src2[i];

	dst1[i] = cos_theta*x + sin_theta*y;
	dst2[i] = -sin_theta*x + cos_theta*y;
    }
}


struct rot2_params {
    int nfreq = 0;
    int nt_chunk = 0;

    // Rotation to be applied
    float theta = 0.0;
    
    // Indices of input and output buffers.
    int ix_in0 = 0;
    int ix_in1 = 0;
    int ix_out0 = 0;
    int ix_out1 = 0;

    bool create_in0 = false;
    bool create_in1 = false;
    bool create_out0 = false;
    bool create_out1 = false;

    unsigned int seed_in0 = 1;   // only used if create_in0 == true
    unsigned int seed_in1 = 1;   // only used if create_in1 == true

    void validate() const
    {
	rf_assert(nfreq > 0);
	rf_assert(nt_chunk > 0);
	rf_assert(ix_in0 >= 0);
	rf_assert(ix_in1 >= 0);
	rf_assert(ix_in0 != ix_in1);
	rf_assert(ix_out0 >= 0);
	rf_assert(ix_out1 >= 0);
	rf_assert(ix_out0 != ix_out1);
    }

    bool can_be_first() const
    {
	return (create_in0 && create_in1);
    }

    static rot2_params make_random(std::mt19937 &rng, int nfreq_, vector<bool> &bflag)
    {
	rot2_params ret;
	
	int nbuf = bflag.size();
	rf_assert(nbuf >= 2);
	
	ret.nfreq = nfreq_;
	ret.nt_chunk = randint(rng, 1, 26);
	ret.theta = uniform_rand(rng, 0, 2*M_PI);
	ret.seed_in0 = randint(rng, 1, 100000);
	ret.seed_in1 = randint(rng, 1, 100000);

	do {
	    ret.ix_in0 = randint(rng, 0, nbuf);
	    ret.ix_in1 = randint(rng, 0, nbuf);
	} while (ret.ix_in0 == ret.ix_in1);
	
	do {
	    ret.ix_out0 = randint(rng, 0, nbuf);
	    ret.ix_out1 = randint(rng, 0, nbuf);
	} while (ret.ix_out0 == ret.ix_out1);

	ret.create_in0 = !bflag[ret.ix_in0];
	ret.create_in1 = !bflag[ret.ix_in1];
	bflag[ret.ix_in0] = bflag[ret.ix_in1] = true;

	ret.create_out0 = !bflag[ret.ix_out0];
	ret.create_out1 = !bflag[ret.ix_out1];
	bflag[ret.ix_out0] = bflag[ret.ix_out1] = true;

	return ret;
    }
};


struct rot2 : public chunked_pipeline_object {
    const rot2_params params;
    const ssize_t nt_end = 0;

    std::mt19937 rng_in0;
    std::mt19937 rng_in1;
    
    float cos_theta = 0.0;
    float sin_theta = 0.0;

    // Initialized in _bind()
    shared_ptr<ring_buffer> buf_in0;
    shared_ptr<ring_buffer> buf_in1;
    shared_ptr<ring_buffer> buf_out0;
    shared_ptr<ring_buffer> buf_out1;

    // Initialized in _bind()
    int mode_in0 = ring_buffer::ACCESS_NONE;
    int mode_in1 = ring_buffer::ACCESS_NONE;
    int mode_out0 = ring_buffer::ACCESS_NONE;
    int mode_out1 = ring_buffer::ACCESS_NONE;

    ring_buffer_subarray arr_in0;
    ring_buffer_subarray arr_in1;
    ring_buffer_subarray _arr_out0;
    ring_buffer_subarray _arr_out1;

    ring_buffer_subarray *arr_out0 = &_arr_out0;
    ring_buffer_subarray *arr_out1 = &_arr_out1;

    
    rot2(const rot2_params &params_, ssize_t nt_end_) :
	chunked_pipeline_object("rot2", params_.can_be_first()),
	params(params_),
	nt_end(nt_end_),
	rng_in0(params_.seed_in0),
	rng_in1(params_.seed_in1)
    {
	this->nt_chunk = params_.nt_chunk;

	params.validate();
	cos_theta = cos(params.theta);
	sin_theta = sin(params.theta);
    }

    
    // Helper for _bindc()
    void _bind_input(ring_buffer_dict &rb_dict, int ix, bool create, shared_ptr<ring_buffer> &buf, int &mode)
    {
	if (create) {
	    buf = this->create_buffer(rb_dict, bufname_from_index(ix), { params.nfreq }, 1);
	    mode = ring_buffer::ACCESS_APPEND;
	}
	else {
	    buf = this->get_buffer(rb_dict, bufname_from_index(ix));
	    mode = ring_buffer::ACCESS_READ;
	}

	rf_assert(buf->csize == params.nfreq);
	rf_assert(buf->nds == 1);
    }


    // Helper for _bindc()
    void _bind_output(ring_buffer_dict &rb_dict, int ix, bool create, shared_ptr<ring_buffer> &buf, int &mode, ring_buffer_subarray* &arr)
    {
	if (ix == params.ix_in0) {
	    rf_assert(!create);
	    mode_in0 |= ring_buffer::ACCESS_WRITE;
	    arr = &arr_in0;
	}
	else if (ix == params.ix_in1) {
	    rf_assert(!create);
	    mode_in1 |= ring_buffer::ACCESS_WRITE;
	    arr = &arr_in1;
	}
	else if (create) {
	    buf = this->create_buffer(rb_dict, bufname_from_index(ix), { params.nfreq }, 1);
	    mode = ring_buffer::ACCESS_APPEND;
	}
	else {
	    buf = this->get_buffer(rb_dict, bufname_from_index(ix));
	    mode = ring_buffer::ACCESS_WRITE;
	}
    }
    

    virtual void _bindc(ring_buffer_dict &rb_dict, Json::Value &json_attrs) override
    {
	_bind_input(rb_dict, params.ix_in0, params.create_in0, buf_in0, mode_in0);
	_bind_input(rb_dict, params.ix_in1, params.create_in1, buf_in1, mode_in1);
	_bind_output(rb_dict, params.ix_out0, params.create_out0, buf_out0, mode_out0, arr_out0);
	_bind_output(rb_dict, params.ix_out1, params.create_out1, buf_out1, mode_out1, arr_out1);
    }


    virtual bool _process_chunk(ssize_t pos) override
    {
	arr_in0.get(buf_in0, pos, pos + nt_chunk, mode_in0);
	arr_in1.get(buf_in1, pos, pos + nt_chunk, mode_in1);

	if (buf_out0)
	    _arr_out0.get(buf_out0, pos, pos + nt_chunk, mode_out0);
	if (buf_out1)
	    _arr_out1.get(buf_out1, pos, pos + nt_chunk, mode_out1);

	if (params.create_in0)
	    _randomize(arr_in0, rng_in0);

	if (params.create_in1)
	    _randomize(arr_in1, rng_in1);

	for (ssize_t ic = 0; ic < params.nfreq; ic++) {
	    _rotate(arr_out0->data + ic * arr_out0->stride,
		    arr_out1->data + ic * arr_out1->stride,
		    arr_in0.data + ic * arr_in0.stride,
		    arr_in1.data + ic * arr_in1.stride,
		    nt_chunk, cos_theta, sin_theta);
	}

	arr_in0.reset();
	arr_in1.reset();
	_arr_out0.reset();
	_arr_out1.reset();

	return (pos < nt_end);
    }


    virtual Json::Value jsonize() const override
    {
	Json::Value ret;
	ret["class_name"] = "rot2";
	ret["ix_in0"] = params.ix_in0;
	ret["ix_in1"] = params.ix_in1;
	ret["ix_out0"] = params.ix_out0;
	ret["ix_out1"] = params.ix_out1;
	ret["seed_in0"] = params.seed_in0;
	ret["seed_in1"] = params.seed_in1;
	ret["theta"] = params.theta;

	return ret;
    }
};


// Called in reference pipeline
static void _randomize(vector<float> &v, std::mt19937 &rng, ssize_t nc, ssize_t nt)
{
    rf_assert(v.size() == 0);
    v.resize(nc * nt, 0.0);

    // Note loop ordering!
    for (ssize_t it = 0; it < nt; it++)
	for (ssize_t ic = 0; ic < nc; ic++)
	    v[it*nc + ic] = uniform_rand(rng, -1.0, 1.0);
}


static void _rotate(vector<float> &dst1, vector<float> &dst2, vector<float> &src1, vector<float> &src2, float cos_theta, float sin_theta)
{
    rf_assert(src1.size() > 0);
    rf_assert(src2.size() == src1.size());
    rf_assert(dst1.size() == src1.size());
    rf_assert(dst2.size() == src1.size());

    _rotate(&dst1[0], &dst2[0], &src1[0], &src2[0], src1.size(), cos_theta, sin_theta);
}


struct reference_rot2 : public reference_pipeline_object {
    const rot2_params params;
    float cos_theta = 0.0;
    float sin_theta = 0.0;

    std::mt19937 rng_in0;
    std::mt19937 rng_in1;


    reference_rot2(const rot2_params &params_) :
	params(params_),
	rng_in0(params_.seed_in0),
	rng_in1(params_.seed_in1)
    {
	this->csize = params.nfreq;
	this->nds = 1;
	this->cos_theta = cos(params.theta);
	this->sin_theta = sin(params.theta);
    }


    virtual void apply_reference(vector<vector<float>> &v, ssize_t nt_end) override
    {
	int ix_max = -1;
	ix_max = max(ix_max, params.ix_in0);
	ix_max = max(ix_max, params.ix_in1);
	ix_max = max(ix_max, params.ix_out0);
	ix_max = max(ix_max, params.ix_out1);

	if (ix_max >= (int)v.size())
	    v.resize(ix_max+1);
	
	if (params.create_in0)
	    _randomize(v[params.ix_in0], rng_in0, params.nfreq, nt_end);

	if (params.create_in1)
	    _randomize(v[params.ix_in1], rng_in1, params.nfreq, nt_end);

	if (params.create_out0) {
	    rf_assert(v[params.ix_out0].size() == 0);
	    v[params.ix_out0].resize(params.nfreq * nt_end, 0.0);
	}

	if (params.create_out1) {
	    rf_assert(v[params.ix_out1].size() == 0);
	    v[params.ix_out1].resize(params.nfreq * nt_end, 0.0);
	}
	
	_rotate(v[params.ix_out0], v[params.ix_out1], v[params.ix_in0], v[params.ix_in1], cos_theta, sin_theta);
    }


    virtual shared_ptr<pipeline_object> make_real_pipeline_object(ssize_t nt_end) override
    {
	return make_shared<rot2> (this->params, nt_end);
    }
};


// -------------------------------------------------------------------------------------------------


static shared_ptr<reference_pipeline_object>
make_random_pipeline(std::mt19937 &rng, int nfreq, vector<bool> &bflag, int ntransforms)
{
    rf_assert(ntransforms > 0);
    
    if ((ntransforms > 1) || (uniform_rand(rng) < 0.1)) {
	// Make pipeline.
	// Length of pipeline chain.
	int np = randint(rng, 1, min(ntransforms+1,9));

	// Number of transforms in each step in chain.
	vector<int> nt(np, 1);
	for (int i = np; i < ntransforms; i++)
	    nt[randint(rng,0,np)]++;

	auto ret = make_shared<reference_pipeline> ();
	for (int i = 0; i < np; i++)
	    ret->add(make_random_pipeline(rng, nfreq, bflag, nt[i]));

	return ret;
    }

    // Make rot2
    rot2_params p = rot2_params::make_random(rng, nfreq, bflag);
    return make_shared<reference_rot2> (p);
}


static shared_ptr<reference_pipeline_object>
make_random_pipeline(std::mt19937 &rng, int nfreq, int nbuffers, int ntransforms)
{
    vector<bool> bflag(nbuffers, false);
    return make_random_pipeline(rng, nfreq, bflag, ntransforms);
}



int main(int argc, char **argv)
{
    const int niter = 1000;

    std::random_device rd;
    std::mt19937 rng(rd());
    
    for (int iter = 0; iter < niter; iter++) {
	if (iter % 50 == 0)
	    cout << "test-core-pipeline-logic: iteration " << iter << "/" << niter << endl;

	ssize_t nfreq = randint(rng, 1, 20);
	ssize_t nbuffers = randint(rng, 2, 10);
	ssize_t ntransforms = randint(rng, 1, 30);
	ssize_t nt_end = randint(rng, 1000, 10000);

	auto p = make_random_pipeline(rng, nfreq, nbuffers, ntransforms);
	// Json::Value j = p->jsonize();
	// Json::StyledWriter w;
	// cout << w.write(j);
	
	p->run_test(nt_end);
    }

    cout << "test-core-pipeline-logic: pass" << endl;
    return 0;
}
