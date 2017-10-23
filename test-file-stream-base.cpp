#include "rf_pipelines_internals.hpp"

using namespace std;
using namespace rf_pipelines;


// -------------------------------------------------------------------------------------------------
//
// test_stream


class test_stream : public chime_file_stream_base {
public:
    // Parameters of the test_stream not specified in the constructor (e.g. nfreq) will be randomly generated.
    test_stream(std::mt19937 &rng, const vector<string> &filename_list, int nt_chunk, int noise_source_align, int nsamples_per_noise_source_switch, int it_initial);

    // Factory function which additionally randomizes (nfiles, nt_chunk, noise_source_align).
    static shared_ptr<test_stream> make_random(std::mt19937 &rng);

    // These fields begin with underscores, to distinguish them from fields of the 
    // base classes (wi_stream, chime_file_stream_base) with the same names.

    int _nfiles = 0;
    int _curr_ifile = -1;
    int _it0 = 0;
    int _it1 = 0;

    vector<int> _file_it;
    vector<int> _file_nt;
    vector<bool> _file_freq_inc;   // 'frequencies_are_increasing' flag

    int _nfreq = 0;
    double _freq_lo_MHz = 0.0;
    double _freq_hi_MHz = 0.0;
    double _dt_sample = 0.0;
    int _noise_source_align = 0;

    // An arbitrary but nonrandom (intensity, weight) for every (ifreq, it).
    static inline float _intensity(int ifreq, int it) { return sin(1.83*ifreq + 0.329*it); }
    static inline float _weight(int ifreq, int it) { return 1.0 + cos(0.71*ifreq - 0.87*it); }

    bool is_masked(int it) const;

protected:
    // Override virtual member functions in base class chime_file_stream_base.
    virtual void load_file(const std::string &filename) override;
    virtual void set_params_from_file() override;
    virtual void check_file_consistency() const override;
    virtual void read_data(float *dst_int, float *dst_wt, ssize_t it_file, ssize_t n, ssize_t dst_istride, ssize_t dst_wstride) const override;
    virtual void close_file() override;
};


// static member function
shared_ptr<test_stream> test_stream::make_random(std::mt19937 &rng)
{
    int nfiles = randint(rng, 1, 11);
    int nt_chunk = randint(rng, 1, 11);
    int noise_source_align = randint(rng,0,2) ? (1 << randint(rng,1,5)) : 0;
    int nsamples_per_noise_source_switch = 1 << randint(rng,4,7);
    int it_initial = randint(rng, 0, 100);

    vector<string> filename_list(nfiles);
    for (int i = 0; i < nfiles; i++)
	filename_list[i] = to_string(i);

    return make_shared<test_stream> (rng, filename_list, nt_chunk, noise_source_align, nsamples_per_noise_source_switch, it_initial);
}


test_stream::test_stream(std::mt19937 &rng, const vector<string> &filename_list_, int nt_chunk_, int noise_source_align_, int nsamples_per_noise_source_switch_, int it_initial_) :
    chime_file_stream_base("test_stream", filename_list_, nt_chunk_, noise_source_align_)
{
    this->_nfiles = filename_list_.size();
    this->_curr_ifile = -1;
    this->_noise_source_align = noise_source_align_;
    this->_it0 = it_initial_;

    this->_file_it.resize(_nfiles);
    this->_file_nt.resize(_nfiles);
    this->_file_freq_inc.resize(_nfiles);
	
    // This loop generates a list of "files" (i.e. time ranges)
    // whose total length is at least 'noise_source_align'.

    int f = 0;
    do {
	this->_file_it[0] = _it0;
	this->_file_nt[0] = f + randint(rng,0,10);  // randint(rng,f,f+10) gives spurious gcc5 warning!
	
	for (int i = 1; i < _nfiles; i++) {
	    this->_file_it[i] = _file_it[i-1] + _file_nt[i-1] + randint(rng,0,11);
	    this->_file_nt[i] = randint(rng, 1, 11);
	    this->_file_freq_inc[i] = randint(rng, 0, 2);
	}
	
	this->_it1 = _file_it[_nfiles-1] + _file_nt[_nfiles-1];
	f++;
    } while (_it1 <= _it0 + noise_source_align_);
    
    this->_nfreq = randint(rng, 1, 11);
    this->_freq_lo_MHz = uniform_rand(rng, 200.0, 600.0);
    this->_freq_hi_MHz = uniform_rand(rng, 600.0, 1000.0);
    this->_dt_sample = (1 << 23) * chime_seconds_per_fpga_count / double(nsamples_per_noise_source_switch_);
}


bool test_stream::is_masked(int it) const
{
    for (int i = 0; i < _nfiles; i++)
	if ((it >= _file_it[i]) && (it < _file_it[i] + _file_nt[i]))
	    return false;

    return true;
}

// virtual override
void test_stream::load_file(const std::string &filename)
{
    int ifile = lexical_cast<int> (filename);

    rf_assert(_curr_ifile < 0);
    rf_assert(ifile >= 0 && ifile < _nfiles);
    this->_curr_ifile = ifile;
}

// virtual override
void test_stream::close_file()
{
    rf_assert(_curr_ifile >= 0);
    this->_curr_ifile = -1;
}

// virtual override
void test_stream::set_params_from_file()
{
    rf_assert(_curr_ifile >= 0);
    rf_assert(_curr_ifile < _nfiles);

    this->nfreq = _nfreq;
    this->freq_lo_MHz = _freq_lo_MHz;
    this->freq_hi_MHz = _freq_hi_MHz;
    this->dt_sample = _dt_sample;
    this->time_lo = _dt_sample * (_file_it[_curr_ifile]);
    this->time_hi = _dt_sample * (_file_it[_curr_ifile] + _file_nt[_curr_ifile]);
    this->nt_file = _file_nt[_curr_ifile];
    this->frequencies_are_increasing = _file_freq_inc[_curr_ifile];
}

// virtual override
void test_stream::check_file_consistency() const
{
    rf_assert(nfreq == _nfreq);
    rf_assert(fabs(freq_lo_MHz - _freq_lo_MHz) < 1.0e-3);
    rf_assert(fabs(freq_hi_MHz - _freq_hi_MHz) < 1.0e-3);
    rf_assert(fabs(dt_sample - _dt_sample) < 1.0e-3);
}

// virtual override
void test_stream::read_data(float *dst_int, float *dst_wt, ssize_t it_file, ssize_t n, ssize_t dst_istride, ssize_t dst_wstride) const
{
    rf_assert(dst_int != nullptr);
    rf_assert(dst_wt != nullptr);
    rf_assert(_curr_ifile >= 0);
    rf_assert(_curr_ifile < _nfiles);

    rf_assert(n > 0);
    rf_assert(it_file >= 0);
    rf_assert(it_file + n <= _file_nt[_curr_ifile]);

    for (int i = 0; i < _nfreq; i++) {
	for (int j = 0; j < n; j++) {
	    int ifreq = _file_freq_inc[_curr_ifile] ? (_nfreq-1-i) : i;
	    int it = _file_it[_curr_ifile] + it_file + j;

	    dst_int[i*dst_istride + j] = _intensity(ifreq, it);
	    dst_wt[i*dst_wstride + j] = _weight(ifreq, it);
	}
    }
}


// -------------------------------------------------------------------------------------------------
//
// test_transform


struct test_transform : public wi_transform {
    shared_ptr<test_stream> s;
    int it_curr = 0;

    test_transform(std::mt19937 &rng, const shared_ptr<test_stream> &s_);

    virtual void _process_chunk(float *intensity, ssize_t istride, float *weights, ssize_t wstride, ssize_t pos) override;
};


test_transform::test_transform(std::mt19937 &rng, const shared_ptr<test_stream> &s_) :
    wi_transform("test_transform")
{
    this->s = s_;
    this->nt_chunk = randint(rng, 1, 11);

    int n = s->_noise_source_align;
    this->it_curr = (n > 0) ? round_up(s->_it0,n) : s->_it0;
}


// Virtual override
void test_transform::_process_chunk(float *intensity, ssize_t istride, float *weights, ssize_t wstride, ssize_t pos)
{
    rf_assert(this->nfreq == s->_nfreq);
    
    for (int j = 0; j < nt_chunk; j++) {
	if (s->is_masked(it_curr+j)) {
	    for (int i = 0; i < nfreq; i++)
		rf_assert(weights[i*wstride+j] == 0.0);
	}
	else {
	    for (int i = 0; i < nfreq; i++) {
		rf_assert(intensity[i*istride+j] == test_stream::_intensity(i,it_curr+j));
		rf_assert(weights[i*wstride+j] == test_stream::_weight(i,it_curr+j));
	    }
	}
    }

    this->it_curr += nt_chunk;
}


// -------------------------------------------------------------------------------------------------


int main(int argc, char **argv)
{
    std::random_device rd;
    std::mt19937 rng(rd());

    cerr << "test-file-stream-base..";

    for (int iter = 0; iter < 1000; iter++) {
	if (iter % 10 == 0)
	    cerr << ".";

	shared_ptr<test_stream> sp = test_stream::make_random(rng);
	shared_ptr<test_transform> tp = make_shared<test_transform> (rng, sp);

	// run two-element pipeline with no outputs
	vector<shared_ptr<pipeline_object>> v = { sp, tp };
	auto p = make_shared<pipeline> (v, "test_pipeline");

	run_params params;
	params.outdir = "";
	params.verbosity = 0;
	params.debug = true;
	p->run(params);

	rf_assert(tp->it_curr >= sp->_it1);
    }

    cerr << "pass\n";
    return 0;
}
