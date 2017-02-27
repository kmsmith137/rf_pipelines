#include <random>
#include "rf_pipelines_internals.hpp"
#include "chime_file_stream_base.hpp"

using namespace std;
using namespace rf_pipelines;


// -------------------------------------------------------------------------------------------------
//
// test_stream


class test_stream : public chime_file_stream_base {
public:
    // Parameters of the test_stream not specified in the constructor (e.g. nfreq) will be randomly generated.
    test_stream(const vector<string> &filename_list, int nt_chunk, int noise_source_align, int nsamples_per_noise_source_switch, int it_initial);

    // Factory function which additionally randomizes (nfiles, nt_chunk, noise_source_align).
    static shared_ptr<test_stream> make_random();

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
    virtual void read_data(float *dst_int, float *dst_wt, ssize_t it_file, ssize_t n, ssize_t dst_stride) const override;
    virtual void close_file() override;
};


// static member function
shared_ptr<test_stream> test_stream::make_random()
{
    int nfiles = randint(1, 11);
    int nt_chunk = randint(1, 11);
    int noise_source_align = randint(0,2) ? (1 << randint(1,5)) : 0;
    int nsamples_per_noise_source_switch = 1 << randint(4,7);
    int it_initial = randint(0, 100);

    vector<string> filename_list(nfiles);
    for (int i = 0; i < nfiles; i++)
	filename_list[i] = to_string(i);

    return make_shared<test_stream> (filename_list, nt_chunk, noise_source_align, nsamples_per_noise_source_switch, it_initial);
}


test_stream::test_stream(const vector<string> &filename_list_, int nt_chunk_, int noise_source_align_, int nsamples_per_noise_source_switch_, int it_initial_) :
    chime_file_stream_base(filename_list_, nt_chunk_, noise_source_align_)
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
	this->_file_nt[0] = randint(f, f+10);
	
	for (int i = 1; i < _nfiles; i++) {
	    this->_file_it[i] = _file_it[i-1] + _file_nt[i-1] + randint(0,11);
	    this->_file_nt[i] = randint(1, 11);
	    this->_file_freq_inc[i] = randint(0, 2);
	}
	
	this->_it1 = _file_it[_nfiles-1] + _file_nt[_nfiles-1];
	f++;
    } while (_it1 <= _it0 + noise_source_align_);
    
    this->_nfreq = randint(1, 11);
    this->_freq_lo_MHz = uniform_rand(200.0, 600.0);
    this->_freq_hi_MHz = uniform_rand(600.0, 1000.0);
    this->_dt_sample = (1 << 23) * constants::chime_seconds_per_fpga_count / double(nsamples_per_noise_source_switch_);
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
    this->nt = _file_nt[_curr_ifile];
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
void test_stream::read_data(float *dst_int, float *dst_wt, ssize_t it_file, ssize_t n, ssize_t dst_stride) const
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

	    dst_int[i*dst_stride + j] = _intensity(ifreq, it);
	    dst_wt[i*dst_stride + j] = _weight(ifreq, it);
	}
    }
}


// -------------------------------------------------------------------------------------------------
//
// test_transform


struct test_transform : public wi_transform {
    shared_ptr<test_stream> s;
    int it_curr = 0;

    test_transform(const shared_ptr<test_stream> &s_);

    virtual void set_stream(const wi_stream &stream) override;
    virtual void start_substream(int isubstream, double t0) override;
    virtual void process_chunk(double t0, double t1, float *intensity, float *weights, ssize_t stride, float *pp_intensity, float *pp_weights, ssize_t pp_stride) override;
    virtual void end_substream() override;
};


test_transform::test_transform(const shared_ptr<test_stream> &s_)
{
    this->s = s_;
    this->nt_chunk = randint(1, 11);
    this->name = "test_transform";

    int n = s->_noise_source_align;
    this->it_curr = (n > 0) ? round_up(s->_it0,n) : s->_it0;
}

// virtual override
void test_transform::set_stream(const wi_stream &stream)
{
    rf_assert(stream.nfreq == s->_nfreq);
    rf_assert(stream.freq_lo_MHz == s->_freq_lo_MHz);
    rf_assert(stream.freq_hi_MHz == s->_freq_hi_MHz);
    rf_assert(stream.dt_sample == s->_dt_sample);

    this->nfreq = stream.nfreq;
}

// virtual override
void test_transform::start_substream(int isubstream, double t0)
{
    double expected_t0 = s->_dt_sample * it_curr;
    rf_assert(fabs(t0 - expected_t0) < 1.0e-3);
}

// virtual override
void test_transform::process_chunk(double t0, double t1, float *intensity, float *weights, ssize_t stride, float *pp_intensity, float *pp_weights, ssize_t pp_stride)
{
    double expected_t0 = s->_dt_sample * it_curr;
    double expected_t1 = s->_dt_sample * (it_curr + nt_chunk);

    rf_assert(fabs(t0 - expected_t0) < 1.0e-3);
    rf_assert(fabs(t1 - expected_t1) < 1.0e-3);

    for (int j = 0; j < nt_chunk; j++) {
	if (s->is_masked(it_curr+j)) {
	    for (int i = 0; i < nfreq; i++)
		rf_assert(weights[i*stride+j] == 0.0);
	}
	else {
	    for (int i = 0; i < nfreq; i++) {
		rf_assert(intensity[i*stride+j] == test_stream::_intensity(i,it_curr+j));
		rf_assert(weights[i*stride+j] == test_stream::_weight(i,it_curr+j));
	    }
	}
    }

    this->it_curr += nt_chunk;
}

// virtual override
void test_transform::end_substream()
{
    rf_assert(this->it_curr >= s->_it1);
}


// -------------------------------------------------------------------------------------------------


int main(int argc, char **argv)
{
    cerr << "test-file-stream-base..";

    for (int iter = 0; iter < 1000; iter++) {
	if (iter % 10 == 0)
	    cerr << ".";

	shared_ptr<test_stream> sp = test_stream::make_random();
	shared_ptr<test_transform> tp = make_shared<test_transform> (sp);

	// run pipeline with no outputs
	sp->run({tp}, "", nullptr, 0);
    }

    cerr << "pass\n";
    return 0;
}
