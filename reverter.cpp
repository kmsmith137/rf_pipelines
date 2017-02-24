#include <iostream>
//#include <deque>
#include "rf_pipelines_internals.hpp"
#include "rf_pipelines.hpp"
#include "reverter.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif

Saver::Saver(ssize_t nt_chunk) :
    wi_transform(),
    intensity(NULL),
    weight(NULL)
{
    // base class members
    this->nt_chunk = nt_chunk;
    this->nfreq = 0;
    this->name = "Saver";
}

Saver::~Saver() {
    free(intensity);
    free(weight);
}

void Saver::start_substream(int isubstream, double t0) {}
void Saver::end_substream() {}

void Saver::set_stream(const wi_stream &stream) {
    this->nfreq = stream.nfreq;
    free(intensity);
    free(weight);
    intensity = aligned_alloc<float>(nfreq * nt_chunk);
    weight    = aligned_alloc<float>(nfreq * nt_chunk);
}

void Saver::process_chunk(double t0, double t1,
                          float* ii, float* ww, ssize_t stride,
                          float* pp_ii, float* pp_ww, ssize_t pp_stride) {
    cout << "Saver: process_chunk, t " << t0 << " to " << t1 << endl;
    /*
     if (filled)
     throw runtime_error("Saver: process_chunk called but I already have a chunk buffered!");
     */
    for (ssize_t i=0; i<nfreq; i++)
        memcpy(intensity + i*nt_chunk, ii + i*stride, sizeof(float)*nt_chunk);
    for (ssize_t i=0; i<nfreq; i++)
        memcpy(weight + i*nt_chunk,    ww + i*stride, sizeof(float)*nt_chunk);
    this->t0 = t0;
    this->t1 = t1;
}

void Saver::revert_chunk(double t0, double t1,
                         float* ii, float* ww, ssize_t stride) {
    /*
     if (!filled)
     throw runtime_error("Saver: revert_chunk called but I do not have a chunk buffered!");
     */
    if (t0 != this->t0 || t1 != this->t1)
        throw runtime_error("Saver: revert_chunk called for time range " + to_string(t0) + " to " + to_string(t1) + " but I have range " + to_string(this->t0) + " to " + to_string(this->t1) + " buffered!");
            
    for (ssize_t i=0; i<nfreq; i++)
        memcpy(ii + i*stride, intensity + i*nt_chunk, sizeof(float)*nt_chunk);
    for (ssize_t i=0; i<nfreq; i++)
        memcpy(ww + i*stride, weight + i*nt_chunk,    sizeof(float)*nt_chunk);
    //filled = false;
}

Reverter::Reverter(shared_ptr<Saver> s) :
    wi_transform(),
    saver(s)
{
    // base class members
    this->nt_chunk = s->nt_chunk;
    this->nfreq = 0;
    this->name = "Reverter";
}

Reverter::~Reverter() {}

void Reverter::set_stream(const wi_stream &stream) {
    this->nfreq = stream.nfreq;
}

void Reverter::start_substream(int isubstream, double t0) {}
void Reverter::end_substream() {}
    
void Reverter::process_chunk(double t0, double t1,
                             float* ii, float* ww, ssize_t stride,
                             float* pp_ii, float* pp_ww, ssize_t pp_stride) {
    cout << "Reverter: process_chunk, t " << t0 << " to " << t1 << endl;
    saver->revert_chunk(t0, t1, ii, ww, stride);
}

pair<shared_ptr<wi_transform>, shared_ptr<wi_transform> >
make_reverter(ssize_t nt_chunk) {
    shared_ptr<Saver> s = make_shared<Saver>(nt_chunk);
    shared_ptr<Reverter> r = make_shared<Reverter>(s);
    return make_pair(s, r);
}

std::shared_ptr<Saver> make_saver(ssize_t nt_chunk) {
    return make_shared<Saver>(nt_chunk);
}

}
