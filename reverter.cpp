#include <iostream>
//#include <deque>
#include "rf_pipelines_internals.hpp"
#include "rf_pipelines.hpp"
using namespace std;

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif


pair<shared_ptr<wi_transform>, shared_ptr<wi_transform> >
make_reverter(ssize_t nfreq, ssize_t nt_chunk);

class Reverter;

class Saver : public wi_transform {
    friend class Reverter;
public:
    Saver(ssize_t nt_chunk) :
        wi_transform(),
        intensity(NULL),
        weight(NULL),
        filled(false)
    {
        // base class members
        this->nt_chunk = nt_chunk;
        this->nfreq = 0;
        this->name = "Saver";
    }

    virtual ~Saver() {
        free(intensity);
        free(weight);
    }

    virtual void start_substream(int isubstream, double t0) {}
    virtual void end_substream() {}

    virtual void set_stream(const wi_stream &stream) {
        this->nfreq = stream.nfreq;
        free(intensity);
        free(weight);
        intensity = aligned_alloc<float>(nfreq * nt_chunk);
        weight    = aligned_alloc<float>(nfreq * nt_chunk);
    }

    virtual void process_chunk(double t0, double t1, float* ii, float* ww,
                               ssize_t stride, float* pp_ii, float* pp_ww, ssize_t pp_stride) {
        cout << "Saver: process_chunk, t " << t0 << " to " << t1 << endl;
        if (filled)
            throw runtime_error("Saver: process_chunk called but I already have a chunk buffered!");
        for (ssize_t i=0; i<nfreq; i++)
            memcpy(intensity + i*nt_chunk, ii + i*stride, sizeof(float)*nt_chunk);
        for (ssize_t i=0; i<nfreq; i++)
            memcpy(weight + i*nt_chunk,    ww + i*stride, sizeof(float)*nt_chunk);
        filled = true;
    }

protected:
    //ssize_t nfreq;
    //ssize_t nt_chunk;

    float* intensity;
    float* weight;
    
    bool filled;

    void revert_chunk(float* ii, float* ww, ssize_t stride) {
        if (!filled)
            throw runtime_error("Saver: revert_chunk called but I do not have a chunk buffered!");
        for (ssize_t i=0; i<nfreq; i++)
            memcpy(ii + i*stride, intensity + i*nt_chunk, sizeof(float)*nt_chunk);
        for (ssize_t i=0; i<nfreq; i++)
            memcpy(ww + i*stride, weight + i*nt_chunk,    sizeof(float)*nt_chunk);
        filled = false;
    }

};

class Reverter : public wi_transform {
public:
    Reverter(shared_ptr<Saver> s) :
        wi_transform(),
        saver(s)
    {
        // base class members
        this->nt_chunk = s->nt_chunk;
        this->nfreq = 0;
        this->name = "Reverter";
    }

    virtual ~Reverter() {}

    virtual void set_stream(const wi_stream &stream) {
        this->nfreq = stream.nfreq;
    }
    virtual void start_substream(int isubstream, double t0) {}
    virtual void end_substream() {}
    
    virtual void process_chunk(double t0, double t1, float* ii, float* ww,
                               ssize_t stride, float* pp_ii, float* pp_ww, ssize_t pp_stride) {
        cout << "Reverter: process_chunk, t " << t0 << " to " << t1 << endl;
        saver->revert_chunk(ii, ww, stride);
    }

protected:
    shared_ptr<Saver> saver;
};


pair<shared_ptr<wi_transform>, shared_ptr<wi_transform> >
make_reverter(ssize_t nt_chunk) {
    shared_ptr<Saver> s = make_shared<Saver>(nt_chunk);
    shared_ptr<Reverter> r = make_shared<Reverter>(s);
    return make_pair(s, r);
}

}
