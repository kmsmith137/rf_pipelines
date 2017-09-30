#include "rf_pipelines_internals.hpp"

#include "simpulse.hpp"

using namespace std;
using namespace simpulse;

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif

class pulse_adder : public wi_transform {
public:
    pulse_adder(ssize_t nt, std::vector<shared_ptr<simpulse::single_pulse> > &pulses, double weight=1.0);
    virtual ~pulse_adder() {}
    virtual void set_stream(const wi_stream &stream);
    virtual void start_substream(int isubstream, double t0) { }
    virtual void process_chunk(double t0, double t1, float *intensity, float *weights, ssize_t stride, float *pp_intensity, float *pp_weights, ssize_t pp_stride);
    virtual void end_substream() { }
protected:
    std::vector<std::shared_ptr<simpulse::single_pulse> > pulses;
    double pulse_weight;
};


pulse_adder::pulse_adder(ssize_t nt, std::vector<shared_ptr<simpulse::single_pulse> > &thepulses, double weight) :
    pulses(thepulses),
    pulse_weight(weight)
{
    this->name = "pulse_adder";
    this->nt_chunk = nt;
}

void pulse_adder::set_stream(const wi_stream &stream) {
    this->nfreq = stream.nfreq;
}


void pulse_adder::process_chunk(double t0, double t1,
                                float *intensity, float *weights, ssize_t stride,
                                float *pp_intensity, float *pp_weights, ssize_t pp_stride)
{
    // hi to lo frequencies -- negative stride, positive offset
    int pulse_stride = -stride;
    float* pulse_array = intensity + stride * (this->nfreq - 1);
        
    for (shared_ptr<single_pulse> pulse : pulses) {
        pulse->add_to_timestream(pulse_array, t0, t1, this->nt_chunk, pulse_stride, pulse_weight);
        double pt0, pt1;
        pulse->get_endpoints(pt0, pt1);
        cout << "Adding pulse with time endpoints " << pt0 << " to " << pt1 << " to chunk with time range " << t0 << " to " << t1 << endl;
    }
}

// Externally visible factory function declared in rf_transforms.hpp
shared_ptr<wi_transform> make_pulse_adder(ssize_t nt, std::vector<shared_ptr<simpulse::single_pulse> > &pulses, double weight) {
    return make_shared<pulse_adder>(nt, pulses, weight);
}

}  // namespace rf_pipelines
