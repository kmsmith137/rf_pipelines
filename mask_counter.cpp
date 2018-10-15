#include <rf_kernels/mask_counter.hpp>

#include "rf_pipelines_internals.hpp"
#include "rf_pipelines_inventory.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif

mask_counter_transform::mask_counter_transform(int nt_chunk_, string where_,
                                               string class_name_) :
    wi_transform(class_name_),
    where(where_)
{	
    ringbuf = make_shared<mask_measurements_ringbuf>();

    stringstream ss;
    ss << class_name_ << "(nt_chunk=" << nt_chunk_ << ", where=" << where << ")";
    this->name = ss.str();
    this->nt_chunk = nt_chunk_;
    this->nds = 0;  // allows us to run in a wi_sub_pipeline

    if (nt_chunk == 0)
        throw runtime_error("rf_pipelines::mask_counter: nt_chunk must be specified");
}

void mask_counter_transform::_process_chunk(float *intensity, ssize_t istride, float *weights, ssize_t wstride, ssize_t pos)
{
    int nt = nt_chunk/nds;
    mask_measurements meas;
    init_measurements(meas);
    meas.pos = pos;

    rf_kernels::mask_counter_data d;
    d.nfreq = nfreq;
    d.nt_chunk = nt;      // not 'nt_chunk'
    d.in = weights;       // not 'intensity'
    d.istride = wstride;  // not 'istride'
    d.out_fcounts = meas.freqs_unmasked.get();
    
    // Run mask-counting kernel.
    meas.nsamples_unmasked = d.mask_count();

    // cout << "mask_counter " << where << ", pos " << pos 
    // << ": N samples masked: " << (meas.nsamples - meas.nsamples_unmasked)
    // << "/" << meas.nsamples << endl;

    process_measurement(meas);
}

void mask_counter_transform::process_measurement(mask_measurements& meas)
{
    ringbuf->add(meas);
}

void mask_counter_transform::init_measurements(mask_measurements& meas) {
    int nt = nt_chunk/nds;
    meas.nf = nfreq;
    meas.nt = nt;
    meas.nsamples = nfreq * nt;
    meas.freqs_unmasked = make_sptr<int> (nfreq);
}

std::shared_ptr<mask_measurements_ringbuf>
mask_counter_transform::get_ringbuf() {
    return ringbuf;
}

Json::Value mask_counter_transform::jsonize() const
    {
        Json::Value ret;

        ret["class_name"] = "mask_counter";
        ret["nt_chunk"] = int(this->get_prebind_nt_chunk());
        ret["where"] = where;
        return ret;
    }


shared_ptr<mask_counter_transform>
mask_counter_transform::from_json(const Json::Value &j)
    {
        ssize_t nt_chunk = ssize_t_from_json(j, "nt_chunk");
        string where = string_from_json(j, "where");
        return make_shared<mask_counter_transform> (nt_chunk, where);
    }

namespace {
    struct _init {
        _init() {
            pipeline_object::register_json_deserializer("mask_counter", mask_counter_transform::from_json);
        }
    } init;
}

// Externally callable
shared_ptr<wi_transform> make_mask_counter(int nt_chunk, string where)
{
    return make_shared<mask_counter_transform> (nt_chunk, where);
}


}  // namespace rf_pipelines
