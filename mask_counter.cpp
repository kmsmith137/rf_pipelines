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
    stringstream ss;
    ss << class_name_ << "(nt_chunk=" << nt_chunk_ << ", where=" << where << ")";
    this->name = ss.str();
    this->nt_chunk = nt_chunk_;
    this->nds = 0;  // allows us to run in a wi_sub_pipeline

    if (nt_chunk == 0)
        throw runtime_error("rf_pipelines::mask_counter: nt_chunk must be specified");
	
    // Can't construct the kernel yet, since 'nfreq' is not known until set_stream()
}

void mask_counter_transform::_process_chunk(float *intensity, ssize_t istride, float *weights, ssize_t wstride, ssize_t pos)
{
    if (callbacks.size() == 0) {
        cout << "mask_counter " << where << ": no callbacks; not counting" << endl;
        return;
    }
    int nt = nt_chunk/nds;

    mask_measurements meas;
    meas.pos = pos;
    meas.nsamples = nfreq*nt;
    meas.nsamples_masked = 0;
    meas.nt = nt;
    meas.nt_masked = 0;
    meas.nf = nfreq;
    meas.nf_masked = 0;
    meas.freqs_masked = shared_ptr<uint16_t>((uint16_t*)calloc(nfreq, sizeof(uint16_t)), free);
    meas.times_masked = shared_ptr<uint16_t>((uint16_t*)calloc(nt,    sizeof(uint16_t)), free);

    uint16_t* fm = meas.freqs_masked.get();
    uint16_t* tm = meas.times_masked.get();

    for (int i_f=0; i_f<nfreq; i_f++) {
        for (int i_t=0; i_t<nt; i_t++) {
            if (weights[i_f*wstride + i_t] == 0) {
                meas.nsamples_masked++;
                fm[i_f]++;
                tm[i_t]++;
            }
        }
    }

    for (int i_f=0; i_f<nfreq; i_f++)
        if (fm[i_f] == nt)
            meas.nf_masked++;
    for (int i_t=0; i_t<nt; i_t++)
        if (tm[i_t] == nfreq)
            meas.nt_masked++;
            
    cout << "mask_counter " << where << ", pos " << pos << ": N samples masked: " << meas.nsamples_masked << "/" << (meas.nsamples) << "; n times " << meas.nt_masked << "/" << meas.nt << "; n freqs " << meas.nf_masked << "/" << meas.nf << endl;
    //cout << "mask_counter: calling " << callbacks.size() << " callbacks" << endl;
    for (const auto &cb : callbacks)
        cb->mask_count(meas);
}

void mask_counter_transform::process_measurement()
{
    // add to ringbuf
}

void mask_counter_transform::init_measurements(mask_measurements& meas) {
    // alloc arrays, etc
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

void mask_counter_transform::add_callback(const std::shared_ptr<mask_counter_callback> cb) {
    callbacks.push_back(cb);
}

void mask_counter_transform::remove_callback(const std::shared_ptr<mask_counter_callback> cb) {
    for (auto it=callbacks.begin(); it!=callbacks.end(); it++) {
        if (*it == cb) {
            callbacks.erase(it);
            break;
        }
    }
};


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
