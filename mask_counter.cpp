#include "rf_pipelines_internals.hpp"
#include "rf_pipelines_inventory.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif

mask_counter_transform::mask_counter_transform(int nt_chunk_, string where_) :
    wi_transform("mask_counter"),
    where(where_)
    {	
        stringstream ss;
        ss << "mask_counter(nt_chunk=" << nt_chunk_ << ", where=" << where << ")";
        this->name = ss.str();
        this->nt_chunk = nt_chunk_;
        this->nds = 0;  // allows us to run in a wi_sub_pipeline

        if (nt_chunk == 0)
            throw runtime_error("rf_pipelines::mask_counter: nt_chunk must be specified");
	
        // Can't construct the kernel yet, since 'nfreq' is not known until set_stream()
    }

// Called after (nfreq, nds) are initialized.
void mask_counter_transform::_bind_transform(Json::Value &json_attrs)
{
    this->any_unmasked_t = unique_ptr<bool[]>(new bool[nt_chunk/nds]);
}

void mask_counter_transform::_process_chunk(float *intensity, ssize_t istride, float *weights, ssize_t wstride, ssize_t pos)
    {
        // _process_chunk(intensity, istride, weights, wstride, pos)
        //
        //    This is the "core" method which is responsible for processing the 'intensity' and 'weights'
        //    arrays over timestamp range [pos,pos+nt_chunk).  The array strides are istride/wstride, 
        //    i.e. the (i,j)-th element of the intensity  array is intensity[i*istride+j], and similarly 
        //    for the weights.
        //
        //    Transforms which support downsampling (i.e. nds > 1) should note that 'pos'
        //    and 'nt_chunk' do not have the downsampling factor applied, but array dimensions
        //    do.  That is, the 'intensity' and 'weights' arrays have shape (nfreq, nt_chunk/nds),
        //    not shape (nfreq, nt_chunk), and 'pos' increases by nt_chunk (not nt_chunk/nds)
        //    in each call to _process_chunk();
        int nt = nt_chunk/nds;
        memset(any_unmasked_t.get(), 0, nt);

        mask_counter_measurements meas;
        memset(&meas, 0, sizeof(mask_counter_measurements));
        meas.pos = pos;
        meas.nsamples = nfreq*nt;
        meas.nt = nt;
        meas.nf = nfreq;

        meas.nsamples_masked = 0;
        meas.nf_masked = 0;

        for (int i_f=0; i_f<nfreq; i_f++) {
            bool allmasked = true;
            for (int i_t=0; i_t<nt; i_t++) {
                if (weights[i_f*wstride + i_t] == 0) {
                    meas.nsamples_masked++;
                } else {
                    allmasked = false;
                    any_unmasked_t[i_t] = true;
                }
            }
            if (allmasked)
                meas.nf_masked++;
        }

        meas.nt_masked = 0;
        for (int i=0; i<nt; i++)
            if (!any_unmasked_t[i])
                meas.nt_masked++;

        cout << "pos " << pos << ": N samples masked: " << meas.nsamples_masked << "/" << (meas.nsamples) << "; n times " << meas.nt_masked << "/" << meas.nt << "; n freqs " << meas.nf_masked << "/" << meas.nf << endl;

        cout << "Calling callback on " << callbacks.size() << " objects" << endl;
        for (const auto &cb : callbacks)
            cb->mask_count(meas);
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
