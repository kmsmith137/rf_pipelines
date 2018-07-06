#include "rf_pipelines_internals.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif

struct mask_counter_measurements {
    ssize_t pos;
    int nsamples;
    int nsamples_masked;
    int nt;
    int nt_masked;
    int nf;
    int nf_masked;
};

struct mask_counter_transform : public wi_transform {

    unique_ptr<bool[]> any_unmasked_t;
    // Ring buffer of measurements
    int max_measurements;
    vector<mask_counter_measurements> measurements;
    int imeasurement;
    int nmeasurements;

    mask_counter_transform(int nt_chunk_, max_measurements_) :
        wi_transform("mask_counter"),
        max_measurements(max_measurements_),
        measurements(max_measurements),
        imeasurement(0)
    {	
        stringstream ss;
        ss << "mask_counter(nt_chunk=" << nt_chunk_ << ")";
        this->name = ss.str();
        this->nt_chunk = nt_chunk_;
        this->nds = 0;  // allows us to run in a wi_sub_pipeline

        if (nt_chunk == 0)
            throw runtime_error("rf_pipelines::mask_counter: nt_chunk must be specified");
	
        // Can't construct the kernel yet, since 'nfreq' is not known until set_stream()
    }

    virtual ~mask_counter_transform() { }

    // Called after (nfreq, nds) are initialized.
    virtual void _bind_transform(Json::Value &json_attrs) override
    {
        this->any_unmasked_t = unique_ptr<bool[]>(new bool[nt_chunk/nds]);
    }

    virtual void _process_chunk(float *intensity, ssize_t istride, float *weights, ssize_t wstride, ssize_t pos) override
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

        int nmasked = 0;
        int nfmasked = 0;

        for (int i_f=0; i_f<nfreq; i_f++) {
            bool allmasked = true;
            for (int i_t=0; i_t<nt; i_t++) {
                if (weights[i_f*wstride + i_t] == 0) {
                    nmasked++;
                } else {
                    allmasked = false;
                    any_unmasked_t[i_t] = true;
                }
            }
            if (allmasked)
                nfmasked++;
        }

        int ntmasked = 0;
        for (int i=0; i<nt; i++)
            if (!any_unmasked_t[i])
                ntmasked++;

        cout << "pos " << pos << ": N samples masked: " << nmasked << "/" << (nfreq*nt) << "; n times " << ntmasked << "/" << nt << "; n freqs " << nfmasked << "/" << nfreq << endl;
    }

    virtual Json::Value jsonize() const override
    {
        Json::Value ret;

        ret["class_name"] = "mask_counter";
        ret["nt_chunk"] = int(this->get_prebind_nt_chunk());
	
        return ret;
    }

    static shared_ptr<mask_counter_transform> from_json(const Json::Value &j)
    {
        ssize_t nt_chunk = ssize_t_from_json(j, "nt_chunk");
        int max_measurements = int_from_json(j, "max_measurements");
        return make_shared<mask_counter_transform> (nt_chunk, max_measurements);
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
shared_ptr<wi_transform> make_mask_counter(int nt_chunk)
{
    return make_shared<mask_counter_transform> (nt_chunk);
}


}  // namespace rf_pipelines
