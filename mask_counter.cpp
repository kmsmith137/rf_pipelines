#include "rf_pipelines_internals.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif


struct mask_counter_transform : public wi_transform {

    mask_counter_transform(int nt_chunk_) :
        wi_transform("mask_counter")
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
        //this->kernel = make_unique<rf_kernels::mask_counter> (nfreq, xdiv(nt_chunk,nds), axis, sigma, Df, Dt, two_pass);
    }

    virtual void _process_chunk(float *intensity, ssize_t istride, float *weights, ssize_t wstride, ssize_t pos) override
    {
        cout << "mask_counters process_chunk" << endl;
        //this->kernel->clip(intensity, istride, weights, wstride);
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
        bool*  any_unmasked_t = (bool*)malloc(nt);
        memset(any_unmasked_t, 0, nt);

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

        free(any_unmasked_t);
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
        return make_shared<mask_counter_transform> (nt_chunk);
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
