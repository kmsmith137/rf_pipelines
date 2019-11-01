#include <assert.h>
#include "rf_kernels/spline_detrender.hpp"
#include "rf_pipelines_internals.hpp"
#include "ch_frb_io.hpp"
//#include "ch_frb_io/chlog.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif


class chime_spline_detrender : public spline_detrender, public chime_wi_transform {
public:
    chime_spline_detrender(int nt_chunk_, rf_kernels::axis_type axis_, int nbins_, double epsilon_);
    virtual void _handle_spline(ssize_t pos, float *coeffs, ssize_t cstride, float *intensity, ssize_t istride, float *weights, ssize_t wstride) override;
    virtual void _bind_transform(Json::Value &json_attrs) override;
    virtual Json::Value jsonize() const override;
    static std::shared_ptr<chime_spline_detrender> from_json(const Json::Value &j);
    
};


chime_spline_detrender::chime_spline_detrender(int nt_chunk_, rf_kernels::axis_type axis_, int nbins_, double epsilon_) :
    wi_transform("chime_spline_detrender"),
    spline_detrender(nt_chunk_, axis_, nbins_, epsilon_)
{}

void chime_spline_detrender::_bind_transform(Json::Value &json_attrs) {
    // borrrring
    spline_detrender::_bind_transform(json_attrs);
    chime_wi_transform::_bind_transform(json_attrs);
}

void chime_spline_detrender::_handle_spline(ssize_t pos, float *coeffs, ssize_t cstride, float *intensity, ssize_t istride, float *weights, ssize_t wstride) {
    if (this->chime_stream) {
        auto chunk = this->assembled_chunk_for_pos(pos);
        if (chunk) {
            assert(!chunk->has_detrend_f);
            assert(chunk->detrend_params_f);
            int nco = (this->nbins+1)*2;
            assert(nco == chunk->n_detrend_f);
            assert(this->nt_chunk == ch_frb_io::constants::nt_per_assembled_chunk);
            for (ssize_t i=0; i<nco; i++)
                memcpy(chunk->detrend_params_f + i*this->nt_chunk,
                       coeffs + i*cstride,
                       this->nt_chunk * sizeof(float));
            chunk->detrend_f_type = "SPLINE";
            chunk->has_detrend_f = true;
            cout << "chime_spline_detrender: saved coefficients in assembled chunk " << chunk->ichunk << endl;
            chime_stream->updated_assembled_chunk(chunk);
        }
    }
}

Json::Value chime_spline_detrender::jsonize() const
{
    Json::Value ret = spline_detrender::jsonize();
    ret["class_name"] = "chime_spline_detrender";
    return ret;
}

shared_ptr<chime_spline_detrender> chime_spline_detrender::from_json(const Json::Value &j)
{
    // copied from spline_detrender!
    int nbins = int_from_json(j, "nbins");
    ssize_t nt_chunk = ssize_t_from_json(j, "nt_chunk");
    double epsilon = double_from_json(j, "epsilon");
    rf_kernels::axis_type axis = axis_type_from_json(j, "axis");
    return make_shared<chime_spline_detrender> (nt_chunk, axis, nbins, epsilon);
}

namespace {
    struct _init {
	_init() {
	    pipeline_object::register_json_deserializer("chime_spline_detrender", chime_spline_detrender::from_json);
	}
    } init;
}

/*
 // Externally callable factory function
 shared_ptr<wi_transform> make_spline_detrender(int nt_chunk, rf_kernels::axis_type axis, int nbins, double epsilon)
 {
 return make_shared<spline_detrender> (nt_chunk, axis, nbins, epsilon);
 }
 */



};
