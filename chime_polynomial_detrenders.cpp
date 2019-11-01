#include "rf_kernels/spline_detrender.hpp"
#include "rf_pipelines_internals.hpp"
#include "ch_frb_io.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif


class chime_polynomial_detrender : public polynomial_detrender, public chime_wi_transform {
public:
    chime_polynomial_detrender(rf_kernels::axis_type axis_, int nt_chunk_, int nbins_, double epsilon_);
    virtual void _handle_coeffs(ssize_t pos, float *coeffs, float *intensity, ssize_t istride, float *weights, ssize_t wstride) override;
    virtual void _bind_transform(Json::Value &json_attrs) override;
    virtual Json::Value jsonize() const override;
    static std::shared_ptr<chime_polynomial_detrender> from_json(const Json::Value &j);
};


chime_polynomial_detrender::chime_polynomial_detrender(rf_kernels::axis_type axis_, int nt_chunk_,  int ndeg_, double epsilon_) :
    wi_transform("chime_polynomial_detrender"),
    polynomial_detrender(axis_, nt_chunk_, ndeg_, epsilon_)
{}

void chime_polynomial_detrender::_bind_transform(Json::Value &json_attrs) {
    // borrrring
    polynomial_detrender::_bind_transform(json_attrs);
    chime_wi_transform::_bind_transform(json_attrs);
}

void chime_polynomial_detrender::_handle_coeffs(ssize_t pos, float *coeffs, float *intensity, ssize_t istride, float *weights, ssize_t wstride) {
    if (this->chime_stream && kernel.axis == rf_kernels::AXIS_TIME) {
        auto chunk = this->assembled_chunk_for_pos(pos);
        if (chunk) {
            assert(!chunk->has_detrend_t);
            assert(chunk->detrend_params_t);
            int nco = (this->kernel.polydeg+1);
            assert(nco == chunk->n_detrend_t);
            assert(this->nt_chunk == ch_frb_io::constants::nt_per_assembled_chunk);
            memcpy(chunk->detrend_params_t,
                   coeffs,
                   this->nfreq * nco * sizeof(float));
            chunk->detrend_t_type = "POLYNOMIAL";
            chunk->has_detrend_t = true;
            cout << "chime_polynomial_detrender: saved coefficients in assembled chunk " << chunk->ichunk << endl;
            chime_stream->updated_assembled_chunk(chunk);
        }
    }
}

Json::Value chime_polynomial_detrender::jsonize() const
{
    Json::Value ret = polynomial_detrender::jsonize();
    ret["class_name"] = "chime_polynomial_detrender";
    return ret;
}

shared_ptr<chime_polynomial_detrender> chime_polynomial_detrender::from_json(const Json::Value &x)
{
    // copied from polynomial_detrender!
    rf_kernels::axis_type axis = axis_type_from_json(x, "axis");
    int nt_chunk = int_from_json(x, "nt_chunk");
    double polydeg = double_from_json(x, "polydeg");
    double epsilon = double_from_json(x, "epsilon");
    return make_shared<chime_polynomial_detrender>(axis, nt_chunk, polydeg, epsilon);
}

namespace {
    struct _init {
	_init() {
	    pipeline_object::register_json_deserializer("chime_polynomial_detrender", chime_polynomial_detrender::from_json);
	}
    } init;
}

/*
 // Externally callable factory function
 shared_ptr<wi_transform> make_polynomial_detrender(int nt_chunk, rf_kernels::axis_type axis, int nbins, double epsilon)
 {
 return make_shared<polynomial_detrender> (nt_chunk, axis, nbins, epsilon);
 }
 */



};
