#include "rf_pipelines_internals.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif


struct mask_expander : public chunked_pipeline_object
{
    // Initialized in constructor.
    // Reminder: 'nt_chunk' is inherited from base class.
    const double width;
    const double threshold;
    const double alpha;
    const string prev_wname;
    const rf_kernels::axis_type axis;

    // Initialized in _bindc().
    shared_ptr<ring_buffer> rb_prev_weights;
    shared_ptr<ring_buffer> rb_curr_weights;
    ssize_t nfreq = 0;

    // The following parameters are convenient in the kernel
    // Initialized in _bindc().

    float a;    // decay constant for moving average: exp(-1/(width*nfreq))
    float b;    // threshold for moving average: (1-threshold) / (1-a)
    float vp;   // value assigned to previously masked samples
    float vb;   // value assigned to boundary: vp/(1-a)

    // Temporary buffers, initialized in _allocate().
    uptr<float> tmp1;
    uptr<float> tmp2;


    mask_expander(rf_kernels::axis_type axis_, const string &prev_wname_, double width_, double threshold_, double alpha_ = 0.0, ssize_t nt_chunk_ = 0) :
	chunked_pipeline_object("mask_expander", false),   // can_be_first=false
	width(width_),
	threshold(threshold_),
	alpha(alpha_),
	prev_wname(prev_wname_),
	axis(axis_)
    {
	if (axis != rf_kernels::AXIS_FREQ)
	    _throw("only AXIS_FREQ is currently implemented");
	if (prev_wname.size() == 0)
	    _throw("'prev_wname' must be a nonempty string");
	if (prev_wname == "WEIGHTS")
	    _throw("'prev_wname' cannot be \"WEIGHTS\"");
	if ((width < 1.0e-5) || (width > 100.0))
	    _throw("'width' must be between 10^-5 and 100");
	if ((threshold <= 0.0) || (threshold >= 1.0))
	    _throw("'threshold' must be between 0 and 1");
	if ((alpha < -1.0) || (alpha > 1.0))
	    _throw("'alpha' must be between -1 and 1");

	this->nt_chunk = nt_chunk_;

	stringstream ss;

	ss << "mask_expander(" << rf_kernels::axis_type_to_string(axis)
	   << ",prev_wname='" << prev_wname 
	   << "',width=" << width 
	   << ",threshold=" << threshold;

	if (alpha != 0.0)
	    ss << ",alpha=" << alpha;

	if (nt_chunk_ != 0)
	    ss << ",nt_chunk=" << nt_chunk_;

	ss << ")";

	this->name = ss.str();
    }


    // Called after this->nfreq is initialized.
    virtual void _bindc(ring_buffer_dict &rb_dict, Json::Value &json_attrs) override
    {
	this->rb_prev_weights = this->get_buffer(rb_dict, prev_wname);
	this->rb_curr_weights = this->get_buffer(rb_dict, "WEIGHTS");

	if (rb_prev_weights->cdims != rb_curr_weights->cdims)
	    _throw("'intensity' and 'weights' buffers have different dimensions");
	if (rb_prev_weights->nds != rb_curr_weights->nds)
	    _throw("'intensity' and 'weights' buffers have different downsampling");
	if (rb_curr_weights->cdims.size() != 1)
	    _throw("expected intensity/weights arrays to be two-dimensional");

	double t = this->threshold;
	double pm = this->alpha;

	this->nfreq = rb_curr_weights->cdims[0];
	this->a = exp(-1.0 / (width * nfreq));
	this->b = (1-t) / (1-a);
	this->vp = (pm >= 0.0) ? (1-t+pm*t) : ((1+pm) * (1-t));
	this->vb = vp / (1-a);
    }


    virtual void _allocate() override
    {
	this->tmp1 = make_uptr<float> (nfreq);
	this->tmp2 = make_uptr<float> (nfreq);
    }


    virtual bool _process_chunk(ssize_t pos) override
    {
	ring_buffer_subarray wprev(rb_prev_weights, pos, pos + nt_chunk, ring_buffer::ACCESS_READ);
	ring_buffer_subarray wcurr(rb_curr_weights, pos, pos + nt_chunk, ring_buffer::ACCESS_RW);

	for (int it = 0; it < nt_chunk; it++) {
	    float t = vb;

	    for (int ifreq = 0; ifreq < nfreq; ifreq++) {
		if (wcurr.data[ifreq * wcurr.stride + it] > 0.0)
		    tmp1[ifreq] = 1.0;  // currently unmasked
		else if (wprev.data[ifreq * wprev.stride + it] > 0.0)
		    tmp1[ifreq] = 0.0;  // currently masked, previously unmasked
		else
		    tmp1[ifreq] = vp;   // previously (and currently) masked

		t = tmp1[ifreq] + a*t;
		tmp2[ifreq] = t;
	    }

	    t = vb;
	    for (int ifreq = nfreq-1; ifreq >= 0; ifreq--) {
		t = tmp1[ifreq] + a*t;
		if ((t < b) && (tmp2[ifreq] < b))
		    wcurr.data[ifreq*wprev.stride + it] = 0.0;
	    }
	}

	return true;
    }


    virtual void _deallocate() override
    {
	this->tmp1.reset();
	this->tmp2.reset();
    }


    virtual void _unbindc() override
    {
	this->rb_prev_weights.reset();
	this->rb_curr_weights.reset();
    }


    virtual Json::Value jsonize() const override
    {
	Json::Value ret;

	ret["class_name"] = "mask_expander";
	ret["width"] = this->width;
	ret["threshold"] = this->threshold;
	ret["alpha"] = this->alpha;
	ret["prev_wname"] = this->prev_wname;
	ret["nt_chunk"] = int(this->get_prebind_nt_chunk());
	ret["axis"] = rf_kernels::axis_type_to_string(rf_kernels::AXIS_FREQ);

	return ret;
    }


    static shared_ptr<mask_expander> from_json(const Json::Value &j)
    {
	double width = double_from_json(j, "width");
	double threshold = double_from_json(j, "threshold");
	double alpha = double_from_json(j, "alpha");
	string prev_wname = string_from_json(j, "prev_wname");
	ssize_t nt_chunk = int_from_json(j, "nt_chunk");
	rf_kernels::axis_type axis = axis_type_from_json(j, "axis");

	return make_shared<mask_expander> (axis, prev_wname, width, threshold, alpha, nt_chunk);
    }
};


namespace {
    struct _init {
	_init() {
	    pipeline_object::register_json_deserializer("mask_expander", mask_expander::from_json);
	}
    } init;
}


// Externally callable factory function
shared_ptr<pipeline_object> make_mask_expander(rf_kernels::axis_type axis, const string &prev_wname, double width, double threshold, double alpha, ssize_t nt_chunk)
{
    return make_shared<mask_expander> (axis, prev_wname, width, threshold, alpha, nt_chunk);
}


}  // namespace rf_pipelines
