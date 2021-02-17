#include "rf_pipelines_internals.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif


struct pipeline_fork : public pipeline_object
{
    struct element {
	// Initialized in constructor.
	string input_bufname;
	string output_bufname;
	
	// Initialized in bind().
	shared_ptr<ring_buffer> input_buffer;
	shared_ptr<ring_buffer> output_buffer;
	ssize_t csize = 0;
    };

    vector<element> elements;


    pipeline_fork(const vector<pair<string,string>> &bufnames) :
	pipeline_object("pipeline_fork")
    {
	unordered_set<string> all_names;

	for (const auto &p: bufnames) {
	    for (const string &s: { p.first, p.second }) {
		if (s.size() == 0)
		    _throw("empty bufname string was specified");
		if (all_names.count(s) > 0)
		    _throw("duplicate bufname '" + s + "' was specified");
	    }
	    
	    element e;
	    e.input_bufname = p.first;
	    e.output_bufname = p.second;

	    this->elements.push_back(e);
	}
    }


    virtual void _bind(ring_buffer_dict &rb_dict, Json::Value &json_attrs) override
    {
	this->nt_chunk_out = nt_chunk_in;
	this->nt_contig = nt_chunk_in;
	this->nt_maxgap = 0;

	for (element &e: this->elements) {
	    e.input_buffer = this->get_buffer(rb_dict, e.input_bufname);
	    e.output_buffer = this->create_buffer(rb_dict, e.output_bufname, e.input_buffer->cdims, e.input_buffer->nds);
	    e.csize = prod(e.input_buffer->cdims);
	}
    }


    virtual ssize_t _advance() override
    {
	// FIXME noticed in passing: this implementation looks wrong!  The pipeline_object is only
	// allowed to request up to 'nt_contig' elements in each ring_buffer_subarray, right?
	// See pipeline_spool.cpp::_advance() for an example.
	//
	// FIXME another problem: downsampling factor 'nds' ignored, right?

	for (element &e: this->elements) {
	    ring_buffer_subarray src(e.input_buffer, pos_lo, pos_hi, ring_buffer::ACCESS_READ);
	    ring_buffer_subarray dst(e.output_buffer, pos_lo, pos_hi, ring_buffer::ACCESS_APPEND);

	    for (ssize_t i = 0; i < e.csize; i++)
		memcpy(dst.data + i*dst.stride, src.data + i*src.stride, (pos_hi - pos_lo) * sizeof(float));
	}

	this->pos_lo = pos_hi.load();
	return SSIZE_MAX;
    }


    virtual Json::Value jsonize() const override
    {
	Json::Value ret;

	ret["class_name"] = "pipeline_fork";
	ret["bufnames"] = Json::Value(Json::arrayValue);

	for (const element &e: this->elements) {
	    Json::Value je(Json::arrayValue);
	    je.append(e.input_bufname);
	    je.append(e.output_bufname);

	    ret["bufnames"].append(je);
	}

	return ret;
    }

    static shared_ptr<pipeline_fork> from_json(const Json::Value &j)
    {
	const Json::Value &jb = array_from_json(j, "bufnames");
	vector<pair<string,string>> bufnames;

	for (int i = 0; i < int(jb.size()); i++) {
	    if (!jb[i].isArray() || (jb[i].size() != 2) || !jb[i][0].isString() || !jb[i][1].isString())
		throw runtime_error("pipeline_fork::from_json: expected each element of 'bufnames' array to be a pair of strings");

	    string input_bufname = jb[i][0].asString();
	    string output_bufname = jb[i][1].asString();
	    bufnames.push_back(pair<string,string> (input_bufname, output_bufname));
	}

	return make_shared<pipeline_fork> (bufnames);
    }
};


namespace {
    struct _init {
	_init() {
	    pipeline_object::register_json_deserializer("pipeline_fork", pipeline_fork::from_json);
	}
    } init;
}


// Externally callable factory function
shared_ptr<pipeline_object> make_pipeline_fork(const vector<pair<string,string>> &bufnames)
{
    return make_shared<pipeline_fork> (bufnames);
}


}  // namespace rf_pipelines
