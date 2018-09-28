#include "rf_pipelines_internals.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}  // namespace rf_pipelines
#endif


pipeline::pipeline(const string &name_) :
    pipeline_object("pipeline", name_)
{ }


pipeline::pipeline(const vector<shared_ptr<pipeline_object>> &elements_, const string &name_) :
    pipeline_object("pipeline", name_),
    elements(elements_)
{
    for (const auto &p: elements)
	if (p.get() == nullptr)
	    _throw("null pointer in pipeline constructor");
}

// For subclasses (e.g. wi_sub_pipeline)
pipeline::pipeline(const string &class_name_, const string &name_) :
    pipeline_object(class_name_, name_)
{ }

void pipeline::add(const shared_ptr<pipeline_object> &p)
{
    if (p.get() == nullptr)
	_throw("null pointer in pipeline constructor");
    if (this->state != UNBOUND)
	_throw("pipeline::add() was called after bind()");
    
    elements.push_back(p);    
}


void pipeline::_bind(ring_buffer_dict &rb_dict, Json::Value &json_attrs)
{
    if (elements.size() == 0)
	_throw("pipeline is empty (length-zero)");

    // Initial values, updated in loop below
    this->nt_chunk_out = nt_chunk_in;
    this->nt_maxgap = 0;
    this->nt_contig = 1;

    for (size_t i = 0; i < elements.size(); i++) {
	auto p = elements[i];

	run_params params = this->get_params();
	params.container_depth++;
	params.container_index = i;

	elements[i]->bind(params, rb_dict, nt_chunk_out, nt_maxlag + nt_maxgap, json_attrs, this->out_mp);
	this->nt_chunk_out = p->nt_chunk_out;
	this->nt_maxgap += p->nt_maxgap;
    }
}


ssize_t pipeline::_advance()
{
    // Initial values, updated in loop below
    ssize_t ret = SSIZE_MAX;
    this->pos_lo = pos_hi;

    for (size_t i = 0; i < elements.size(); i++) {
	auto p = elements[i];

	ssize_t m = p->pos_lo;  // only needed for debug print
	ssize_t n = p->advance(pos_lo, pos_max);
	ret = min(ret, n);
	pos_lo = p->pos_lo;

	if (_params.noisy())
	    p->_print("advance " + to_string(m) + " -> " + to_string(pos_lo));
    }
    
    return ret;
}


ssize_t pipeline::get_preferred_chunk_size()
{
    if (elements.size() == 0)
	_throw("pipeline is empty (length-zero)");
    
    return elements[0]->get_preferred_chunk_size();
}


Json::Value pipeline::jsonize() const
{
    if (elements.size() == 0)
	_throw("pipeline is empty (length-zero)");
    
    Json::Value ret;
    ret["class_name"] = "pipeline";
    ret["name"] = name;

    for (size_t i = 0; i < elements.size(); i++)
	ret["elements"].append(elements[i]->jsonize());

    return ret;
}


shared_ptr<pipeline> pipeline::from_json(const Json::Value &x)
{
    if (string_from_json(x,"class_name") != "pipeline")
	throw runtime_error("rf_pipelines: expected class_name=\"pipeline\" in pipeline json constructor");

    // Old versions of jsonize() didn't always include the 'name'.
    string name = x.isMember("name") ? string_from_json(x, "name") : "";

    vector<shared_ptr<pipeline_object>> elements;
    for (const auto &s: array_from_json(x, "elements"))
	elements.push_back(pipeline_object::from_json(s));

    return make_shared<pipeline> (elements, name);
}


void pipeline::_allocate()
{
    for (auto &p: this->elements)
	p->allocate();
}

void pipeline::_deallocate()
{
    for (auto &p: this->elements)
	p->deallocate();
}

void pipeline::_start_pipeline(Json::Value &json_attrs)
{
    for (auto &p: this->elements)
	p->start_pipeline(json_attrs);
}

void pipeline::_end_pipeline(Json::Value &json_output)
{
    for (int i = 0; i < (int)elements.size(); i++) {
	json_output["pipeline"].append(Json::Value(Json::objectValue));
	elements[i]->end_pipeline(json_output["pipeline"][i]);
    }
}

void pipeline::_reset()
{
    for (auto &p: this->elements)
	p->reset();
}

void pipeline::_unbind()
{
    for (auto &p: this->elements)
	p->unbind();
}

void pipeline::_get_info(Json::Value &j)
{
    j["pipeline"] = Json::Value(Json::arrayValue);
    double mb_cumul = 0.0;

    for (auto &p: this->elements) {
	Json::Value jr = p->get_info();
	double mb = double_from_json(jr, "mb_cumul");

	j["pipeline"].append(jr);
	mb_cumul += mb;
    }

    j["mb_cumul"] = mb_cumul;
}

void pipeline::visit_pipeline(std::function<void(pipeline_object*,int)> f, int depth)
{
    f(this, depth);
    
    for (auto &p: this->elements)
	p->visit_pipeline(f, depth+1);
}


namespace {
    struct _init {
	_init() {
	    pipeline_object::register_json_deserializer("pipeline", pipeline::from_json);
	}
    } init;
}


}  // namespace rf_pipelines
