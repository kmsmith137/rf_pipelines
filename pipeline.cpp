#include "rf_pipelines_internals.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}  // namespace rf_pipelines
#endif


pipeline::pipeline(const string &name_) :
    pipeline_object(name_)
{
    this->_update_name();
}


pipeline::pipeline(const vector<shared_ptr<pipeline_object>> &elements_, const string &name_) :
    pipeline_object(name_),
    elements(elements_)
{
    for (const auto &p: elements)
	if (p.get() == nullptr)
	    _throw("null pointer in pipeline constructor");
    
    this->_update_name();
}


void pipeline::add(const shared_ptr<pipeline_object> &p)
{
    if (p.get() == nullptr)
	_throw("null pointer in pipeline constructor");
    if (this->is_bound())
	_throw("pipeline::add() was called after bind()");
    
    elements.push_back(p);    
    this->_update_name();
}


void pipeline::_update_name()
{
    if (name == _default_name)
	name = _default_name = "pipeline<" + to_string(elements.size()) + ">";
}


void pipeline::_bind(ring_buffer_dict &rb_dict, Json::Value &json_data)
{
    if (elements.size() == 0)
	_throw("pipeline is empty (length-zero)");

    // Initial values, updated in loop below
    this->nt_chunk_out = nt_chunk_in;
    this->nt_maxgap = 0;
    this->nt_contig = 1;

    for (auto &p: this->elements) {
	p->bind(rb_dict, nt_chunk_out, nt_maxlag + nt_maxgap, json_data);
	this->nt_chunk_out = p->nt_chunk_out;
	this->nt_maxgap += p->nt_maxgap;
    }
}


ssize_t pipeline::_advance()
{
    // Initial values, updated in loop below
    ssize_t ret = SSIZE_MAX;
    this->pos_lo = pos_hi;

    for (auto &p: this->elements) {
	ssize_t n = p->advance(pos_lo, pos_max);
	ret = min(ret, n);
	pos_lo = p->pos_lo;
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

    if (name != _default_name)
	ret["name"] = name;

    for (size_t i = 0; i < elements.size(); i++)
	ret["elements"].append(elements[i]->jsonize());

    return ret;
}


shared_ptr<pipeline> pipeline::from_json(const Json::Value &x)
{
    if (string_from_json(x,"class_name") != "pipeline")
	throw runtime_error("rf_pipelines: expected class_name=\"pipeline\" in pipeline json constructor");

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

void pipeline::_start_pipeline(Json::Value &j)
{
    for (auto &p: this->elements)
	p->start_pipeline(this->out_mp, j);
}

void pipeline::_end_pipeline(Json::Value &j)
{
    for (int i = 0; i < (int)elements.size(); i++) {
	j["pipeline"].append(Json::Value(Json::objectValue));
	elements[i]->end_pipeline(j["pipeline"][i]);
    }
}


namespace {
    struct _init {
	_init() {
	    pipeline_object::register_json_constructor("pipeline", pipeline::from_json);
	}
    } init;
}


}  // namespace rf_pipelines
