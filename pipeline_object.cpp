#include <fstream>
#include <algorithm>
#include "rf_pipelines_internals.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}  // emacs pacifier
#endif


// Global registry (class_name -> json_deserializer).
//
// Note: We use a pointer here, rather than simply declaring a static global unordered_map<...>, 
// because we initialize entries with other static global declarations (see example at the end
// of pipeline.cpp), and there is no way to ensure that the unordered_map<> constructor call
// occurs first.  
//
// In contrast, if we use a pointer, then the first call to register_json_deserializer() will initialize 
// the regsistry (see code below).  However, this means that all users of the registry (e.g. from_json(),
// _show_registered_json_deserializers()) must handle the corner case where no constructors have been
// registered yet, and the pointer is still NULL.

using json_registry_t = unordered_map<string, pipeline_object::json_deserializer_t>;
static json_registry_t *json_registry = nullptr;   // global


pipeline_object::pipeline_object(const string &class_name_, const string &name_) : 
    state(UNBOUND),
    class_name(class_name_),
    name((name_.size() > 0) ? name_ : class_name_)
{
    if (class_name.size() == 0)
	throw runtime_error("rf_pipelines::pipeline_object constructor: 'class_name' must be a nonempty string");
}

pipeline_object::~pipeline_object() { }


void pipeline_object::_throw(const string &msg) const
{
    string prefix = (name.size() > 0) ? ("rf_pipelines: " + name + ": ") : "rf_pipelines: ";
    throw runtime_error(prefix + msg);
}


void pipeline_object::_print(const string &msg) const
{
    if (state != UNBOUND) {
	cout << string(_params.container_depth*4, ' ');

	if (_params.container_index >= 0)
	    cout << "[" << _params.container_index << "]: ";
    }

    cout << msg << ": " << this->name << "\n";
}


// -------------------------------------------------------------------------------------------------
//
// bind(), unbind().


// This version of bind() is only called on a top-level pipeline.
void pipeline_object::bind(const run_params &params)
{
    if (params.verbosity >= 2)
	cout << "rf_pipelines: bind() called\n";

    if (this->state >= BOUND) {
	string s = _params.mismatch(params);
	if (s.size() == 0)
	    return;
	_throw("double call to bind() with mismatched " + s + ".  Maybe you forgot to call unbind()?");
    }

    ssize_t n = this->get_preferred_chunk_size();
    if (n <= 0)
	_throw("this object cannot be first in pipeline, since its get_preferred_chunk_size() method is undefined");

    this->json_attrs1 = Json::Value(Json::objectValue);

    ring_buffer_dict rb_dict;
    shared_ptr<outdir_manager> mp = make_shared<outdir_manager> (params.outdir, params.clobber);
    
    this->bind(params, rb_dict, n, n, json_attrs1, mp);
}

// The non-virtual function bind() wraps the pure virtual function _bind().
void pipeline_object::bind(const run_params &params, ring_buffer_dict &rb_dict, ssize_t nt_chunk_in_, ssize_t nt_maxlag_, Json::Value &json_attrs, const shared_ptr<outdir_manager> &mp)
{
    rf_assert(nt_chunk_in_ > 0);
    rf_assert(nt_maxlag_ > 0);
    rf_assert(json_attrs.isObject());
    rf_assert(mp);

    params.check();

    if (name.size() == 0)
	throw runtime_error("rf_pipelines: pipeline_object did not initialize 'name' field in its constructor");

    if (this->state != UNBOUND)
	_throw("Double call to pipeline_object::bind(), suspect pipeline_object is reused in a pipeline, or used in multiple pipelines.");

    this->_params = params;
    this->nt_chunk_in = nt_chunk_in_;
    this->nt_maxlag = nt_maxlag_;
    this->out_mp = mp;
    this->state = BINDING;

    this->_bind(rb_dict, json_attrs);

    rf_assert(nt_chunk_in == nt_chunk_in_);
    rf_assert(nt_maxlag == nt_maxlag_);

    if (nt_maxgap < 0)
	_throw("_bind() failed to initialize nt_maxgap");
    if (nt_chunk_out <= 0)
	_throw("_bind() failed to initialize nt_chunk_out");
    if (nt_contig <= 0)
	_throw("_bind() failed to initialize nt_contig");

    for (auto &rb: this->all_ring_buffers)
	rb->update_params(nt_contig, nt_maxlag + nt_maxgap);
    for (auto &zt: this->zoomable_tilesets)
	zt->update_params(nt_contig, nt_maxlag + nt_maxgap);

    this->state = BOUND;
}


const run_params &pipeline_object::get_params() const
{
    if (state == UNBOUND)
	_throw("get_params() called on unbound pipeline_object");
    
    return this->_params;
}


// Should be called from _bind().
shared_ptr<ring_buffer> pipeline_object::get_buffer(ring_buffer_dict &rb_dict, const string &bufname)
{
    if (this->state != BINDING)
	_throw("pipeline_object::get_buffer() called outside bind()");
    if (!has_key(rb_dict, bufname))
	_throw("buffer '" + bufname + "' does not exist in pipeline");
    if (_params.verbosity >= 3)
	cout << "    bind(): get_buffer(" << bufname << "): " << this->name << "\n";

    auto ret = rb_dict[bufname];
    all_ring_buffers.push_back(ret);

    return ret;
}


// Should be called from _bind().
shared_ptr<ring_buffer> pipeline_object::create_buffer(ring_buffer_dict &rb_dict, const string &bufname, const vector<ssize_t> &cdims, ssize_t nds)
{
    if (this->state != BINDING)
	_throw("pipeline_object::create_buffer() called outside bind()");
    if (has_key(rb_dict, bufname))
	_throw("buffer '" + bufname + "' already exists in pipeline");
    if (_params.verbosity >= 3)
	cout << "    bind(): create_buffer(" << bufname << "): " << this->name << "\n";

    auto ret = make_shared<ring_buffer> (cdims, nds, _params.debug, bufname);

    rb_dict[bufname] = ret;
    all_ring_buffers.push_back(ret);
    new_ring_buffers.push_back(ret);
    
    return ret;
}


// Called from _bind().
vector<shared_ptr<ring_buffer>> 
pipeline_object::add_zoomable_tileset(const shared_ptr<zoomable_tileset> &zt, const Json::Value &json_attrs)
{
    if (this->state != BINDING)
	_throw("add_zoomable_tileset() called outside bind()");
    if (!this->out_mp)
	_throw("add_zoomable_tileset() internal error: 'out_mp' is an empty pointer.");
    if (out_mp->outdir.size() == 0)
	_throw("add_zoomable_tileset() should not be called if there is no pipeline outdir");
    if (_params.verbosity >= 3)
	cout << "    bind(): add_zoomable_tileset(" << zt->img_prefix << "): " << this->name << "\n";

    auto ret = make_shared<zoomable_tileset_state> (zt, *this);

    this->zoomable_tilesets.push_back(ret);
    return ret->ring_buffers[0];
}


void pipeline_object::unbind()
{
    if (this->state == UNBOUND)
	return;

    if (this->state >= ALLOCATED)
	this->deallocate();

    rf_assert(state == BOUND);

    this->_unbind();

    this->nt_chunk_in = 0;
    this->nt_maxlag = 0;
    this->nt_chunk_out = 0;
    this->nt_maxgap = -1;
    this->nt_contig = 0;

    this->all_ring_buffers.clear();
    this->new_ring_buffers.clear();
    this->zoomable_tilesets.clear();
    this->json_attrs1 = Json::Value();
    this->json_attrs2 = Json::Value();
    this->out_mp.reset();

    this->state = UNBOUND;
}


// Default virtual (see comment in rf_pipelines.hpp for discussion)
ssize_t pipeline_object::get_preferred_chunk_size()
{
    return 0;
}

// Default virtual
void pipeline_object::_unbind() { }


// -------------------------------------------------------------------------------------------------
//
// allocate(), deallocate()


void pipeline_object::allocate()
{
    if (this->state >= ALLOCATED)
	return;
    if (this->state < BOUND)
	_throw("allocate() called before bind()");
    if ((_params.verbosity >= 2) && (_params.container_depth == 0))
	cout << "rf_pipelines: allocate() called\n";

    for (auto &p: this->new_ring_buffers)
	p->allocate();
    for (auto &p: this->zoomable_tilesets)
	p->allocate();

    this->_allocate();
    this->state = ALLOCATED;
}


void pipeline_object::deallocate()
{
    if (this->state < ALLOCATED)
	return;

    if (this->state >= RUNNING)
	this->reset();

    rf_assert(state == ALLOCATED);
	
    this->_deallocate();
    
    for (auto &p: this->new_ring_buffers)
	p->deallocate();
    for (auto &p: this->zoomable_tilesets)
	p->deallocate();

    this->state = BOUND;
}


// Default virtuals do nothing.
void pipeline_object::_allocate() { }
void pipeline_object::_deallocate() { }
    

// -------------------------------------------------------------------------------------------------
//
// run() and friends


Json::Value pipeline_object::run(const run_params &params, const callback_t &callback)
{
    params.check();

    if (this->state >= RUNNING)
	_throw("Double call to pipeline_object::run(), maybe you're missing a call to reset(), deallocate(), or unbind()?");

    this->bind(params);
    this->allocate();

    this->json_attrs2 = Json::Value(Json::objectValue);
    this->start_pipeline(json_attrs2);

    rf_assert(this->state == RUNNING);

    // We wrap the advance() loop in try..except, so that if an exception is thrown, we
    // still call end_pipeline() to clean up, and write partially complete output files.

    if (params.verbosity >= 2)
	cout << "rf_pipelines: entering main advance loop\n";

    bool exception_thrown = false;
    string exception_text;

    try {
	ssize_t nt_end = SSIZE_MAX;
    
	while (this->pos_lo < nt_end) {
	    if (params.verbosity >= 3)
		cout << "rf_pipelines: main advance loop " << pos_hi << " -> " << (pos_hi + nt_chunk_in) << "\n";

	    ssize_t m = pos_hi + nt_chunk_in;
	    ssize_t n = this->advance(m, m);
	    nt_end = min(nt_end, n);

	    if (callback)
		callback(pos_lo, pos_hi);
	}
    } catch (std::exception &e) {
	exception_text = e.what();
	exception_thrown = true;
	// fall through...
    }

    if (params.verbosity >= 2) {
	string s = exception_thrown ? "exception thrown" : "normal termination";
	cout << "rf_pipelines: exiting advance() loop, pos=" << pos_lo << ", " << s << "\n";
    }

    Json::Value json_output(Json::objectValue);
    json_output["success"] = !exception_thrown;
    json_output["error_message"] = exception_text;

    // I decided to put all run_params in the json output, even
    // those whose usefulness is dubious!
    json_output["outdir"] = _params.outdir;
    json_output["clobber"] = _params.clobber;
    json_output["img_nzoom"] = Json::Int64(_params.img_nzoom);
    json_output["img_nds"] = Json::Int64(_params.img_nds);
    json_output["img_nx"] = Json::Int64(_params.img_nx);
    json_output["verbosity"] = _params.verbosity;
    json_output["debug"] = _params.debug;

    add_json_object(json_output, this->json_attrs1);  // bind() attributes
    add_json_object(json_output, this->json_attrs2);  // start_pipeline() attributes

    if (params.verbosity >= 2)
	cout << "rf_pipelines: calling end_pipeline()\n";

    this->end_pipeline(json_output);

    // Try to write json file, even if exception was thrown.
    if (params.outdir.size() > 0) {
	string json_filename = params.outdir + "/rf_pipeline_0.json";
	ofstream f(json_filename);

	if (f.fail())
	    _throw("couldn't open output file " + json_filename);

	Json::StyledWriter w;
	f << w.write(json_output);

	if (params.verbosity >= 2)
	    cout << "wrote " << json_filename << endl;
    }

    // FIXME add boolean flag to deallocate on pipeline exit.

    // FIXME is there a better way to save and re-throw the exception?
    // The naive approach, saving a copy of the std::exception, doesn't preserve the exception_text.

    if (exception_thrown)
	throw runtime_error(exception_text);

    return json_output;
}


// The non-virtual function advance() wraps the pure virtual function _advance().
ssize_t pipeline_object::advance(ssize_t pos_hi_, ssize_t pos_max_)
{
    struct timeval tv0 = get_time();

    rf_assert(state == RUNNING);
    rf_assert(nt_chunk_in > 0);
    rf_assert(nt_chunk_out > 0);
    
    rf_assert(pos_hi <= pos_hi_);
    rf_assert(pos_hi_ <= pos_max_);
    rf_assert(pos_max_ <= pos_hi + nt_maxlag);
    rf_assert(pos_hi_ % nt_chunk_in == 0);

    this->pos_hi = pos_hi_;
    this->pos_max = pos_max_;    

    ssize_t ret = this->_advance();

    if (pos_hi != pos_hi_)
	_throw("internal error: value of pos_hi was modified in advance()");
    if (pos_lo % nt_chunk_out)
	_throw("internal error: pos_lo is not a multiple of nt_chunk_out after advance()");	
    if (pos_lo > pos_hi)
	_throw("internal error: pos_lo > pos_hi after advance()");
    if (pos_hi - pos_lo > nt_maxgap)
	_throw("internal error: (pos_hi-pos_lo) > nt_maxgap after advance().");

    for (const auto &p: this->zoomable_tilesets)
	p->advance(this->pos_lo);

    this->time_spent_in_transform += time_diff(tv0, get_time());

    return ret;
}


void pipeline_object::start_pipeline(Json::Value &json_attrs)
{
    rf_assert(this->state == ALLOCATED);

    if ((_params.verbosity >= 2) && (_params.container_depth == 0))
	cout << "rf_pipelines: start_pipeline() called\n";

    this->plot_groups.clear();
    this->time_spent_in_transform = 0.0;
    
    this->pos_lo = 0;
    this->pos_hi = 0;
    this->pos_max = 0;
    
    for (auto &p: this->all_ring_buffers)
	p->reset();

    this->_start_pipeline(json_attrs);

    this->state = RUNNING;
}


void pipeline_object::end_pipeline(Json::Value &json_output)
{
    rf_assert(this->state == RUNNING);

    if ((_params.verbosity >= 2) && (_params.container_depth == 0))
	cout << "rf_pipelines: end_pipeline() called\n";

    if (!json_output.isObject())
	_throw("end_pipeline(): internal error: Json::Value was not an Object as expected");

    this->_end_pipeline(json_output);

    for (const auto &p: this->zoomable_tilesets)
	p->flush();

    this->state = DONE;

    if (!json_output.isMember("class_name"))
	json_output["class_name"] = this->class_name;
    if (!json_output.isMember("name"))
	json_output["name"] = this->name;
    if (!json_output.isMember("cpu_time"))
	json_output["cpu_time"] = this->time_spent_in_transform;

    if (plot_groups.size() > 0) {
	for (const auto &g: plot_groups) {
	    if (g.is_empty)
		continue;
	    
	    Json::Value jp;
            jp["name"] = g.name;
            jp["nt_per_pix"] = g.nt_per_pix;
            jp["ny"] = g.ny;
            jp["it0"] = Json::Value::Int64(g.curr_it0);
            jp["it1"] = Json::Value::Int64(g.curr_it1);
            jp["files"].append(g.files);

            json_output["plots"].append(jp);
	}
    }

    for (const auto &p: this->zoomable_tilesets) {
	// p->json_output is a JSON array.
	for (const auto &j: p->json_output)
	    json_output["plots"].append(j);
    }
}


void pipeline_object::reset()
{
    if (this->state <= ALLOCATED)
	return;

    this->pos_lo = 0;
    this->pos_hi = 0;
    this->pos_max = 0;
    this->plot_groups.clear();
    this->time_spent_in_transform = 0.0;
    this->json_attrs2 = Json::Value();
    
    for (auto &p: this->all_ring_buffers)
	p->reset();
    for (auto &p: this->zoomable_tilesets)
	p->reset();

    this->_reset();
    this->state = ALLOCATED;
}


Json::Value pipeline_object::get_info() 
{
    if (state < BINDING)
	_throw("rf_pipelines::pipeline_object::get_info() must be called after bind()\n");

    Json::Value j(Json::objectValue);
    double mb = 0.0;

    j["class_name"] = this->class_name;
    j["name"] = this->name;
    j["nt_chunk_in"] = Json::Int64(this->nt_chunk_in);
    j["nt_maxlag"] = Json::Int64(this->nt_maxlag);
    j["nt_chunk_out"] = Json::Int64(this->nt_chunk_out);
    j["nt_maxgap"] = Json::Int64(this->nt_maxgap);
    j["nt_contig"] = Json::Int64(this->nt_contig);
    j["ring_buffers"] = Json::Value(Json::arrayValue);

    for (auto &p: this->new_ring_buffers) {
	Json::Value jr = p->get_info();
	j["ring_buffers"].append(jr);
	mb += double_from_json(jr, "mb");
    }

    j["mb_local"] = mb;
    j["mb_cumul"] = mb;
    
    this->_get_info(j);

    return j;
}


// Default virtuals
void pipeline_object::_start_pipeline(Json::Value &j) { }
void pipeline_object::_end_pipeline(Json::Value &j) { }
void pipeline_object::_get_info(Json::Value &j) { }
void pipeline_object::_reset() { }


// -------------------------------------------------------------------------------------------------
//
// json serialization/deserialization


// default virtual
Json::Value pipeline_object::jsonize() const
{
    _throw("jsonize() not implemented");
    return Json::Value();  // compiler pacifier
}


// static member function
void pipeline_object::register_json_deserializer(const string &class_name, const json_deserializer_t &f)
{
    if (class_name.size() == 0)
	throw runtime_error("rf_pipelines::pipeline_object::register_json_deserializer(): class_name must be a nonempty string");

    // First call to register_json_deserializer() will initialize the regsistry.
    if (!json_registry)
	json_registry = new json_registry_t;

    auto p = json_registry->find(class_name);

    if (p != json_registry->end())
	throw runtime_error("rf_pipelines::pipeline_object::register_json_deserializer(): duplicate registration for class_name='" + class_name + "'");

    (*json_registry)[class_name] = f;
}


// Static member function, for debugging.
void pipeline_object::_show_registered_json_deserializers()
{
    // All users of the json registry must handle the corner case where no constructors have
    // been registered yet, and the pointer is still NULL,

    if (!json_registry)
	return;

    vector<string> all_class_names;

    for (const auto &p: *json_registry)
	all_class_names.push_back(p.first);

    std::sort(all_class_names.begin(), all_class_names.end());
    
    cout << "[";
    for (const string &class_name : all_class_names)
	cout << " " << class_name;
    cout << " ]";
}


// static member function
pipeline_object::json_deserializer_t pipeline_object::_find_json_deserializer(const string &class_name)
{
    // All users of the json registry must handle the corner case where no constructors have
    // been registered yet, and the pointer is still NULL,

    if (!json_registry)
	return NULL;

    auto p = json_registry->find(class_name);
    return (p != json_registry->end()) ? p->second : NULL;
}


// static member function
shared_ptr<pipeline_object> pipeline_object::from_json(const Json::Value &x)
{
    if (!x.isObject())
	throw runtime_error("rf_pipelines: pipeline_object::from_json(): expected json argument to be an Object");

    // throws exception if 'class_name' not found
    string class_name = string_from_json(x, "class_name");

    json_deserializer_t f = _find_json_deserializer(class_name);

    if (f == NULL)
	throw runtime_error("rf_pipelines::pipeline_object::from_json(): class_name='" + class_name + "' not found.  Maybe you're trying to use a python transform in C++ code?");

    shared_ptr<pipeline_object> ret = f(x);

    if (!ret)
	throw runtime_error("rf_pipelines::pipeline_object::from_json(): json_deserializer for class_name='" + class_name + "' returned empty pointer");

    return ret;
}


// -------------------------------------------------------------------------------------------------
//
// Output file management (including plots)


// Returns group id
int pipeline_object::add_plot_group(const string &name, int nt_per_pix, int ny)
{
    if ((this->state != ALLOCATED) || (!this->out_mp))
	_throw("add_plot_group() must be called from _start_pipeline()");
    if (this->out_mp->outdir.size() == 0)
	_throw("add_plot_group() was called in pipeline with no outdir (probably need to check out_mp->outdir.size() before calling add_plot_group)");
    if (nt_per_pix < 1)
	_throw("add_plot_group(): nt_per_pix must be >= 1");
    if (ny < 1)
	_throw("add_plot_group(): ny must be >= 1");

    plot_group g;
    g.name = name;
    g.nt_per_pix = nt_per_pix;
    g.ny = ny;

    this->plot_groups.push_back(g);
    return plot_groups.size()-1;
}


string pipeline_object::add_plot(const string &basename, int64_t it0, int nt, int nx, int ny, int group_id)
{
    if (this->state != RUNNING)
	_throw("add_plot() must be called from _advance()");
    if (this->plot_groups.size() == 0)
	_throw("add_plot() called but no plot_groups defined, maybe you forgot to call add_plot_group()?");

    if ((group_id < 0) || (group_id >= (int)plot_groups.size()))
	_throw("add_plot(): bad group_id specified");

    plot_group &g = plot_groups[group_id];

    if (nt != g.nt_per_pix * nx)
	_throw("add_plot(): requirement (nt == nx*nt_per_pix) failed");
    if (ny != g.ny)
	_throw("add_plot(): ny doesn't match value specified in add_plot_group()");

    if (g.is_empty) {
	g.is_empty = false;
	g.curr_it0 = it0;
    }
    else if (it0 != g.curr_it1)
	_throw("add_plot(): plot time ranges are not contiguous");

    string filename = this->add_file(basename);

    Json::Value file;
    file["filename"] = basename;
    file["it0"] = Json::Value::Int64(it0);
    file["nx"] = nx;

    g.curr_it1 = it0 + nt;
    g.files.append(file);

    return filename;
}


string pipeline_object::add_file(const string &basename)
{
    if (this->state != RUNNING)
	_throw("add_plot() must be called from _advance()");
    if (!out_mp)
	_throw("internal error: no outdir_manager in pipeline_object::add_file()");
    if (this->out_mp->outdir.size() == 0)
	_throw("add_file() was called in pipeline with no outdir (probably need to check out_mp->outdir.size() before calling add_plot_group)");

    return this->out_mp->add_file(basename);
}


}  // namespace rf_pipelines
