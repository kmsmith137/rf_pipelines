// FIXME the big missing feature here is pythonizing 'class ring_buffer', and
// extending the pythonization of 'pipeline_object' and 'chunked_pipeline_object'.
//
// FIXME: currently we can define wi_streams and wi_transforms from python, but not more
// general pipeline_objects.
//
// FIXME (minor): wi_stream, wi_transform constructors should have optional 'name' arg,
// for consistency with C++ API.

#include <pyclops.hpp>
#include <rf_kernels.hpp>
#include "rf_pipelines_internals.hpp"

using namespace std;
using namespace pyclops;
using namespace rf_pipelines;


// -------------------------------------------------------------------------------------------------
//
// C++ typedefs, converters


// Used frequently in the wrappers below.
using in_arr2d = in_ncarray<float,2,1>;
using io_arr2d = io_ncarray<float,2,1>;


namespace pyclops {

    template<>
    struct converter<rf_kernels::axis_type> {
	// The rf_kernels::axis converter allows a flexible sytnax.
	static rf_kernels::axis_type from_python(const py_object &x, const char *where = nullptr)
	{
	    if (x.is_none())
		return rf_kernels::AXIS_NONE;

	    if (x.is_integer()) {
		ssize_t i = converter<ssize_t>::from_python(x);
		if (i == 0)
		    return rf_kernels::AXIS_FREQ;
		if (i == 1)
		    return rf_kernels::AXIS_TIME;
	    }

	    if (x.is_string()) {
		string s = converter<string>::from_python(x);
		const char *c = s.c_str();

		if (!strcasecmp(c,"none") || !strcasecmp(c,"axis_none"))
		    return rf_kernels::AXIS_NONE;
		if (!strcasecmp(c,"time") || !strcasecmp(c,"axis_time"))
		    return rf_kernels::AXIS_TIME;
		if (!strcasecmp(c,"freq") || !strcasecmp(c,"axis_freq") || !strcasecmp(c,"frequency"))
		    return rf_kernels::AXIS_FREQ;
	    }

	    throw runtime_error((where ? string(where) : string("rf_kernels")) + ": bad 'axis' parameter");
	}
    };

    // FIXME need systematic framework for converting STL containers.
    // For now, just doing it on an ad hoc basis!

    template<typename T>
    struct predicated_converter<vector<T>, typename enable_if<converts_from_python<T>::value>::type> {
	static vector<T> from_python(const py_object &x, const char *where = nullptr)
	{
	    // This check is important, otherwise x='string' could be interpreted as the
	    // list-of-strings [ 's', 't', 'r', 'i', 'n', 'g' ].
	    if (PyString_Check(x.ptr))
		throw runtime_error((where ? string(where) : string("")) + ": expected list or iterator, got single string");

	    PyObject *iter_p = PyObject_GetIter(x.ptr);
	    if (!iter_p)
		throw runtime_error((where ? string(where) : string("")) + ": expected list or iterator, got non-iterable object");

	    py_object iter = py_object::new_reference(iter_p);
	    vector<T> ret;
	    
	    for (;;) {
		PyObject *item_p = PyIter_Next(iter_p);
		if (!item_p)
		    break;  // PyErr_Occurred() will be checked below
		
		py_object item = py_object::new_reference(item_p);
		T s = converter<T>::from_python(item, where);
		ret.push_back(s);
	    }

	    if (PyErr_Occurred())
		throw pyerr_occurred();

	    return ret;
	}
    };


    template<typename T, typename U>
    struct predicated_converter<pair<T,U>, typename enable_if<(converts_from_python<T>::value && converts_from_python<U>::value)>::type> {
	static pair<T,U> from_python(const py_object &x, const char *where = nullptr)
	{
	    // FIXME cut-and-paste from vector converter
	    // FIXME this could use general cleanup!

	    // This check is important, otherwise x='string' could be interpreted as the
	    // list-of-strings [ 's', 't', 'r', 'i', 'n', 'g' ].
	    if (PyString_Check(x.ptr))
		throw runtime_error((where ? string(where) : string("")) + ": expected list or iterator, got single string");

	    PyObject *iter_p = PyObject_GetIter(x.ptr);
	    if (!iter_p)
		throw runtime_error((where ? string(where) : string("")) + ": expected list or iterator, got non-iterable object");

	    py_object iter = py_object::new_reference(iter_p);
	    PyObject *item0_p = PyIter_Next(iter_p);

	    if (!item0_p) {
		if (PyErr_Occurred())
		    throw pyerr_occurred();
		throw runtime_error((where ? string(where) : string("")) + ": expected length-2 iterable object, got length-0");
	    }

	    py_object item0 = py_object::new_reference(item0_p);
	    PyObject *item1_p = PyIter_Next(iter_p);

	    if (!item1_p) {
		if (PyErr_Occurred())
		    throw pyerr_occurred();
		throw runtime_error((where ? string(where) : string("")) + ": expected length-2 iterable object, got length-1");
	    }

	    py_object item1 = py_object::new_reference(item1_p);
	    PyObject *item2_p = PyIter_Next(iter_p);

	    if (item2_p) {
		Py_DECREF(item2_p);
		throw runtime_error((where ? string(where) : string("")) + ": expected length-2 iterable object, got length >= 3");
	    }

	    if (PyErr_Occurred())
		throw pyerr_occurred();

	    return make_pair<T,U> (converter<T>::from_python(item0,where), converter<U>::from_python(item1,where));
	}
    };


    // "Deep" converter for Json::Value.
    // FIXME put somewhere more general (pyclops/addons/jsoncpp.hpp?)
    // FIXME doesn't check for circular reference, so can lead to infinite recursion (fixable)
    // FIXME how to convert uint64_t Json::Value to python?

    template<> 
    struct converter<Json::Value> {
	// Use case: pipeline_object::jsonize() is called from C++, and underlying object is implemented in python.
	static Json::Value from_python(const py_object &x, const char *where=nullptr)
	{
	    if (x.is_none())
		return Json::Value(Json::nullValue);
	    if (x.ptr == Py_True)
		return Json::Value(true);
	    if (x.ptr == Py_False)
		return Json::Value(false);
	    if (x.is_integer())
		return Json::Value((Json::Int64) converter<ssize_t>::from_python(x, where));
	    if (x.is_floating_point())
		return Json::Value(converter<double>::from_python(x, where));
	    if (x.is_string())
		return Json::Value(converter<string>::from_python(x, where));
	    
	    if (x.is_list()) {
		py_list l = x;
		ssize_t n = l.size();

		Json::Value ret(Json::arrayValue);
		for (ssize_t i = 0; i < n; i++)
		    ret.append(converter<Json::Value>::from_python(l.get_item(i), where));

		return ret;
	    }

	    if (x.is_dict()) {
		Json::Value ret(Json::objectValue);

		for (const auto &p: py_dict(x)) {
		    string key = converter<string>::from_python(p.first);
		    Json::Value val = converter<Json::Value>::from_python(p.second);
		    ret[key] = val;
		}

		return ret;
	    }

	    throw runtime_error(string(where ? where : "pyclops") + ": couldn't convert object of type '" + x.type_name() + "' to Json::Value");
	}
	    
	// Example use case: pipeline_object::jsonize() is called from python, and transform is implemented in C++.
	static py_object to_python(const Json::Value &x, const char *where=nullptr)
	{
	    if (x.isNull())
		return py_object();  // Py_None
	    if (x.isBool())
		return converter<bool>::to_python(x.asBool());
	    if (x.isIntegral())
		return converter<ssize_t>::to_python(x.asLargestInt());
	    if (x.isDouble())
		return converter<double>::to_python(x.asDouble());
	    if (x.isString())
		return converter<string>::to_python(x.asString());

	    if (x.isArray()) {
		py_list ret;
		for (const auto &p: x)
		    ret.append(converter<Json::Value>::to_python(p));
		return ret;
	    }

	    if (x.isObject()) {
		py_dict ret;
		for (auto p = x.begin(); p != x.end(); p++) {
		    string key = p.key().asString();
		    py_object val = converter<Json::Value>::to_python(*p);
		    ret.set_item(key, val);
		}
		return ret;
	    }

	    // If this exception ever gets thrown, there's something I don't understand about jsoncpp!
	    throw runtime_error("pyclops: internal error when converting jsoncpp value to python (couldn't classify object)");
	}
    };
} // namespace pyclops


// For constructing docstrings throughout this source file.
static string doc_axis_none = "    None or 'None' for AXIS_NONE\n";
static string doc_axis_freq = "    0 or 'freq' for AXIS_FREQ\n";
static string doc_axis_time = "    1 or 'time' for AXIS_TIME\n";


// -------------------------------------------------------------------------------------------------
//
// wrap_j(): specialized wrapper, for wrapping a member function of the form
//
//   void C::f(Json::Value &j)
//
// For example, pipeline_object::_bind(), _start_pipeline(), _end_pipeline().
//
// FIXME is there a good way to make this more general?  Would need to generalize
// the from-python converter to include a function which is called after the
// wrapped function has returned.


template<class C>
static std::function<py_object(C*, py_tuple, py_dict)> wrap_j(void (C::*f)(Json::Value &))
{
    std::function<void(C *, py_dict)>
	ret = [f](C *self, py_dict j)
	{
	    // We "convert" the (Json::Value &) by making a deep copy, passing it to
	    // the C++ function where it can be modified, then copying it back.
	    // This approach is bulletproof but potentially slow.  It should be OK
	    // in practice, since we only use wrap_j() in one-time initializations
	    // (and where the Json::Value is small anyway).
	    
	    Json::Value cpp_json = converter<Json::Value>::from_python(j);

	    (self->*f)(cpp_json);

	    for (auto p = cpp_json.begin(); p != cpp_json.end(); p++) {
		string key = p.key().asString();
		py_object val = converter<Json::Value>::to_python(*p);
		j.set_item(key, val);
	    }
	};

    return wrap_method(ret, "j");
}


// -------------------------------------------------------------------------------------------------
//
// Upcalling helpers.
//
// FIXME generalize and move to pyclops.


// "_nodef": no default virtual (or rather, default virtual does nothing)
template<typename T, typename B, typename TT>
static inline void _upcall_nodef(const extension_type<T,B> &type, TT *self, const char *method_name)
{
    virtual_function<void> vf(type, self, method_name);
    if (vf.exists)
	vf.upcall();
}


// "_nodef_j": virtual takes a (Json::Value &) and returns void.
template<typename T, typename B, typename TT>
static inline void _upcall_nodef_j(const extension_type<T,B> &type, TT *self, const char *method_name, Json::Value &j)
{
    virtual_function<void> vf(type, self, method_name);

    if (!vf.exists)
	return;  // no need to call default virtual, since it does nothing.

    
    // We "convert" the (Json::Value &) by making a deep copy, passing it to
    // the python function where it can be modified, then copying it back.
    // This approach is bulletproof but potentially slow.  It should be OK in
    // practice, since we only use _upcall_nodef_j() in one-time initializations
    // (and where the Json::Value is small anyway).

    py_object py_json = converter<Json::Value>::to_python(j);
    vf.upcall(py_json);
    j = converter<Json::Value>::from_python(py_json);
}


// -------------------------------------------------------------------------------------------------
//
// All PyTypeObjects


static string doc_pipeline_object =
    ("pipeline_object: rf_pipelines base class.\n"
     "\n"
     "Current pythonization is partial: can call run(), bind(), allocate(), jsonize(), etc.\n"
     "but new pipeline_objects can't be defined from python.  This is because 'class ring_buffer'\n"
     "isn't pythonized yet (this is on my todo list!)\n"
     "\n"
     "To run a pipeline, just call pipeline.run()!\n");

static string doc_chunked_pipeline_object =
    ("chunked_pipeline_object: represents any pipeline_object which processes in fixed-size chunks.\n"
     "\n"
     "Current pythonization is partial: new chunked_pipeline_objects can't be defined from python.\n"
     "This is because 'class ring_buffer' isn't pythonized yet (this is on my todo list!)\n");

static string doc_wi_stream =
    ("wi_stream: represents any pipeline_object which outputs intensity/weights arrays in fixed-size chunks");

static string doc_wi_transform =
    ("wi_transform: represents any pipeline_object which processes intensity/weights arrays in fixed-size chunks");

static string doc_pipeline =
    ("pipeline: this ubiquitous container class is used to chain pipeline_objects together.\n"
     "\n"
     "Constructor syntax:\n"
     "   p = pipeline(object_list=[], name=\"\")");

static string doc_wi_sub_pipeline =
    ("wi_sub_pipeline: this container class is used to run a \"sub-pipeline\" at\n"
     "lower (freqency, time) resolution, then upsample and apply the resulting mask.\n"
     "\n"
     "Constructor syntax:\n"
     "    p = wi_sub_pipeline(sub_pipeline, w_cutoff=0.0, nt_chunk=0, nfreq_out=0, nds_out=0, Df=0, Dt=0)\n"
     "\n"
     " 	  nfreq_out = number of frequency channels after downsampling to sub-pipeline\n"
     "    nds_out = cumulative time downsampling (relative to input data) after downsampling to sub-pipeline\n"
     "    Df = frequency downsampling factor (between input pipeline and sub-pipeline)\n"
     "    Dt = time downsampling factor (between input pipeline and sub-pipeline)\n"
     "\n"
     "The initializer allows a flexible syntax where some fields can be specified (i.e. nonzero)\n"
     "and others unspecified (i.e. zero).  For example:\n"
     "\n"
     "    - If 'nfreq_out' is specified and 'Df' is not, then Df will be set to (nfreq_in / nfreq_out).\n"
     "    - If 'nfreq_out' is unspecified and 'Dt' is specified, then nfreq_out will be set to (nfreq_in / Df).\n"
     "    - If 'nfreq_out' and 'Df' are both specified, then an exception will be raised unless nfreq_in = (nfreq_out * Df)\n"
     "    - If neither 'nfreq_out' nor 'Df' are specified, then an exception will be raised.\n"
     "\n"
     "The parameter pair (nds_out, Dt) behaves similarly.");


static extension_type<pipeline_object>
pipeline_object_type("pipeline_object", doc_pipeline_object);

static extension_type<chunked_pipeline_object, pipeline_object>
chunked_pipeline_object_type("chunked_pipeline_object", doc_chunked_pipeline_object, pipeline_object_type);

static extension_type<wi_stream, chunked_pipeline_object>
wi_stream_type("wi_stream", doc_wi_stream, chunked_pipeline_object_type);

static extension_type<wi_transform, chunked_pipeline_object>
wi_transform_type("wi_transform", doc_wi_transform, chunked_pipeline_object_type);

static extension_type<pipeline, pipeline_object>
pipeline_type("pipeline", doc_pipeline, pipeline_object_type);

static extension_type<wi_sub_pipeline, pipeline>
wi_sub_pipeline_type("wi_sub_pipeline", doc_wi_sub_pipeline, pipeline_type);


namespace pyclops {
    template<> struct xconverter<pipeline_object>  { static constexpr auto *type = &pipeline_object_type; };    
    template<> struct xconverter<chunked_pipeline_object>  { static constexpr auto *type = &chunked_pipeline_object_type; };
    
    template<> struct xconverter<wi_stream>        { static constexpr auto *type = &wi_stream_type; };
    template<> struct xconverter<wi_transform>     { static constexpr auto *type = &wi_transform_type; };
    template<> struct xconverter<pipeline>         { static constexpr auto *type = &pipeline_type; };
    template<> struct xconverter<wi_sub_pipeline>  { static constexpr auto *type = &wi_sub_pipeline_type; };
}


// -------------------------------------------------------------------------------------------------
//
// wrap_pipeline_object()
//
// Current pythonization is partial: can call run(), bind(), allocate(), jsonize(), etc.
// but new pipeline_objects can't be defined from python.  This is because 'class ring_buffer'
// isn't pythonized yet (this is on my todo list!)
//
// Note: the "upcalling" class
//   struct py_pipeline_object : pipeline_object { ... };
//
// is not needed, since subclassing from python isn't implemented yet.


// Helper function, passed as 'callback' in pipeline_object::run(), so that the pipeline
// will periodically check for control-C, and throw the appropriate exception.
static void check_signals(ssize_t pos_lo, ssize_t pos_hi)
{
    if (PyErr_Occurred() || PyErr_CheckSignals())
	throw pyerr_occurred();
}


// Used to wrap pipeline_object::bind() and pipeline_object::run().
static run_params make_run_params(const py_object &outdir, bool clobber, ssize_t img_nzoom, ssize_t img_nds, ssize_t img_nx, int verbosity, bool debug, const py_object &extra_attrs)
{
    run_params ret;

    // Allow outdir=None (equivalent to outdir="")
    if (outdir.is_none())
	ret.outdir = "";
    else if (outdir.is_string())
	ret.outdir = converter<string>::from_python(outdir, "outdir");
    else
	throw runtime_error("expected 'outdir' to be a string or None");

    // Allow extra_attrs=None (equivalent to extra_attrs={})
    if (extra_attrs.is_none())
	ret.extra_attrs = Json::Value(Json::objectValue);
    else
	ret.extra_attrs = converter<Json::Value>::from_python(extra_attrs, "extra_attrs");

    ret.clobber = clobber;    
    ret.img_nzoom = img_nzoom;
    ret.img_nds = img_nds;
    ret.img_nx = img_nx;
    ret.verbosity = verbosity;
    ret.debug = debug;

    return ret;
}


static void wrap_pipeline_object(extension_module &m)
{
    std::function<pipeline_object* ()>
	_init = []() -> pipeline_object* { throw runtime_error("rf_pipelines.pipeline_object cannot be directly subclassed from python (yet)"); };

    std::function<string& (pipeline_object *)>
	_name = [](pipeline_object *self) -> string& { return self->name; };

    std::function<string& (pipeline_object *)>
	_class_name = [](pipeline_object *self) -> string& { return self->class_name; };

    std::function<void (pipeline_object *, const py_object &, bool, ssize_t, ssize_t, ssize_t, int, bool, const py_object &)>
	_bind = [](pipeline_object *self, const py_object &outdir, bool clobber, ssize_t img_nzoom, ssize_t img_nds, ssize_t img_nx, int verbosity, bool debug, const py_object &extra_attrs)
	{
	    run_params p = make_run_params(outdir, clobber, img_nzoom, img_nds, img_nx, verbosity, debug, extra_attrs);
	    return self->bind(p);
	};

    std::function<Json::Value (pipeline_object *, const py_object &, bool, ssize_t, ssize_t, ssize_t, int, bool, const py_object &)>
	_run = [](pipeline_object *self, const py_object &outdir, bool clobber, ssize_t img_nzoom, ssize_t img_nds, ssize_t img_nx, int verbosity, bool debug, const py_object &extra_attrs)
	{
	    // FIXME for completeness, should allow python caller to specify a callback function.
	    // (This should be called via a C++ wrapper which also calls check_signals().)

	    run_params p = make_run_params(outdir, clobber, img_nzoom, img_nds, img_nx, verbosity, debug, extra_attrs);
	    return self->run(p, check_signals);
	};

    // __str__()
    // FIXME doesn't actually work yet!
    // FIXME define __repr__() which returns stringified json?
    std::function<string (pipeline_object *)>
	_str = [](pipeline_object *self) { return self->name; };

    // Wrapped as staticmethod pipeline_object.register_json_deserializer().
    std::function<void (const string &, const py_object &)>
	_register_json_deserializer = [](const string &class_name, const py_object &f)
	{
	    if (!f.is_callable())
		throw runtime_error("rf_pipelines.pipeline_object.register_json_deserializer(): Argument 'f' must be callable");

	    // FIXME define generic std::function from-python converter?
	    std::function<shared_ptr<pipeline_object> (const Json::Value &)>
	    _constructor = [f](const Json::Value &json_data)
	    {
		py_object py_json = converter<Json::Value>::to_python(json_data);
		py_tuple py_args = py_tuple::make(py_json);
		py_object py_ret = f.call(py_args);
		return converter<shared_ptr<pipeline_object>>::from_python(py_ret);
	    };

	    pipeline_object::register_json_deserializer(class_name, _constructor);
	};

    // doc_rp1, doc_rp2 are building blocks for doc_bind, doc_run, which both take a run_params.
    string doc_rp1 = ("outdir='.', clobber=True, img_nzoom=4, img_nds=16, img_nx=256, verbosity=2, debug=False, extra_attrs=None");

    string doc_rp2 = ("'outdir' is the rf_pipelines output directory, where the rf_pipelines json file will\n"
		      "be written, in addition to other transform-specific output files such as plots.\n"
		      "\n"
		      "If 'outdir' is an empty string, then the json file will not be written, and\n"
		      "any transform which tries to write an output file (such as a plotter_transform)\n"
		      "will throw an exception.\n"
		      "\n"
		      "If 'clobber' is false, then an exception will be thrown if the pipeline tries to\n"
		      "overwrite an old rf_pipelines.json file.\n"
		      "\n"
		      "Plot-related params:\n"
		      "   img_nzoom = number of zoom levels (FIXME currently hardcoded, should switch to adaptive)\n"
		      "   img_nds = time downsampling factor of plots at lowest zoom level\n"
		      "   img_nx = number of x-pixels (i.e. time axis) in each plot tile\n"
		      "\n"
		      "The meaning of the 'verbosity' argument is:\n"
		      "    0 = no output\n"
		      "    1 = high-level summary output (names of transforms, number of samples processed etc.)\n"
		      "    2 = show all output files\n"
		      "    3+ = debug trace through pipeline (larger value means that debug messages are printed to higher depth)\n"
		      "\n"
		      "If 'debug' is true, some extra debug tests are implemented.  This slows down\n"
		      "pipeline processing, so should only be specified for debugging/testing.\n"
		      "\n"
		      "If specified, 'extra_attrs' should be a dict containing extra json attributes for the pipeline run.\n"
		      "These attributes will be passed to _bind() and _start_pipeline(), and also end up in the pipeline json output.\n");

    string doc_bind = ("bind(" + doc_rp1 + ")\n"
		       "\n"
		       "Does global pipeline initializations, including ring buffer sizes.  Called automatically by run().\n"
		       "\n" + doc_rp2);

    string doc_run = ("run(" + doc_rp1 + ") -> json\n"
		      "\n"
		      "High-level API: to run a pipeline, just call run().\n"
		      "\n" + doc_rp2);

    string doc_cs = ("get_preferred_chunk_size() -> integer\n"
		     "\n"
		     "This is only called on the first pipeline_object in the pipeline\n"
		     "(subsequent chunk sizes are determined automatically).  By default,\n"
		     "this function returns 0, which results in an exception\n"
		     "   \"...: this pipeline_object cannot appear first in pipeline\".\n"
		     "\n"
		     "Stream-type pipeline_objects which can appear first in a pipeline should override\n"
		     "get_preferred_chunk_size() to return a nonzero value.");

    string doc_fj = ("from_json(json_data) -> pipeline_object\n"
		     "\n"
		     "This staticmethod is the counterpart of the non-static method jsonize().\n"
		     "It \"de-serializes\" json data, and returns a pipeline_object.");

    string doc_rjc = ("register_json_deserializer(class_name, f)\n"
		      "\n"
		      "The 'f' argument should be a function f(x) whose single argument 'x' is a json-serialized\n"
		      "object (e.g. output of the jsonize() method), and returns a pipeline_object instance.\n"
		      "\n"
		      "By convention, f() is usually chosen to be a staticmethod 'from_json' of the pipeline_object\n"
		      "ssubclass, e.g.:\n"
		      "\n"
		      "    class X(rf_pipelines.pipeline_object):\n"
		      "            ...\n"
		      "        @staticmethod\n"
		      "        def from_json(json_data):\n"
		      "            ...\n"
		      "            return X(...)\n"
		      "\n"
		      "    rf_pipeline.pipeline_object.register_json_converter('X', X.from_json)");

    string doc_add_pg = ("add_plot_group(name, nt_per_pix, ny) -> group_id [integer]\n"
			 "\n"
			 "This helper function is used by pipeline_objects which make output plots.\n"
			 "The output plots are divided into one or more \"plot groups\".\n"
			 "For example, the bonsai_dedisperser can write one plot group per internally defined tree.\n"
			 "The 'nt_per_pix' arg is the number of pipeline time samples per x-pixel in the plot.\n"
			 "The 'ny' arg is the number of y-pixels (assumed to be the same for all plots in the group).\n"
			 "\n"
			 "The return value is an integer group_id, which will be passed to pipeline_object.add_plot()\n"
			 "whenever a new plot is written.\n"
			 "\n"
			 "Note: add_plot_group() should be called in _start_pipeline(), not in the pipeline_object\n"
			 "constructor or in _bind().  This is because plot_groups are \"reset\" between pipeline runs.");

    string doc_add_plot = ("add_plot(basename, it0, nt, nx, ny, group_id=0) -> fullpath [string]n"
			   "\n"
			   "This helper function is used by pipeline_objects which make output plots.\n"
			   "It should be called just before writing each plot.\n"
			   "\n"
			   "The range of time samples in the plot is [it0:it0+nt), where (it0,nt) are \"upsampled\" time sample indices\n"
			   "at the \"native\" pipeline time resolution.  The pixel dimensions of the plot are (nx,ny).\n"
			   "\n"
			   "Note that nx is redundant, since nt should always equal (nt_per_pix * nds * nx), where nt_per_pix\n"
			   "was specified in add_plot_group() and nds is the current time downsampling factor in the pipeline.\n"
			   "Similarly, ny is redundant since it must match the value specified in add_plot_group(). Currently,\n"
			   "we require the caller to specify these values so that they can be used for error checking.\n"
			   "\n"
			   "The return value is the full pathname ('basename' with the pipeline output_dir prepended)");

    string doc_add_file = ("add_file(basename) -> fullpath [string]\n"
			   "\n"
			   "This helper function is used by pipeline_objects which write output files which are not plots.\n"
			   "(For plots, use the helper functions add_plot_group(), add_plot().)\n"
			   "\n"
			   "It should be called just before writing each file, to check for filename collisions.\n"
			   "The return value is the full pathname ('basename' with the pipeline output_dir prepended)");
			   
    auto _add_pg = wrap_method(&pipeline_object::add_plot_group, "name", "nt_per_pix", "ny");
    auto _add_plot = wrap_method(&pipeline_object::add_plot, "basename", "it0", "nt", "nx", "ny", kwarg("group_id",0));
    auto _add_file = wrap_method(&pipeline_object::add_file, "basename");
    
    pipeline_object_type.add_constructor(wrap_constructor(_init));
    pipeline_object_type.add_property("name", "Name of pipeline_object", _name);
    pipeline_object_type.add_property("class_name", "Name of pipeline_object subclass", _class_name);

    pipeline_object_type.add_method("run", doc_run, wrap_method(_run, kwarg("outdir",py_object()), kwarg("clobber",true), kwarg("img_nzoom",4), 
								kwarg("img_nds",16), kwarg("img_nx",256), kwarg("verbosity",2), kwarg("debug",false),
								kwarg("extra_attrs",py_object())));

    pipeline_object_type.add_method("bind", doc_bind, wrap_method(_bind, kwarg("outdir",py_object()), kwarg("clobber",true), kwarg("img_nzoom",4), 
								  kwarg("img_nds",16), kwarg("img_nx",256), kwarg("verbosity",2), kwarg("debug",false),
								  kwarg("extra_attrs",py_object())));
							 
    pipeline_object_type.add_method("allocate", "Allocates all pipeline buffers", wrap_method(&pipeline_object::allocate));
    pipeline_object_type.add_method("deallocate", "Deallocates all pipeline buffers", wrap_method(&pipeline_object::deallocate));
    pipeline_object_type.add_method("reset", "Resets pipeline after run().", wrap_method(&pipeline_object::reset));
    pipeline_object_type.add_method("unbind", "\"Unbinds\" pipeline, allowing its components to be used elsewhere", wrap_method(&pipeline_object::unbind));
    pipeline_object_type.add_method("get_info", "Can be called any time after bind(), returns nested json object", wrap_method(&pipeline_object::get_info));

    pipeline_object_type.add_method("get_preferred_chunk_size", doc_cs, wrap_method(&pipeline_object::get_preferred_chunk_size));

    pipeline_object_type.add_method("jsonize", "jsonize(): returns json serialization of pipeline_object", wrap_method(&pipeline_object::jsonize));
    pipeline_object_type.add_staticmethod("from_json", doc_fj, wrap_func(&pipeline_object::from_json, "json_data"));
    pipeline_object_type.add_staticmethod("register_json_deserializer", doc_rjc, wrap_func(_register_json_deserializer, "class_name", "f"));

    pipeline_object_type.add_method("add_plot_group", doc_add_pg, _add_pg);
    pipeline_object_type.add_method("add_plot", doc_add_plot, _add_plot);
    pipeline_object_type.add_method("add_file", doc_add_file, _add_file);

    // FIXME _bind() unwrapped for now, since ring_buffer class and friends are not python-wrapped
    pipeline_object_type.add_method("_allocate", "_allocate(): optional", wrap_method(&pipeline_object::_allocate));
    pipeline_object_type.add_method("_deallocate", "_deallocate(): optional", wrap_method(&pipeline_object::_deallocate));
    pipeline_object_type.add_method("_start_pipeline", "_start_pipeline(): optional", wrap_j(&pipeline_object::_start_pipeline));
    pipeline_object_type.add_method("_end_pipeline", "_end_pipeline(): optional", wrap_j(&pipeline_object::_end_pipeline));
    pipeline_object_type.add_method("_unbind", "_unbind(): optional", wrap_method(&pipeline_object::_unbind));
    pipeline_object_type.add_method("_reset", "_reset(): optional", wrap_method(&pipeline_object::_reset));
    pipeline_object_type.add_method("_get_info", "_get_info(): optional", wrap_j(&pipeline_object::_get_info));

    // FIXME doesn't really work -- more complicated than I thought!
    pipeline_object_type.add_method("__str__", "", wrap_method(_str));

    m.add_type(pipeline_object_type);
}   


// -------------------------------------------------------------------------------------------------
//
// wrap_chunked_pipeline_object
//
// Current pythonization is partial: new chunked_pipeline_objects can't be defined from python.
// This is because 'class ring_buffer' isn't pythonized yet (this is on my todo list!)


// Note: the "upcalling" class
//   struct py_chunked_pipeline_object : chunked_pipeline_object { ... };
//
// is not needed, since subclassing from python isn't implemented yet.


static void wrap_chunked_pipeline_object(extension_module &m)
{
    std::function<chunked_pipeline_object* ()>
	_init = []() -> chunked_pipeline_object* { throw runtime_error("rf_pipelines.chunked_pipeline_object cannot be directly subclassed from python (yet)"); };

    std::function<ssize_t& (chunked_pipeline_object *)>
	_nt_chunk = [](chunked_pipeline_object *self) -> ssize_t& { return self->nt_chunk; };
    
    chunked_pipeline_object_type.add_constructor(wrap_constructor(_init));
    chunked_pipeline_object_type.add_property("nt_chunk", "Chunk size of chunked_pipeline_object", _nt_chunk);
    chunked_pipeline_object_type.add_method("finalize_nt_chunk", "Initializes nt_chunk to a reasonable default, if not yet initialized", wrap_method(&chunked_pipeline_object::finalize_nt_chunk));

    m.add_type(chunked_pipeline_object_type);    
}


// -------------------------------------------------------------------------------------------------
//
// wrap_wi_stream()


struct py_wi_stream : wi_stream {
    py_wi_stream(const string &class_name_) :
	wi_stream(class_name_)
    { }

    // _fill_chunk() is the only pure virtual.
    virtual bool _fill_chunk(float *intensity, ssize_t istride, float *weights, ssize_t wstride, ssize_t pos) override
    {
	pure_virtual_function<bool> vf(wi_stream_type, this, "_fill_chunk");

	// Conversion of C++ arguments to python arguments is nontrivial here...
	// FIXME: figure out with 100% certainty which numpy flags should be specified in py_array::from_pointer().
	npy_intp sf = npy_intp(sizeof(float));
	npy_intp py_shape[2] = { nfreq, nt_chunk };
	npy_intp py_istrides[2] = { istride * sf, sf };
	npy_intp py_wstrides[2] = { wstride * sf, sf };
	py_array py_intensity = py_array::from_pointer(2, py_shape, py_istrides, sizeof(float), intensity, NPY_FLOAT, NPY_ARRAY_WRITEABLE);
	py_array py_weights = py_array::from_pointer(2, py_shape, py_wstrides, sizeof(float), weights, NPY_FLOAT, NPY_ARRAY_WRITEABLE);

	bool ret = vf.upcall(py_intensity, py_weights, pos);

	// FIXME I'd like to "invalidate" these arrays if a reference was kept.
	if ((py_intensity.get_refcount() > 1) || (py_weights.get_refcount() > 1))
	    throw runtime_error(name + ": python _fill_chunk() method is not allowed to keep references to its 'intensity' or 'weights' array arguments.  You may get a segfault!");
	
	return ret;
    }

    // Non-pure virtuals follow.

    virtual Json::Value jsonize() const override
    {
	virtual_function<Json::Value> vf(wi_stream_type, this, "jsonize");
	return vf.exists ? vf.upcall() : wi_stream::jsonize();
    }

    virtual void _bind_stream(Json::Value &j) override    { _upcall_nodef_j(wi_stream_type, this, "_bind_stream", j); }
    virtual void _start_pipeline(Json::Value &j) override { _upcall_nodef_j(wi_stream_type, this, "_start_pipeline", j); }
    virtual void _end_pipeline(Json::Value &j) override   { _upcall_nodef_j(wi_stream_type, this, "_end_pipeline", j); }
    virtual void _get_info(Json::Value &j) override       { _upcall_nodef_j(wi_stream_type, this, "_get_info", j); }

    virtual void _allocate() override       { _upcall_nodef(wi_stream_type, this, "_allocate"); }
    virtual void _deallocate() override     { _upcall_nodef(wi_stream_type, this, "_deallocate"); }
    virtual void _unbind_stream() override  { _upcall_nodef(wi_stream_type, this, "_unbind_stream"); }
    virtual void _reset() override          { _upcall_nodef(wi_stream_type, this, "_reset"); }
};


static void wrap_wi_stream(extension_module &m)
{
    // constructor for python subclasses of wi_stream
    std::function<wi_stream* (const string &)>
	_init = [](const string &class_name) -> wi_stream*
	{
	    return new py_wi_stream(class_name);
	};

    // nfreq property
    std::function<ssize_t& (wi_stream *)>
	_nfreq = [](wi_stream *self) -> ssize_t& { return self->nfreq; };

    // python-callable wi_stream::_fill_chunk().
    std::function<bool(wi_stream *, const io_arr2d &, const io_arr2d &, ssize_t)>
	_fill_chunk = [](wi_stream *self, const io_arr2d &intensity, const io_arr2d &weights, ssize_t pos) -> bool
	{
	    int istride = xdiv(intensity.stride(0), sizeof(float));
	    int wstride = xdiv(weights.stride(0), sizeof(float));
	    return self->_fill_chunk(intensity.data, istride, weights.data, wstride, pos);
	};

    wi_stream_type.add_constructor(wrap_constructor(_init, "class_name"));
    wi_stream_type.add_property("nfreq", "Number of frequency channels", _nfreq);
    wi_stream_type.add_method("_bind_stream", "_bind_stream(j): optional", wrap_j(&wi_stream::_bind_stream));
    wi_stream_type.add_method("_unbind_stream", "_unbind_stream(): optional", wrap_method(&wi_stream::_unbind_stream));
    
    wi_stream_type.add_method("_fill_chunk", "_fill_chunk(intensity, weights, pos)", 
			      wrap_method(_fill_chunk, "intensity", "weights", "pos"));

    m.add_type(wi_stream_type);
}


// -------------------------------------------------------------------------------------------------
//
// wrap_wi_transform()


struct py_wi_transform : wi_transform {
    py_wi_transform(const string &class_name_) :
	wi_transform(class_name_)
    { }

    // _process_chunk() is the only pure virtual.
    virtual void _process_chunk(float *intensity, ssize_t istride, float *weights, ssize_t wstride, ssize_t pos) override
    {
	pure_virtual_function<void> vf(wi_transform_type, this, "_process_chunk");

	// Conversion of C++ arguments to python arguments is nontrivial here...
	// FIXME: figure out with 100% certainty which numpy flags should be specified in py_array::from_pointer().
	npy_intp sf = npy_intp(sizeof(float));
	npy_intp py_shape[2] = { nfreq, nt_chunk };
	npy_intp py_istrides[2] = { istride * sf, sf };
	npy_intp py_wstrides[2] = { wstride * sf, sf };
	py_array py_intensity = py_array::from_pointer(2, py_shape, py_istrides, sizeof(float), intensity, NPY_FLOAT, NPY_ARRAY_WRITEABLE);
	py_array py_weights = py_array::from_pointer(2, py_shape, py_wstrides, sizeof(float), weights, NPY_FLOAT, NPY_ARRAY_WRITEABLE);

	vf.upcall(py_intensity, py_weights, pos);

	// FIXME I'd like to "invalidate" these arrays if a reference was kept.
	if ((py_intensity.get_refcount() > 1) || (py_weights.get_refcount() > 1))
	    throw runtime_error(name + ": python _fill_chunk() method is not allowed to keep references to its 'intensity' or 'weights' array arguments.  You may get a segfault!");
    }

    // Non-pure virtuals follow.

    virtual Json::Value jsonize() const override
    {
	virtual_function<Json::Value> vf(wi_transform_type, this, "jsonize");
	return vf.exists ? vf.upcall() : wi_transform::jsonize();
    }

    virtual void _bind_transform(Json::Value &j) override { _upcall_nodef_j(wi_transform_type, this, "_bind_transform", j); }
    virtual void _start_pipeline(Json::Value &j) override { _upcall_nodef_j(wi_transform_type, this, "_start_pipeline", j); }
    virtual void _end_pipeline(Json::Value &j) override   { _upcall_nodef_j(wi_transform_type, this, "_end_pipeline", j); }
    virtual void _get_info(Json::Value &j) override       { _upcall_nodef_j(wi_transform_type, this, "_get_info", j); }

    virtual void _allocate() override { _upcall_nodef(wi_transform_type, this, "_allocate"); }
    virtual void _deallocate() override { _upcall_nodef(wi_transform_type, this, "_deallocate"); }
};


static void wrap_wi_transform(extension_module &m)
{
    // Properties
    std::function<ssize_t& (wi_transform *)> _nfreq = [](wi_transform *self) -> ssize_t& { return self->nfreq; };
    std::function<ssize_t& (wi_transform *)> _nds = [](wi_transform *self) -> ssize_t& { return self->nds; };

    // constructor for python subclasses of wi_transform
    std::function<wi_transform* (const string&)>
	_init = [](const string &class_name)
	{
	    return new py_wi_transform(class_name);
	};

    // python-callable wi_transform::_process_chunk()
    std::function<void(wi_transform *, const io_arr2d &, const io_arr2d &, ssize_t)>
	_process_chunk = [](wi_transform *self, const io_arr2d &intensity, const io_arr2d &weights, ssize_t pos)
	{
	    int istride = xdiv(intensity.stride(0), sizeof(float));
	    int wstride = xdiv(weights.stride(0), sizeof(float));
	    self->_process_chunk(intensity.data, istride, weights.data, wstride, pos);
	};

    wi_transform_type.add_constructor(wrap_constructor(_init, "class_name"));
    wi_transform_type.add_property("kernel_chunk_size", "Kernel chunk size (optional)", _nfreq);
    wi_transform_type.add_property("nfreq", "Number of frequency channels", _nfreq);
    wi_transform_type.add_property("nds", "Time downsampling factor", _nds);

    wi_transform_type.add_method("_bind_transform", "_bind_transform(): optional", wrap_j(&wi_transform::_bind_transform));
    wi_transform_type.add_method("_unbind_transform", "_unbind_transform(): optional", wrap_method(&wi_transform::_unbind_transform));

    wi_transform_type.add_method("_process_chunk", "_process_chunk(intensity, weights, pos)",
				 wrap_method(_process_chunk, "intensity", "weights", "pos"));

    m.add_type(wi_transform_type);
}


// -------------------------------------------------------------------------------------------------
//
// wrap_containers()


static void wrap_containers(extension_module &m)
{
    using chain_t = vector<shared_ptr<pipeline_object>>;
	
    std::function<pipeline* (const chain_t &, const string &)>
	_pipeline_init = [](const chain_t &stages, const string &name) { return new pipeline(stages,name); };

    std::function<int (const pipeline *)>
	_pipeline_size = [](const pipeline *p) { return p->size(); };

    std::string doc_add = ("add(pipeline_object)\n"
			   "\n"
			   "Appends a pipeline_object to the pipeline.");
    
    pipeline_type.add_constructor(wrap_constructor(_pipeline_init, kwarg("object_list",chain_t()), kwarg("name",string())));
    pipeline_type.add_method("add", doc_add, wrap_method(&pipeline::add, "object"));
    pipeline_type.add_property("size", "Number of pipeline_objects in pipeline container", _pipeline_size);

    std::function<wi_sub_pipeline* (const shared_ptr<pipeline_object> &ds_pipeline, double, ssize_t, ssize_t, ssize_t, ssize_t, ssize_t)>
	_ws_init = [](const shared_ptr<pipeline_object> &ds_pipeline, double w_cutoff, ssize_t nt_chunk, ssize_t nfreq_out, ssize_t nds_out, ssize_t Df, ssize_t Dt)
	{
	    wi_sub_pipeline::initializer ini_params;
	    ini_params.w_cutoff = w_cutoff;
	    ini_params.nt_chunk = nt_chunk;
	    ini_params.nfreq_out = nfreq_out;
	    ini_params.nds_out = nds_out;
	    ini_params.Df = Df;
	    ini_params.Dt = Dt;

	    return new wi_sub_pipeline(ds_pipeline, ini_params);
	};

    auto ws_init = wrap_constructor(_ws_init, "sub_pipeline", kwarg("w_cutoff",0.0), kwarg("nt_chunk",0), 
				    kwarg("nfreq_out",0), kwarg("nds_out",0), kwarg("Df",0), kwarg("Dt",0));

    wi_sub_pipeline_type.add_constructor(ws_init);

    m.add_type(pipeline_type);
    m.add_type(wi_sub_pipeline_type);
}


// -------------------------------------------------------------------------------------------------
//
// wrap_utility_classes()


static void wrap_utility_classes(extension_module &m)
{
    m.add_function("mask_expander",
		   "mask_expander(axis, prev_wname, width, threshold, alpha=0.0, nt_chunk=0) -> pipeline_object\n"
		   "\n"
		   "Experimental: expands the RFI mask, in a way which is intended to \"fill gaps\".\n"
		   "\n"
		   "It is assuemd that the caller has saved the weights at a previous point in the pipeline\n"
		   "(using pipeline_fork, see the 'pipeline_fork' docstring).  We use the term \"prev_mask\" to\n"
		   "mean the RFI mask at this previous point in the pipeline, and \"delta_mask\" to mean the set\n"
		   "of pixels which are currently masked in the pipeline, but were not masked in the prev_mask.\n"
		   "\n"
		   "By default, the mask_expander actually expands the delta-mask, but this behavior\n"
		   "can be modified (see the 'alpha' parameter below).\n"
		   "\n"
		   "The expansion is done by computing exponential moving averages of the delta-mask in\n"
		   "both directions, and masking pixels when both averages are above a threshold.  This\n"
		   "will be written up in more detail later!\n"
		   "\n"
		   "Constructor arguments\n"
		   "---------------------\n"
		   "\n"
		   "'axis': currently, only AXIS_FREQ is implemented.  AXIS_TIME is coming soon!\n"
		   "\n"
		   "'prev_wname': pipeline bufname of the saved weights (a string).  Note that in order\n"
		   "  to save the weights at a previous point in the pipeline, you can use a pipeline_fork\n"
		   "  whose input_bufname parameter is \"WEIGHTS\" and whose output_bufname is a string which\n"
		   "  uniquely identifies the saved weights (e.g. \"WEIGHTS_SAVE1\").  The 'prev_wname' argument\n"
		   "  to the mask_expander should be the same as the 'output_bufname' argument of the\n"
		   "  pipeline_fork.\n"
		   "\n"
		   "'width': the decay width of the exponential moving average.  In the AXIS_FREQ case,\n"
		   "  this is expressed as a fraction of the frequency band, i.e. width=0.1 means that\n"
		   "  the characteristic width of the mask_expander is 10% of the full frequency band.\n"
		   "\n"
		   "'threshold': value between 0 and 1 which determines how aggressive the mask_expander is.\n"
		   "  Low values correspond to more masking.  The numerical value can be roughly interpreted\n"
		   "  as the fraction of data which must be delta-masked before mask expansion will\n"
		   "  occur.  For example, if threshold=0.1, then mask expansion will occur in regions of\n"
		   "  the data where ~10% or more of the pixels are delta-masked.\n"
		   "\n"
		   "'alpha': to explain this parameter, we first note that delta-masked pixels are\n"
		   "  \"sources\" for the mask_expander, and unmasked pixels are \"sinks\".  That is,\n"
		   "  mask expansion occurs in regions where the number of delta-masked pixels\n"
		   "  relative to the number of unmasked pixels is above a threshold.\n"
		   "\n"
		   "  The alpha paramaeter determines how the prev_mask is handled by the mask_expander.\n"
		   "  By default (alpha=0), prev-masked pixels are \"neutral\", i.e. they are neither\n"
		   "  sources nor sinks for the mask_expander.\n"
		   "\n"
		   "  If 0 < alpha < 1, then prev-masked pixels are sinks for the mask_expander, i.e.\n"
		   "  they reduce the amount of mask expansion, and the amount of reduction is proportional\n"
		   "  to alpha.  If alpha=1, then prev_masked pixels are equivalent to unmasked pixels.\n"
		   "\n"
		   "  If -1 < alpha < 0, then prev-masked pixels are sources for the mask_expander, i.e.\n"
		   "  they increase the amount of mask expansion, and the amount of reduction is proportional\n"
		   "  to (-alpha).  If alpha=-1, then prev_masked pixels are equivalent to delta_masked pixels.\n",
		   wrap_func(make_mask_expander, "axis", "prev_wname", "width", "threshold", kwarg("alpha",0.0), kwarg("nt_chunk",0)));
		   
    m.add_function("pipeline_fork",
		   "pipeline_fork(bufnames) -> pipeline_object\n"
		   "\n"
		   "Creates one or more new pipeline ring_buffers, by copying existing ring_buffers.\n"
		   "The 'bufnames' argument should be a list of (input_bufname, output_bufname) pairs.\n"
		   "Frequently, the input_bufname will be one of the built-in names \"INTENSITY\" or \"WEIGHTS\".\n",
		   wrap_func(make_pipeline_fork, "bufnames"));
}


// -------------------------------------------------------------------------------------------------
//
// wrap_chime_streams()


static void wrap_chime_streams(extension_module &m)
{
    // For building docstrings.
    string doc_fs = ("The 'nt_chunk' arg is the chunk size used internally when moving data from hdf5 file\n"
		     "into the rf_pipelines buffer.  If unspecified or zero, it will default to a reasonable value.\n"
		     "\n"
		     "If 'noise_source_align' is nonzero, then it should be equal to the DETRENDER chunk size (not the chime_file_stream nt_chunk).\n"
		     "In this case, the stream will align the noise source edges with the detrender chunks, by discarding initial data if necessary.\n"
		     "\n"
		     "Note: functions beginning 'chime_stream..' are HDF5 streams, whereas functions beginning 'chime_frb_stream...' are msgpack.\n"
		     "For example, chime_stream_from_filename() and chime_frb_stream_from_filename() create streams from a single HDF5 or msgpack\n"
		     "file, respectively.");

    m.add_function("chime_stream_from_acqdir",
		   "chime_stream_from_acqdir(dirname, nt_chunk=0, noise_source_align=0, nfiles=0)\n\n"
		   "Makes a CHIME data stream from a directory containing HDF5 files.\n"
		   "The directory is scanned for filenames of the form NNNNNNNN.h5, where N=[0,9].\n"
		   "The 'nfiles' optional argument can be used to limit the acquisition to the first N files.\n\n" + doc_fs,
		   wrap_func(make_chime_stream_from_acqdir, "dirname", kwarg("nt_chunk",0), kwarg("noise_source_align",0), kwarg("nfiles",0)));

    m.add_function("chime_stream_from_filename",
		   "chime_stream_from_filename(filename, nt_chunk=0, noise_source_align=0)\n\n"
		   "Makes a CHIME data stream from a single HDF5 file.\n\n" + doc_fs,
		   wrap_func(make_chime_stream_from_filename, "filename", kwarg("nt_chunk",0), kwarg("noise_source_align",0)));

    m.add_function("chime_stream_from_filename_list",
		   "chime_stream_from_filename_list(filename_list, nt_chunk=0, noise_source_align=0)\n\n"
		   "Makes a CHIME data stream from a python list of HDF5 filenames.\n\n" + doc_fs,
		   wrap_func(make_chime_stream_from_filename_list, "filename_list", kwarg("nt_chunk",0), kwarg("noise_source_align",0)));

    m.add_function("chime_frb_stream_from_glob",
		   "chime_frb_stream_from_glob(glob_pattern, nt_chunk=0, noise_source_align=0)\n\n"
		   "Makes a CHIME data stream from a glob pattern matching msgpack files.\n\n" + doc_fs,
		   wrap_func(make_chime_frb_stream_from_glob, "glob_pattern", kwarg("nt_chunk",0), kwarg("noise_source_align",0)));

    m.add_function("chime_frb_stream_from_filename",
		   "chime_frb_stream_from_filename(filename, nt_chunk=0, noise_source_align=0)\n\n"
		   "Makes a CHIME data stream from a single msgpack file.\n\n" + doc_fs,
		   wrap_func(make_chime_frb_stream_from_filename, "filename", kwarg("nt_chunk",0), kwarg("noise_source_align",0)));

    m.add_function("chime_frb_stream_from_filename_list",
		   "chime_frb_stream_from_filename_list(filename, nt_chunk=0, noise_source_align=0)\n"
		   "Makes a CHIME data stream from a python list of msgpack filenames.\n\n" + doc_fs,
		   wrap_func(make_chime_frb_stream_from_filename_list, "filename_list", kwarg("nt_chunk",0), kwarg("noise_source_align",0)));

    // There are two versions of chime_network_stream(), one which takes an integer udp_port,
    // and another which takes a shared_ptr<ch_frb_io::intensity_network_stream>.  We only python-wrap
    // the former, since 'class ch_frb_io::intensity_network_stream' is not python-wrapped.

    using ns_t = shared_ptr<wi_stream> (*) (int,int,float);

    m.add_function("chime_network_stream",
		   "chime_network_stream(udp_port=0, beam_id=0,prescale=1.0)\n"
		   "\nIf the udp_port is zero, then the default chimefrb port will be used.",
		   wrap_func((ns_t) make_chime_network_stream, kwarg("udp_port",0), kwarg("beam_id",0), kwarg("prescale",1.0)));

    m.add_function("chime_dummy_network_stream",
		   "chime_dummy_network_stream(nt_tot, nupfreq=16, nt_per_packet=16, fpga_counts_per_sample=384, pool_gb=1.0)\n"
		   "\n"
		   "\"Dummy\" CHIME network stream, intended for timing.\n"
		   "Returns a stream which decodes a preallocated ch_frb_io::assembled_chunk pool as the pipeline progresses,\n"
		   "but does not actually receive packets over the network.  This allows the CPU cost of assembled_chunk\n"
		   "decoding to be included in pipeline timings.  Default parameter values are appropriate for full CHIME.",
		   wrap_func(make_dummy_chime_network_stream, "nt_tot", kwarg("nupfreq",16), kwarg("nt_per_packet",16), kwarg("fpga_counts_per_sample",384), kwarg("pool_gb",1.0)));

    m.add_function("chime_16k_spike_mask",
		   "chime_16k_spike_mask(nt_chunk=0)\n"
		   "Experimental: removes \"spikes\" from 16K data.",		   
		   wrap_func(make_chime_16k_spike_mask, kwarg("nt_chunk",0)));

    m.add_function("chime_16k_derippler",
		   "chime_16k_derippler(fudge_factor=1.0, nt_chunk=0)\n"
		   "Experimental: removes \"ripples\" from 16K data.",
		   wrap_func(make_chime_16k_derippler, kwarg("fudge_factor",1.0), kwarg("nt_chunk",0)));

    m.add_function("chime_16k_stripe_analyzer",
		   "chime_16k_stripe_analyzer(Dt1=16, Df2=16, Dt2=16)\n"
		   "Experimental: analyzes 16k-ripples and writes result to HDF5 file for follow-up analysis.",
		   wrap_func(make_chime_16k_stripe_analyzer, kwarg("Dt1",16), kwarg("Df2",16), kwarg("Dt2",16)));
    
    m.add_function("spectrum_analyzer",
		   "spectrum_analyzer(Dt1=16, Dt2=16)\n"
		   "Experimental: writes HDF5 file containing the intensity spectrum.",
		   wrap_func(make_spectrum_analyzer, kwarg("Dt1",16), kwarg("Dt2",16)));
}


// -------------------------------------------------------------------------------------------------
//
// wrap_chime_ostreams()


static void wrap_chime_ostreams(extension_module &m)
{
    string doc_fw = ("chime_file_writer(filename, clobber=False, bitshuffle=2, nt_chunk=0)\n"
		     "\n"
		     "This is a pseudo-transform which doesn't actually modify the data, it just writes it to a file in\n"
		     "CHIME hdf5 format.  (For now, the entire stream is written to a single file, I'll generalize later\n"
		     "to break the stream into multiple files.\n"
		     "\n"
		     "If 'clobber' is false, and the target file already exists, an exception will be thrown rather than clobbering the old file.\n"
		     "If 'nt_chunk' is set to zero, a default chunk size will be chosen.\n"
		     "\n"
		     "The meaning of the 'bitshuffle' arg is:\n"
		     "   0 = no compression\n"
		     "   1 = try to compress, but if plugin fails then just write uncompressed data instead\n"
		     "   2 = try to compress, but if plugin fails then print a warning and write uncompressed data instead\n"
		     "   3 = compression mandatory\n"
		     "\n"
		     "Note: a quick way to inspect a CHIME hdf5 file is using the 'ch-show-intensity-file' and 'ch-plot-intensity-file'\n"
		     "programs, in the ch_frb_io github repo.");

    string doc_cp = ("chime_packetizer(dstname, nfreq_coarse_per_packet, nt_per_chunk, nt_per_packet, wt_cutoff, target_gpbs, beam_id=0)\n"
		     "\n"
		     "Converts a stream to UDP packets in \"CHIME L0_L1\" format, and sends them over the network.\n"
		     "This interface is less general than the low-level interface in ch_frb_io: only one beam can\n"
		     "be sent, and not all boolean options are supported.\n"
		     "\n"
		     "Some artificial restrictions: the stream 'nfreq' value must be a multiple of 1024, and\n"
		     "the stream 'dt_sample' value must be an integer multiple of 2.56e-6 seconds.  This is because\n"
		     "the packet protocol doesn't include a count of total frequency channels, or the fpga clock\n"
		     "rate, so these parameters are frozen to the CHIME instrumental values.\n"
		     "\n"
		     "The 'dstname' argument is a string of the form HOSTNAME:PORT.  For example 'localhost:13178' or\n"
		     "'chimer.physics.ubc.ca:13178'.  If the port is omitted then the default chimefrb port is used.\n"
		     "(Be careful sending packets over the internet since the bandwidth can be very high!)\n"
		     "\n"
		     "The 'wt_cutoff' argument is used to convert the rf_pipelines 'weights' array to a boolean mask.\n"
		     "This conversion is necessary because the CHIME L0_L1 packet format doesn't support a floating-point\n"
		     "weight array.  Samples with weight below the cutoff will be masked.\n"
		     "\n"
		     "If the 'target_gbps' argument is nonzero, then output will be \"throttled\" to the target bandwidth, specified\n"
		     "in Gbps.  If target_gbps=0, then packets will be sent as quickly as possible.\n"
		     "\n"
		     "The nfreq_coarse_per_packet, nt_per_packet arguments define the amount of data sent per packet.\n"
		     "The nt_per_chunk arg just determines an internal chunk size and isn't very important (must be\n"
		     "a multiple of nt_per_packet; suggest a value like 512).\n");
    
    auto f_fw = wrap_func(make_chime_file_writer, "filename", kwarg("clobber",false), kwarg("bitshuffle",2), kwarg("nt_chunk",0));

    auto f_cp = wrap_func(make_chime_packetizer, "dstname", "nfreq_coarse_per_packet", "nt_per_chunk", 
			  "nt_per_packet", "wt_cutoff", "target_gpbs", kwarg("beam_id",0));

    m.add_function("chime_file_writer", doc_fw, f_fw);
    m.add_function("chime_packetizer", doc_cp, f_cp);

    string doc_aw = "chime_assembled_chunk_file_writer(filename, clobber=False)\n";
    auto f_aw = wrap_func(make_chime_assembled_chunk_file_writer, "filename", kwarg("clobber",false), kwarg("beams", std::vector<int>()));
}


// -------------------------------------------------------------------------------------------------
//
// wrap_misc_streams()


static void wrap_misc_streams(extension_module &m)
{
    string doc_gs = ("gaussian_noise_stream(nfreq, nt_tot, freq_lo_MHz, freq_hi_MHz, dt_sample, sample_rms=1.0, nt_chunk=1024, randomize_weights=false)\n"
		     "\n"
		     "A simple stream which simulates Gaussian random noise.\n"
		     "\n"
		     "Constructor arguments:\n"
		     "   nfreq               Number of frequency channels\n"
		     "   nt_tot              Total number of time samples written before stream ends.\n"
		     "   freq_lo_MHz         Lowest frequency in band (e.g. 400 for CHIME)\n"
		     "   freq_hi_MHz         Highest frequency in band (e.g. 800 for CHIME)\n"
		     "   dt_sample           Length of a time sample in seconds\n"
		     "   nt_chunk            Stream block size (if zero, will default to a reasonable value)\n"
		     "   randomize_weights   If true, weights will be uniform random numbers (if false, all weights will be 1.0)\n");
    
    auto f_gs = wrap_func(make_gaussian_noise_stream, "nfreq", "nt_tot", "freq_lo_MHz", "freq_hi_MHz", "dt_sample",
			  kwarg("sample_rms",1.0), kwarg("nt_chunk",1024), kwarg("randomize_weights",false));

    m.add_function("gaussian_noise_stream", doc_gs, f_gs);
}


// -------------------------------------------------------------------------------------------------
//
// wrap_detrenders()


static void wrap_detrenders(extension_module &m)
{
    m.add_function("polynomial_detrender",
		   "polynomial_detrender(nt_chunk, axis, polydeg, epsilon=1.0e-2)\n\n"
		   "Detrends along either the time or frequency axis, by subtracting a best-fit polynomial.\n"
		   "The detrending is independent in every \"row\" (where \"row\" means \"frequency channel\" in the\n"
		   "case of time-axis detrending, or \"time sample\" in the case of frequency-axis detrending).\n\n"
		   "If the fit is poorly conditioned then the entire row will be masked (by setting its weights to zero).\n"
		   "The threshold is controlled by the parameter 'epsilon'.  I think that 1.0e-2 is a reasonable default here,\n"
		   "but haven't experimented systematically.\n\n"
		   "The 'axis' argument should be one of the following:\n" + doc_axis_freq + doc_axis_time,
		   wrap_func(make_polynomial_detrender, "nt_chunk", "axis", "polydeg", kwarg("epsilon",1.0e-2)));

    m.add_function("spline_detrender",
		   "spline_detrender(nt_chunk, axis, nbins, epsilon=3.0e-4)\n\n"
		   "Experimental: spline_detrender.\n"
		   "I suspect this will work better than the polynomial_detrender, and it will definitely be faster!\n"
		   "Currently, the only allowed axis type is AXIS_FREQ, which can be specified with either axis='freq' or axis=0.\n",
		   wrap_func(make_spline_detrender, "nt_chunk", "axis", "nbins", kwarg("epsilon",3.0e-4)));
}


// -------------------------------------------------------------------------------------------------
//
// wrap_clippers()


static void wrap_clippers(extension_module &m)
{
    string doc_bm = ("badchannel_mask(mask_path, mask_ranges=None)\n"
		     "\n"
		     "The badchannel_mask transform sets bad freq channels of a weights array to 0.\n"
		     "\n"
		     "'mask_path' is the full path to a mask file that contains affected freq\n"
		     "intervals, written in rows with the following format: e.g., 420.02,423.03.\n"
		     "If mask_path is an empty string (or None), then no mask file will be read.\n"
		     "\n"
		     "'mask_ranges' is a list of (freq_lo, freq_hi) pairs, which define additional\n"
		     "frequency ranges to be masked.  If mask_ranges is an empty list (or None),\n"
		     "then no additional ranges will be masked.\n");

    string doc_ic = ("intensity_clipper(nt_chunk, axis, sigma, niter=1, iter_sigma=0, Df=1, Dt=1, two_pass=false)\n"
		     "\n"
		     "intensity_clipper: this \"clips\" an array by masking outlier intensities.\n"
		     "The masking is performed by setting elements of the weights array to zero.\n"
		     "\n"
		     "The 'sigma' argument is the threshold (in sigmas from the mean) for clipping.  Note\n"
		     "that the weights are used when calculating both the mean and rms intensity.\n"
		     "\n"
		     "The (Df,Dt) args are downsampling factors on the frequency/time axes.\n"
		     "If no downsampling is desired, set Df=Dt=1.\n"
		     "\n"
		     "The 'axis' argument has the following meaning:\n"
		     "   axis=0 or 'freq'   clip along frequency axis, with an outer loop over time samples\n"
		     "   axis=1 or 'time'   clip along time axis, with an outer loop over frequency samples\n"
		     "   axis=None          2-d clipper\n"
		     "\n"
		     "If niter > 1, then the mean/rms intensity will be computed using iterated clipping,\n"
		     "with threshold 'iter_sigma'.  If the 'iter_sigma' argument is zero, then it defaults\n"
		     "to 'sigma', but the two thresholds need not be the same.\n"
		     "\n"
		     "If the 'two_pass' flag is set, a more numerically stable but slightly slower algorithm will be used.");

    
    string doc_sdc = ("std_dev_clipper()\n"
		      "\n"
		      "std_dev_clipper: this \"clips\" an array by masking rows/columns whose standard deviation is an outlier.\n"
		      "\n"
		      "The 'axis' argument has the following meaning:\n"
		      "   axis=0 or 'freq'   clip time samples whose variance in frequency is high\n"
		      "   axis=1 or 'time'   clip frequency channels whose variance in time is high\n"
		      "\n"
		      "The (Df,Dt) args are downsampling factors on the frequency/time axes.\n"
		      "If no downsampling is desired, set Df=Dt=1.\n"
		      "\n"
		      "The 'sigma' argument is the threshold (in sigmas from the mean) for clipping.\n"
		      "\n"
		      "If the 'two_pass' flag is set, a more numerically stable but slightly slower algorithm will be used.");

    
    // python-callable make_badchannel_mask()
    std::function<shared_ptr<pipeline_object>(py_object, py_object)>
	_make_bm = [](py_object py_mask_path, py_object py_mask_ranges) -> shared_ptr<pipeline_object>
	{
	    string cpp_mask_path;
	    vector<pair<double,double>> cpp_mask_ranges;

	    if (!py_mask_path.is_none())
		cpp_mask_path = converter<string>::from_python(py_mask_path, "badchannel_mask: 'mask_path'");

	    if (!py_mask_ranges.is_none())
		cpp_mask_ranges = converter<vector<pair<double,double>>>::from_python(py_mask_ranges, "badchannel_mask: 'mask_ranges'");

	    return make_badchannel_mask(cpp_mask_path, cpp_mask_ranges);
	};


    auto f_bm = wrap_func(_make_bm, "mask_path", kwarg("mask_ranges",py_object()));
    
    auto f_ic = wrap_func(make_intensity_clipper, "nt_chunk", "axis", "sigma", kwarg("niter",1),
			  kwarg("iter_sigma",0.0), kwarg("Df",1), kwarg("Dt",1), kwarg("two_pass",false));

    auto f_sd = wrap_func(make_std_dev_clipper, "nt_chunk", "axis", "sigma", kwarg("Df",1), kwarg("Dt",1), kwarg("two_pass",false));
    

    m.add_function("badchannel_mask", doc_bm, f_bm);
    m.add_function("intensity_clipper", doc_ic, f_ic);
    m.add_function("std_dev_clipper", doc_sdc, f_sd);
}		   


static void wrap_mask_counters(extension_module &m)
{
    string doc_mc = ("mask_counter(nt_chunk, where)\n"
		     "\n"
                     "mask_counter: this counts how many intensity samples per chunk have been masked out by previous steps in the RFI chain.\n"
                     "The 'where' argument, which must be a unique string (unique within the pipeline), is used for reporting where in the pipeline the measurement is being made.\n"
		     "");

    auto f_mc = wrap_func(make_mask_counter, "nt_chunk", "where");
    m.add_function("mask_counter", doc_mc, f_mc);
}

// -------------------------------------------------------------------------------------------------
//
// wrap_kernels().  (I think these are only used in unit tests now.)
//
// FIXME it would be better to have an rf_kernels python module which exports these!
// FIXME minor loose end: no kernel wrapper for spline_detrender


static void apply_polynomial_detrender(io_arr2d intensity, io_arr2d weights, rf_kernels::axis_type axis, int polydeg, double epsilon=1.0e-2)
{
    if ((intensity.shape(0) != weights.shape(0)) || (intensity.shape(1) != weights.shape(1)))
	throw runtime_error("rf_pipelines.apply_polynomial_detrender: 'intensity' and 'weights' arrays do not have the same shape");

    // Additional argument checks will be performed in kernel constructor.
    rf_kernels::polynomial_detrender kernel(axis, polydeg);
    
    int nfreq = intensity.shape(0);
    int nt = intensity.shape(1);
    int istride = xdiv(intensity.stride(0), sizeof(float));
    int wstride = xdiv(weights.stride(0), sizeof(float));

    kernel.detrend(nfreq, nt, intensity.data, istride, weights.data, wstride, epsilon);
}


static void apply_spline_detrender(io_arr2d intensity, io_arr2d weights, rf_kernels::axis_type axis, int nbins, double epsilon=3.0e-4)
{
    if ((intensity.shape(0) != weights.shape(0)) || (intensity.shape(1) != weights.shape(1)))
	throw runtime_error("rf_pipelines.apply_spline_detrender: 'intensity' and 'weights' arrays do not have the same shape");

    if (axis != rf_kernels::AXIS_FREQ)
	throw runtime_error("rf_pipelines.apply_spline_detrender: currently, only AXIS_FREQ is supported");
    
    int nfreq = intensity.shape(0);
    int nt = intensity.shape(1);
    int istride = xdiv(intensity.stride(0), sizeof(float));
    int wstride = xdiv(weights.stride(0), sizeof(float));

    // Additional argument checks will be performed in kernel constructor.
    rf_kernels::spline_detrender kernel(nfreq, nbins, epsilon);

    kernel.detrend(nt, intensity.data, istride, weights.data, wstride);
}


static void apply_intensity_clipper(in_arr2d intensity, io_arr2d weights, rf_kernels::axis_type axis, double sigma, int niter=1, double iter_sigma=0, int Df=1, int Dt=1, bool two_pass=false)
{
    if ((intensity.shape(0) != weights.shape(0)) || (intensity.shape(1) != weights.shape(1)))
	throw runtime_error("rf_pipelines.apply_intensity_clipper: 'intensity' and 'weights' arrays do not have the same shape");

    int nfreq = intensity.shape(0);
    int nt = intensity.shape(1);
    int istride = xdiv(intensity.stride(0), sizeof(float));
    int wstride = xdiv(weights.stride(0), sizeof(float));
    
    rf_kernels::intensity_clipper ic(nfreq, nt, axis, sigma, Df, Dt, niter, iter_sigma, two_pass);
    ic.clip(intensity.data, istride, weights.data, wstride);
}


static void apply_std_dev_clipper(in_arr2d intensity, io_arr2d weights, rf_kernels::axis_type axis, double sigma, int Df=1, int Dt=1, bool two_pass=false)
{
    if ((intensity.shape(0) != weights.shape(0)) || (intensity.shape(1) != weights.shape(1)))
	throw runtime_error("rf_pipelines.apply_std_dev_clipper: 'intensity' and 'weights' arrays do not have the same shape");

    int nfreq = intensity.shape(0);
    int nt = intensity.shape(1);
    int istride = xdiv(intensity.stride(0), sizeof(float));
    int wstride = xdiv(weights.stride(0), sizeof(float));
    
    rf_kernels::std_dev_clipper sd(nfreq, nt, axis, sigma, Df, Dt, two_pass);
    sd.clip(intensity.data, istride, weights.data, wstride);
}				  


static py_tuple wi_downsample(const in_arr2d &in_intensity, const in_arr2d &in_weights, int Df, int Dt)
{
    if ((in_intensity.shape(0) != in_weights.shape(0)) || (in_intensity.shape(1) != in_weights.shape(1)))
	throw runtime_error("rf_pipelines.wi_downsample: 'intensity' and 'weights' arrays do not have the same shape");

    int in_nfreq = in_intensity.shape(0);
    int in_nt = in_intensity.shape(1);
    int in_istride = xdiv(in_intensity.stride(0), sizeof(float));
    int in_wstride = xdiv(in_weights.stride(0), sizeof(float));
    
    if ((in_nfreq <= 0) || (in_nt <= 0))
	throw runtime_error("rf_pipelines.wi_downsample: (in_nfreq,in_nt)=(" + to_string(in_nfreq) + "," + to_string(in_nt) + ") is invalid");

    if ((Df <= 0) || (Dt <= 0))
	throw runtime_error("rf_pipelines.wi_downsample: (Df,Dt)=(" + to_string(Df) + "," + to_string(Dt) + ") is invalid");

    if (in_nfreq % Df)
	throw runtime_error("rf_pipelines.wi_downsample: in_nfreq=" + to_string(in_nfreq) + " is not divisible by Df=" + to_string(Df));

    if (in_nt % Dt)
	throw runtime_error("rf_pipelines.wi_downsample: in_nt=" + to_string(in_nfreq) + " is not divisible by Dt=" + to_string(Dt));

    int out_nfreq = xdiv(in_nfreq, Df);
    int out_nt = xdiv(in_nt, Dt);

    // FIXME a little awkward.  Maybe pyclops::io_array<float>::make_2d(nfreq,nt) would be a good interface?
    npy_intp out_shape[2] = { out_nfreq, out_nt };
    py_array out_intensity = py_array::make(2, out_shape, npy_type<float>::id);
    py_array out_weights = py_array::make(2, out_shape, npy_type<float>::id);
    
    rf_kernels::wi_downsampler ds(Df, Dt);
    
    ds.downsample(out_nfreq, out_nt,
		  (float *)out_intensity.data(), out_nt,
		  (float *)out_weights.data(), out_nt,
		  in_intensity.data, in_istride,
		  in_weights.data, in_wstride);

    return py_tuple::make(out_intensity, out_weights);
}


static void weight_upsample(io_arr2d w_out, const in_arr2d &w_in, float w_cutoff=0.0)
{
    int nfreq_out = w_out.shape(0);
    int nt_out = w_out.shape(1);
    int nfreq_in = w_in.shape(0);
    int nt_in = w_in.shape(1);

    if (nfreq_out==0 || nt_out==0 || nfreq_in==0 || nt_in==0)
	throw runtime_error("rf_pipelines.weight_upsample: all array dimensions must be nonzero");

    if ((nfreq_out % nfreq_in) || (nt_out % nt_in))
	throw runtime_error("rf_pipelines.weight_upsample: output array dimensions must be multiples of input array dimensions");

    int Df = xdiv(nfreq_out, nfreq_in);
    int Dt = xdiv(nt_out, nt_in);
    int ostride = xdiv(w_out.stride(0), sizeof(float));
    int istride = xdiv(w_in.stride(0), sizeof(float));

    rf_kernels::weight_upsampler u(Df, Dt);
    u.upsample(nfreq_in, nt_in, w_out.data, ostride, w_in.data, istride, w_cutoff);
}


static py_tuple weighted_mean_and_rms(in_arr2d intensity, in_arr2d weights, int niter=1, double sigma=3.0, bool two_pass=false)
{
    if ((intensity.shape(0) != weights.shape(0)) || (intensity.shape(1) != weights.shape(1)))
	throw runtime_error("rf_pipelines.weighted_mean_and_rms: 'intensity' and 'weights' arrays do not have the same shape");

    int nfreq = intensity.shape(0);
    int nt = intensity.shape(1);
    int istride = xdiv(intensity.stride(0), sizeof(float));
    int wstride = xdiv(weights.stride(0), sizeof(float));
    
    rf_kernels::weighted_mean_rms w(nfreq, nt, rf_kernels::AXIS_NONE, 1, 1, niter, sigma, two_pass);
    w.compute_wrms(intensity.data, istride, weights.data, wstride);

    double mean = w.out_mean[0];
    double rms = w.out_rms[0];
    
    return py_tuple::make(mean, rms);
}


static void wrap_kernels(extension_module &m)
{
    string doc_pd = ("apply_polynomial_detrender(intensity, weights, axis, polydeg, epsilon=1.0e-2)\n"
		     "\n"
		     "Detrends along the specified axis by subtracting a best-fit polynomial.\n"
		     "axis=0 means 'detrend in time', axis=1 means 'detrend in frequency'.\n"
		     "\n"
		     "If the fit is poorly conditioned then the entire frequency channel will be masked\n"
		     "(by setting its weights to zero).  The threshold is controlled by the parameter\n"
		     "'epsilon'.  I think that 1.0e-2 is a reasonable default here, but haven't\n"
		     "experimented systematically.\n");

    string doc_sd = ("apply_spline_detrender(intensity, weights, axis, nbins, epsilon=3.0e-4)\n");

    string doc_ic = ("apply_intensity_clipper(intensity, weights, axis, sigma, niter=1, iter_sigma=0.0, Df=1, Dt=1, two_pass=False)\n"
		     "\n"
		     "'Clips' an array by masking outlier intensities.\n"
		     "The masking is performed by setting elements of the weights array to zero.\n"
		     "\n"
		     "The 'sigma' argument is the threshold (in sigmas from the mean) for clipping.  Note\n"
		     "that the weights are used when calculating both the mean and rms intensity.\n"
		     "\n"
		     "The (Df,Dt) args are downsampling factors on the frequency/time axes.\n"
		     "If no downsampling is desired, set Df=Dt=1.\n"
		     "\n"
		     "The 'axis' argument has the following meaning:\n"
		     "   axis=0      clip along frequency axis, with an outer loop over time samples\n"
		     "   axis=1      clip along time axis, with an outer loop over frequency samples\n"
		     "   axis=None   2-d clipper\n"
		     "\n"
		     "If niter > 1, then the mean/rms intensity will be computed using iterated clipping,\n"
		     "with threshold 'iter_sigma'.  If the 'iter_sigma' argument is zero, then it defaults\n"
		     "to 'sigma', but the two thresholds need not be the same.\n"
		     "\n"
		     "If the 'two_pass' flag is set, then a more numerically stable but slightly slower algorithm will be used.\n");

    string doc_sdc = ("apply_std_dev_clipper(intensity, weights, axis, sigma, Df=1, Dt=1, two_pass=False)\n"
		      "\n"
		      "'Clips' an array by masking rows/columns whose standard deviation is an outlier.\n"
		      "The masking is performed by setting elements of the weights array to zero.\n"
		      "\n"
		      "The 'axis' argument has the following meaning:\n"
		      "   axis=0   clip time samples whose variance in frequency is high\n"
		      "   axis=1   clip frequency channels whose variance in time is high\n"
		      "\n"
		      "The (Df,Dt) args are downsampling factors on the frequency/time axes.\n"
		      "If no downsampling is desired, set Df=Dt=1.\n"
		      "\n"
		      "The 'sigma' argument is the threshold (in sigmas from the mean) for clipping.\n"
		      "If the 'two_pass' flag is set, then a more numerically stable but slightly slower algorithm will be used.\n");
    
    string doc_ds = ("wi_downsample(intensity, weights, Df, Dt) -> (intensity, weights)\n"
		     "\n"
		     "Downsamples a weighted intensity array, and returns a new pair (intensity, weights).\n"
		     "The downsampling factors (Df,Dt) must be powers of two.\n"
		     "\n"
		     "Note that the normalization of the downsampled 'weights' array differs\n"
		     "(by a factor of Df*Dt) from the python version of wi_downsample().\n");

    string doc_us = ("weight_upsample(w_out, w_in, w_cutoff=0.0)\n"
		     "\n"
		     "Upsamples the low-res weight array 'w_in', and updates the hi-res weight array 'w_out' in place.\n"
		     "The array dimensions of 'w_in' must evenly divide the dimensions of 'w_out'.\n");
    
    string doc_mr = ("weighted_mean_and_rms(intensity, weights, niter=1, sigma=3.0, two_pass=False)\n"
		     "\n"
		     "Computes weighted mean/rms of a 2D intensity array.\n"
		     "If the 'niter' argument is >1, then the calculation will be iterated, clipping\n"
		     "outlier samples which differ from the mean by the specified number of \"sigmas\".\n"
		     "\n"
		     "If the 'two_pass' flag is set, a more numerically stable but slightly slower algorithm will be used.\n");

    
    m.add_function("apply_polynomial_detrender", doc_pd,
		   wrap_func(apply_polynomial_detrender, "intensity", "weights", "axis", "polydeg", kwarg("epsilon",1.0e-2)));

    m.add_function("apply_spline_detrender", doc_sd,
		   wrap_func(apply_spline_detrender, "intensity", "weights", "axis", "nbins", kwarg("epsilon",3.0e-4)));

    m.add_function("apply_intensity_clipper", doc_ic,
		   wrap_func(apply_intensity_clipper, "intensity", "weights", "axis", "sigma", kwarg("niter",1),
			     kwarg("iter_sigma",0.0), kwarg("Df",1), kwarg("Dt",1), kwarg("two_pass",false)));

    m.add_function("apply_std_dev_clipper", doc_sdc,
		   wrap_func(apply_std_dev_clipper, "intensity", "weights", "axis", "sigma",
			     kwarg("Df",1), kwarg("Dt",1), kwarg("two_pass",false)));

    m.add_function("wi_downsample", doc_ds,
		   wrap_func(wi_downsample, "intensity", "weights", "Df", "Dt"));

    m.add_function("weight_upsample", doc_us,
		   wrap_func(weight_upsample, "intensity", "weights", kwarg("w_cutoff",0.0)));

    m.add_function("weighted_mean_and_rms", doc_mr,
		   wrap_func(weighted_mean_and_rms, "intensity", "weights", kwarg("niter",1),
			     kwarg("sigma",3.0), kwarg("two_pass",false)));
}

// -------------------------------------------------------------------------------------------------
//
// wrap_bonsai().


static void wrap_bonsai(extension_module &m)
{
    string doc_bd = ("bonsai_dedisperser(config_filename, fill_rfi_mask=False, verbosity=1) -> pipeline_object\n"
		     "\n"
		     "A \"transform\" which doesn't actually modify the data, it just runs the bonsai dedisperser (C++ version)\n"
		     "\n"
		     "FIXME: currently, there are two versions of the bonsai_dedisperser, written in python and C++.\n"
		     "From python, they are constructed as 'bonsai_dedisperser' and 'bonsai_dedisperser_cpp' respectively.\n"
		     "In the pipeline json output, they are represented as 'bonsai_dedisperser_python' and 'bonsai_dedisperser_cpp'.\n"
		     "The two versions of the bonsai_dedisperser will be combined eventually!\n");

    std::function<shared_ptr<pipeline_object>(const string &, bool, int)>
	f_bd = [](const string &config_filename, bool fill_rfi_mask, int verbosity) -> shared_ptr<pipeline_object>
	{
	    bonsai_initializer ini_params;
	    ini_params.fill_rfi_mask = fill_rfi_mask;
	    ini_params.verbosity = verbosity;
	    
	    return make_bonsai_dedisperser(config_filename, ini_params);
	};

    m.add_function("bonsai_dedisperser_cpp", doc_bd, wrap_func(f_bd, "config_filename", kwarg("fill_rfi_mask",false), kwarg("verbosity",1)));
}


// -------------------------------------------------------------------------------------------------


PyMODINIT_FUNC initrf_pipelines_c(void)
{
    import_array();

    extension_module m("rf_pipelines_c", "rf_pipelines_c: a C++ library containing low-level rf_pipelines code");

    // rf_pipelines_base_classes.hpp
    wrap_pipeline_object(m);
    wrap_chunked_pipeline_object(m);
    wrap_wi_stream(m);
    wrap_wi_transform(m);

    // rf_pipelines_inventory.hpp
    wrap_containers(m);
    wrap_utility_classes(m);
    wrap_chime_streams(m);
    wrap_chime_ostreams(m);
    wrap_misc_streams(m);
    wrap_detrenders(m);
    wrap_clippers(m);
    wrap_mask_counters(m);
    wrap_kernels(m);
    wrap_bonsai(m);

    m.finalize();
}
