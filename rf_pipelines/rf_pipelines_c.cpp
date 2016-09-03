// Note: I haven't systematically documented the C++ interface to rf_pipelines,
// so the level of documentation will be hit-or-miss.  Also please note that the
// python-wrapping in rf_pipelines_c.cpp is kind of a mess which I hope to improve
// soon.  In the meantime if you want to python-wrap a C++ class, just email me
// and I'll help navigate the mess!

#include "python_extension_helpers.hpp"
#include "rf_pipelines.hpp"

using namespace std;

static PyObject *make_temporary_stream(const rf_pipelines::wi_stream &s);


// -------------------------------------------------------------------------------------------------
//
// wi_transform wrapper class


// "Upcalling" transform whose virtual functions are implemented by python upcalls.
struct upcalling_wi_transform : public rf_pipelines::wi_transform
{
    object weakref;

    upcalling_wi_transform(PyObject *self) :
	weakref(PyWeakref_NewRef(self,NULL), false)
    { }

    virtual ~upcalling_wi_transform() { }

    // returns borrowed reference
    PyObject *get_pyobj() const
    {
	PyObject *ret = PyWeakref_GetObject(weakref.ptr);
	if (!ret)
	    throw python_exception();
	if (ret == Py_None)
	    throw runtime_error("rf_pipelines: internal error: weak reference expired [should never happen]");
	return ret;
    }

    // helper function
    static object array2d_to_python(ssize_t m, ssize_t n, const float *src, ssize_t src_stride)
    {
	if ((m <= 0) || (n <= 0) || !src)
	    throw runtime_error("rf_pipelines: array2d_to_python: internal checks failed [should never happen]");

	npy_intp dims[2] = { m, n };
	npy_intp strides[2] = { src_stride * (ssize_t)sizeof(float), (ssize_t)sizeof(float) };

	// PyArray_New(subtype, nd, dims, type_num, npy_intp* strides, void* data, int itemsize, int flags, PyObject* obj)
	PyObject *p = PyArray_New(&PyArray_Type, 2, dims, NPY_FLOAT, strides, (void *)src, 0, NPY_ARRAY_WRITEABLE, NULL);

	return object(p, false);
    }

    virtual void set_stream(const rf_pipelines::wi_stream &stream) override
    {
	PyObject *sp = make_temporary_stream(stream);
	object s(sp, false);

	PyObject *retp = PyObject_CallMethod(this->get_pyobj(), (char *)"set_stream", (char *)"O", sp);
	object ret(retp, false);  // a convenient way to ensure Py_DECREF gets called, and throw an exception on failure

	if (s.get_refcount() > 1)
	    throw runtime_error("fatal: wi_transform.set_stream() callback kept a reference to the stream");
    }

    virtual void start_substream(int isubstream, double t0) override
    {	
	PyObject *p = PyObject_CallMethod(this->get_pyobj(), (char *)"start_substream", (char *)"id", isubstream, t0);
	object ret(p, false);  // a convenient way to ensure Py_DECREF gets called, and throw an exception on failure
    }

    virtual void process_chunk(double t0, double t1, float *intensity, float *weights, ssize_t stride, float *pp_intensity, float *pp_weights, ssize_t pp_stride) override
    {
	object np_intensity = array2d_to_python(nfreq, nt_chunk + nt_postpad, intensity, stride);
	object np_weights = array2d_to_python(nfreq, nt_chunk + nt_postpad, weights, stride);
	object np_pp_intensity;
	object np_pp_weights;

	if (nt_prepad > 0) {
	    np_pp_intensity = array2d_to_python(nfreq, nt_prepad, pp_intensity, pp_stride);
	    np_pp_weights = array2d_to_python(nfreq, nt_prepad, pp_weights, pp_stride);
	}

	// FIXME: a weird corner case that I'd like to understand more generally: control-C can
	// cause the np_intensity refcount to equal 2 when PyObject_CallMethod() returns.  Does
	// this means it leaks memory?
	PyObject *p = PyObject_CallMethod(this->get_pyobj(), (char *)"process_chunk", (char *)"ddOOOO", t0, t1,
					  np_intensity.ptr, np_weights.ptr, np_pp_intensity.ptr, np_pp_weights.ptr);

	object ret(p, false);  // a convenient way to ensure Py_DECREF gets called, and throw an exception on failure

	if (np_intensity.get_refcount() > 1)
	    throw runtime_error("fatal: wi_transform.process_chunk() callback kept a reference to the 'intensity' array");
	if (np_weights.get_refcount() > 1)
	    throw runtime_error("fatal: wi_transform.process_chunk() callback kept a reference to the 'weights' array");
	if ((nt_prepad > 0) && (np_pp_intensity.get_refcount() > 1))
	    throw runtime_error("fatal: wi_transform.process_chunk() callback kept a reference to the 'pp_intensity' array");
	if ((nt_prepad > 0) && (np_pp_weights.get_refcount() > 1))
	    throw runtime_error("fatal: wi_transform.process_chunk() callback kept a reference to the 'pp_weights' array");
    }

    virtual void end_substream() override
    {
	PyObject *p = PyObject_CallMethod(this->get_pyobj(), (char *)"end_substream", NULL);
	object ret(p, false);
    }

    virtual string get_name() const override
    {
	PyObject *sobj = PyObject_Str(this->get_pyobj());
	object ret(sobj, false);

	// Returns a pointer to an internal buffer, not a copy, so no free() necessary.
	char *s = PyString_AsString(sobj);
	if (!s)
	    throw python_exception();

	return string(s);
    }
};


struct wi_transform_object {
    PyObject_HEAD

    shared_ptr<rf_pipelines::wi_transform> *pshared;

    // forward declarations (these guys need to come after 'wi_transform_type')
    static PyObject *make(const shared_ptr<rf_pipelines::wi_transform> &p);
    static bool isinstance(PyObject *obj);
    
    static PyObject *tp_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
    {
	PyObject *self_ = type->tp_alloc(type, 0);
	if (!self_)
	    return NULL;

	wi_transform_object *self = (wi_transform_object *) self_;
	self->pshared = nullptr;
	return self_;
    }

    static void tp_dealloc(PyObject *self_)
    {
	wi_transform_object *self = (wi_transform_object *)self_;

	delete self->pshared;
	self->pshared = nullptr;
	Py_TYPE(self)->tp_free(self_);
    }

    // Helper for get_pshared(): get a pointer-to-shared-ptr from a (PyObject *) which is known to be a wi_transform object.
    // This is the where the upcalling transform gets created, if it doesn't already exist.
    static inline shared_ptr<rf_pipelines::wi_transform> get_pshared(PyObject *self)
    {
	if (!wi_transform_object::isinstance(self))
	    throw runtime_error("rf_pipelines: 'self' argument to wi_transform method was not an object of type wi_transform");

	wi_transform_object *t = (wi_transform_object *) self;

	if (t->pshared)
	    return *(t->pshared);

	// FIXME this way of instantiating the upcalling transform is awkward and could be improved
	shared_ptr<rf_pipelines::wi_transform> p = make_shared<upcalling_wi_transform> (self);
	t->pshared = new shared_ptr<rf_pipelines::wi_transform> (p);
	return p;
    }

    // Get a bare pointer from a (PyObject *) which is known to be a wi_transform_object
    static inline rf_pipelines::wi_transform *get_pbare(PyObject *self)
    {
	// note: get_pshared() checks that object is a wi_transform
	shared_ptr<rf_pipelines::wi_transform> p = get_pshared(self);
	
	if (!p)
	    throw runtime_error("rf_pipelines: internal error: empty shared_ptr<> in wi_transform_object [should never happen]");

	return p.get();
    }


    static PyObject *nfreq_getter(PyObject *self, void *closure)
    {
	return Py_BuildValue("i", get_pbare(self)->nfreq);
    }

    static int nfreq_setter(PyObject *self, PyObject *value, void *closure)
    {
	get_pbare(self)->nfreq = ssize_t_from_python(value);
	return 0;
    }

    static PyObject *nt_chunk_getter(PyObject *self, void *closure)
    {
	return Py_BuildValue("i", get_pbare(self)->nt_chunk);
    }

    static int nt_chunk_setter(PyObject *self, PyObject *value, void *closure)
    {
	get_pbare(self)->nt_chunk = ssize_t_from_python(value);
	return 0;
    }

    static PyObject *nt_prepad_getter(PyObject *self, void *closure)
    {
	return Py_BuildValue("i", get_pbare(self)->nt_prepad);
    }

    static int nt_prepad_setter(PyObject *self, PyObject *value, void *closure)
    {
	get_pbare(self)->nt_prepad = ssize_t_from_python(value);
	return 0;
    }

    static PyObject *nt_postpad_getter(PyObject *self, void *closure)
    {
	return Py_BuildValue("i", get_pbare(self)->nt_postpad);
    }

    static int nt_postpad_setter(PyObject *self, PyObject *value, void *closure)
    {
	get_pbare(self)->nt_postpad = ssize_t_from_python(value);
	return 0;
    }

    // int add_plot_group(const std::string &name, int nt_per_pix, int ny)
    static PyObject *add_plot_group(PyObject *self, PyObject *args, PyObject *kwds)
    {
	static const char *kwlist[] = { "name", "nt_per_pix", "ny", NULL };

	char *name = nullptr;
	int nt_per_pix = 0;
	int ny = 0;

	if (!PyArg_ParseTupleAndKeywords(args, kwds, "sii", (char **)kwlist, &name, &nt_per_pix, &ny))
	    return NULL;

	int ret = get_pbare(self)->add_plot_group(name, nt_per_pix, ny);
	return Py_BuildValue("i", ret);
    }

    static constexpr const char *add_plot_group_docstring = 
	"Usage: add_plot_group(name, nt_per_pix, ny) -> integer\n\n"
	"    Each transform's output plots are divided into one or more \"plot groups\"\n"
	"    For example, the bonsai dedisperser can write one plot group per internally defined tree.\n"
	"    The 'nt_per_pix' arg is the number of pipeline time samples per x-pixel in the plot.\n"
	"    The 'ny' arg is the number of y-pixels (assumed to be the same for all plots in the group).\n"
	"    The return value is the group_id arg needed in add_plot(), and group_ids always go 0,1,...\n";

    // string add_plot(const string &basename, int64_t it0, int nt, int nx, int ny, int group_id=0)
    static PyObject *add_plot(PyObject *self, PyObject *args, PyObject *kwds)
    {
	static const char *kwlist[] = { "basename", "it0", "nt", "nx", "ny", "group_id", NULL };

	char *basename = nullptr;
	ssize_t it0 = 0;
	int nt = 0;
	int nx = 0;
	int ny = 0;
	int group_id = 0;

	if (!PyArg_ParseTupleAndKeywords(args, kwds, "sniii|i", (char **)kwlist, &basename, &it0, &nt, &nx, &ny, &group_id))
	    return NULL;
	
	string ret = get_pbare(self)->add_plot(basename, it0, nt, nx, ny, group_id);
	return Py_BuildValue("s", ret.c_str());   // Note: Py_BuildValue() copies the string
    }

    static constexpr const char *add_plot_docstring = 
	"Usage: add_plot(basename, int64_t it0, int nt, int nx, int ny, int group_id=0) -> string\n\n"
	"Call just before writing a plot.\n"
	"    The range of time samples in the plot is [it0:it0+nt).\n"
	"    The pixel dimensions of the plot are (nx,ny).  These are redundant since they can be deduced\n"
	"    from (it0,nt) but we use them for error checking.\n"
	"    The return value is the full pathname ('basename' with the stream output_dir prepended)\n";

    // string add_file(const string &basename)
    static PyObject *add_file(PyObject *self, PyObject *args, PyObject *kwds)
    {
	static const char *kwlist[] = { "basename", NULL };

	char *basename = nullptr;

	if (!PyArg_ParseTupleAndKeywords(args, kwds, "s", (char **)kwlist, &basename))
	    return NULL;
	
	string ret = get_pbare(self)->add_file(basename);
	return Py_BuildValue("s", ret.c_str());   // Note: Py_BuildValue() copies the string
    }

    static constexpr const char *add_file_docstring = 
	"Usage: add_file(basename) -> string\n\n"
	"    Call just before writing a non-plot file, to check for filename collisions between transforms.\n"
	"    The return value is the full pathname ('basename' with stream output_dir prepended)\n";


    static constexpr const char *dummy_docstring = 
	"wi_transform is a C++ base class, and transforms written in C++ inherit from it.\n"
	"Transforms written in python will inherit from rf_pipelines.py_wi_transform.\n"
	"For documentation of wi_transform and its methods, see the rf_pipelines.py_wi_transform docstring.";
};


static PyGetSetDef wi_transform_getseters[] = {
    { (char *)"nfreq", 
      tc_wrap_getter<wi_transform_object::nfreq_getter>, 
      tc_wrap_setter<wi_transform_object::nfreq_setter>, 
      (char *)wi_transform_object::dummy_docstring, NULL },

    { (char *)"nt_chunk", 
      tc_wrap_getter<wi_transform_object::nt_chunk_getter>, 
      tc_wrap_setter<wi_transform_object::nt_chunk_setter>,
      (char *)wi_transform_object::dummy_docstring, NULL },

    { (char *)"nt_prepad", 
      tc_wrap_getter<wi_transform_object::nt_prepad_getter>, 
      tc_wrap_setter<wi_transform_object::nt_prepad_setter>,
      (char *)wi_transform_object::dummy_docstring, NULL },

    { (char *)"nt_postpad", 
      tc_wrap_getter<wi_transform_object::nt_postpad_getter>, 
      tc_wrap_setter<wi_transform_object::nt_postpad_setter>,
      (char *)wi_transform_object::dummy_docstring, NULL },

    { NULL, NULL, NULL, NULL, NULL }
};


static PyMethodDef wi_transform_methods[] = {
    { "add_plot_group", (PyCFunction) tc_wrap3<wi_transform_object::add_plot_group>, METH_VARARGS | METH_KEYWORDS, (char *)wi_transform_object::add_plot_group_docstring },
    { "add_plot", (PyCFunction) tc_wrap3<wi_transform_object::add_plot>, METH_VARARGS | METH_KEYWORDS, (char *)wi_transform_object::add_plot_docstring },
    { "add_file", (PyCFunction) tc_wrap3<wi_transform_object::add_file>, METH_VARARGS | METH_KEYWORDS, (char *)wi_transform_object::add_file_docstring },
    { NULL, NULL, 0, NULL }
};


static PyTypeObject wi_transform_type = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "rf_pipelines_c.wi_transform",  /* tp_name */
    sizeof(wi_transform_object),    /* tp_basicsize */
    0,                         /* tp_itemsize */
    wi_transform_object::tp_dealloc,  /* tp_dealloc */
    0,                         /* tp_print */
    0,                         /* tp_getattr */
    0,                         /* tp_setattr */
    0,                         /* tp_reserved */
    0,                         /* tp_repr */
    0,                         /* tp_as_number */
    0,                         /* tp_as_sequence */
    0,                         /* tp_as_mapping */
    0,                         /* tp_hash  */
    0,                         /* tp_call */
    0,                         /* tp_str */
    0,                         /* tp_getattro */
    0,                         /* tp_setattro */
    0,                         /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,   /* tp_flags */
    wi_transform_object::dummy_docstring,       /* tp_doc */
    0,                         /* tp_traverse */
    0,                         /* tp_clear */
    0,                         /* tp_richcompare */
    0,                         /* tp_weaklistoffset */
    0,                         /* tp_iter */
    0,                         /* tp_iternext */
    wi_transform_methods,      /* tp_methods */
    0,                         /* tp_members */
    wi_transform_getseters,    /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    0,                         /* tp_init */
    0,                         /* tp_alloc */
    wi_transform_object::tp_new,  /* tp_new */
};


// static member function
PyObject *wi_transform_object::make(const shared_ptr<rf_pipelines::wi_transform> &ptr)
{
    if (!ptr)
	throw runtime_error("rf_pipelines: internal error: empty pointer passed to make_wi_transform()");

    PyObject *ret_ = wi_transform_object::tp_new(&wi_transform_type, NULL, NULL);
    if (!ret_)
	return NULL;

    wi_transform_object *ret = (wi_transform_object *) (ret_);
    ret->pshared = new shared_ptr<rf_pipelines::wi_transform> (ptr);
    return ret_;
}

// static member function
bool wi_transform_object::isinstance(PyObject *obj)
{
    return PyObject_IsInstance(obj, (PyObject *) &wi_transform_type);
}


// -------------------------------------------------------------------------------------------------
//
// exception_monitor: a dummy transform which regularly checks the python exception state
//   whenever process_chunk() is called.  This ensures that control-C always gets caught,
//   in the case where we're running rf_pipelines through the python interpreter, but all 
//   streams and transforms are C++ classes.


struct exception_monitor : public rf_pipelines::wi_transform
{
    exception_monitor(ssize_t nt_chunk_)
    {
	this->nt_chunk = nt_chunk_;
    }

    virtual ~exception_monitor() { }

    virtual void set_stream(const rf_pipelines::wi_stream &stream) override
    {
	this->nfreq = stream.nfreq;
    }

    virtual void start_substream(int isubstream, double t0) override { }

    virtual void process_chunk(double t0, double t1, float *intensity, float *weights, ssize_t stride, float *pp_intensity, float *pp_weights, ssize_t pp_stride) override
    {
	if (PyErr_Occurred() || PyErr_CheckSignals())
	    throw python_exception();
    }

    virtual void end_substream() override { }
    
    // should never be called
    virtual string get_name() const override { return "exception_monitor"; }
};


// -------------------------------------------------------------------------------------------------
//
// wi_run_state wrapper class


struct wi_run_state_object {
    PyObject_HEAD

    // "Borrowed" reference, cannot be NULL
    rf_pipelines::wi_run_state *pbare;

    // Forward declarations (these guys need to come after 'wi_run_state_type')
    static PyObject *make(rf_pipelines::wi_run_state &rs);
    static bool isinstance(PyObject *obj);


    static inline rf_pipelines::wi_run_state *get_pbare(PyObject *obj)
    {
	if (!wi_run_state_object::isinstance(obj))
	    throw runtime_error("rf_pipelines: 'self' argument to wi_run_state method was not an object of type wi_run_state");

	rf_pipelines::wi_run_state *ret = ((wi_run_state_object *) obj)->pbare;
	if (!ret)
	    throw runtime_error("rf_pipelines: wi_run_state object cannot be constructed directly from Python");

	return ret;
    }
    
    static PyObject *tp_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
    {
	PyObject *self_ = type->tp_alloc(type, 0);
	if (!self_)
	    return NULL;

	wi_run_state_object *self = (wi_run_state_object *) self_;
	self->pbare = nullptr;
	return self_;
    }

    // void start_substream(double t0);
    // Note: this one is METH_O
    static PyObject *start_substream(PyObject *self, PyObject *arg)
    {
	double t0 = double_from_python(arg);
	get_pbare(self)->start_substream(t0);

	Py_INCREF(Py_None);
	return Py_None;
    }

    //
    // From rf_pipelines/__init__.py:
    //    "Comment: the python stream API has an extra copy relative to the C++ API, so python streams may be a little
    //     slower than C++ streams, but it's hard to do anything about this!"
    //
    // void write(intensity, weight, t0=None)
    //
    static PyObject *write(PyObject *self, PyObject *args, PyObject *kwds)
    {
	rf_pipelines::wi_run_state *run_state = get_pbare(self);
	PyObject *src_intensity_obj = Py_None;
	PyObject *src_weight_obj = Py_None;
	PyObject *t0_obj = Py_None;

	static const char *kwlist[] = { "intensity", "weight", "t0", NULL };

	// Note: the object pointers will be borrowed references
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "OO|O", (char **)kwlist, &src_intensity_obj, &src_weight_obj, &t0_obj))
            return NULL;

	// Note that 't0' is meaningful if and only if (t0_obj != Py_None)
	double t0 = (t0_obj != Py_None) ? double_from_python(t0_obj) : 0.0;

	src_intensity_obj = PyArray_FromAny(src_intensity_obj, NULL, 2, 2, 0, NULL);
	object ref1(src_intensity_obj, false);   // manages refcount, throws exception on NULL

	src_weight_obj = PyArray_FromAny(src_weight_obj, NULL, 2, 2, 0, NULL);
	object ref2(src_weight_obj, false);     // manages refcount, throws exception on NULL

	PyArrayObject *src_intensity = (PyArrayObject *) src_intensity_obj;
	PyArrayObject *src_weight = (PyArrayObject *) src_weight_obj;

	// should never happen, since min_depth=max_depth=2 was specified in PyArray_FromAny()
	if ((PyArray_NDIM(src_intensity) != 2) || (PyArray_NDIM(src_weight) != 2))
	    throw runtime_error("ndim != 2 in rf_pipelines.wi_run_state.write()?! (should never happen)");

	npy_intp *src_intensity_shape = PyArray_SHAPE(src_intensity);
	npy_intp *src_weight_shape = PyArray_SHAPE(src_intensity);	

	if ((src_intensity_shape[0] != src_weight_shape[0]) || (src_intensity_shape[1] != src_weight_shape[1]))
	    throw runtime_error("rf_pipelines.wi_run_state.write(): 'intensity' and 'weight' arrays must have the same shape");

	if (src_intensity_shape[0] != run_state->nfreq)
	    throw runtime_error("rf_pipelines.wi_run_state.write(): first dimension of 'intensity' and 'weight' arrays must equal stream's 'nfreq'");
	if (src_intensity_shape[1] <= 0)
	    throw runtime_error("rf_pipelines.wi_run_state.write(): zero-length write (currently treated as an error");
	if (src_intensity_shape[1] > run_state->nt_stream_maxwrite)
	    throw runtime_error("rf_pipelines.wi_run_state.write(): number of time samples written cannot exceed stream's nt_maxwrite");

	ssize_t nt = src_intensity_shape[1];
	float *dst_intensity_ptr = nullptr;
	float *dst_weight_ptr = nullptr;
	ssize_t dst_cstride = 0;
	bool zero_flag = false;

	if (t0_obj == Py_None)
	    run_state->setup_write(nt, dst_intensity_ptr, dst_weight_ptr, dst_cstride, zero_flag);
	else
	    run_state->setup_write(nt, dst_intensity_ptr, dst_weight_ptr, dst_cstride, zero_flag, t0);	    

	// Make destination array objects, in order to call PyArray_CopyInto()	
	// Note syntax is: PyArray_New(subtype, nd, dims, type_num, npy_intp* strides, void* data, int itemsize, int flags, PyObject* obj)

	npy_intp dst_dims[2] = { run_state->nfreq, nt };
	npy_intp dst_strides[2] = { dst_cstride * (ssize_t)sizeof(float), (ssize_t)sizeof(float) };

	PyObject *dst_intensity = PyArray_New(&PyArray_Type, 2, dst_dims, NPY_FLOAT, dst_strides, (void *)dst_intensity_ptr, 0, NPY_ARRAY_WRITEABLE, NULL);
	object ref3(dst_intensity, false);   // manages refcount, throws exception on NULL

	PyObject *dst_weight = PyArray_New(&PyArray_Type, 2, dst_dims, NPY_FLOAT, dst_strides, (void *)dst_weight_ptr, 0, NPY_ARRAY_WRITEABLE, NULL);
	object ref4(dst_weight, false);   // manages refcount, throws exception on NULL

	// Copy data into pipeline buffer

	if (PyArray_CopyInto((PyArrayObject *)dst_intensity, src_intensity))
	    return NULL;
	if (PyArray_CopyInto((PyArrayObject *)dst_weight, src_weight))
	    return NULL;
	
	run_state->finalize_write(nt);

	Py_INCREF(Py_None);
	return Py_None;
    }

    // void end_substream();
    // Note: this one is METH_NOARGS
    static PyObject *end_substream(PyObject *self)
    {
	get_pbare(self)->end_substream();

	Py_INCREF(Py_None);
	return Py_None;
    }

    static constexpr const char *dummy_docstring = 
	"wi_run_state is a C++ helper class which is used to move data from a stream into the rf_pipelines ring buffer.\n"
	"For documentation of class wi_run_state and its methods, see the rf_pipelines.py_wi_stream docstring.";
};


static PyMethodDef wi_run_state_methods[] = {
    { "start_substream", tc_wrap2<wi_run_state_object::start_substream>, METH_O, wi_run_state_object::dummy_docstring },
    { "write", (PyCFunction) tc_wrap3<wi_run_state_object::write>, METH_VARARGS | METH_KEYWORDS, wi_run_state_object::dummy_docstring },
    { "end_substream", (PyCFunction) tc_wrap1<wi_run_state_object::end_substream>, METH_NOARGS, wi_run_state_object::dummy_docstring },
    { NULL, NULL, 0, NULL }
};


static PyTypeObject wi_run_state_type = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "rf_pipelines_c.wi_run_state",  /* tp_name */
    sizeof(wi_run_state_object),    /* tp_basicsize */
    0,                         /* tp_itemsize */
    0,                         /* tp_dealloc */
    0,                         /* tp_print */
    0,                         /* tp_getattr */
    0,                         /* tp_setattr */
    0,                         /* tp_reserved */
    0,                         /* tp_repr */
    0,                         /* tp_as_number */
    0,                         /* tp_as_sequence */
    0,                         /* tp_as_mapping */
    0,                         /* tp_hash  */
    0,                         /* tp_call */
    0,                         /* tp_str */
    0,                         /* tp_getattro */
    0,                         /* tp_setattro */
    0,                         /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,   /* tp_flags */
    wi_run_state_object::dummy_docstring,       /* tp_doc */
    0,                         /* tp_traverse */
    0,                         /* tp_clear */
    0,                         /* tp_richcompare */
    0,                         /* tp_weaklistoffset */
    0,                         /* tp_iter */
    0,                         /* tp_iternext */
    wi_run_state_methods,      /* tp_methods */
    0,                         /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    0,                         /* tp_init */
    0,                         /* tp_alloc */
    wi_run_state_object::tp_new,  /* tp_new */
};

// static member function
PyObject *wi_run_state_object::make(rf_pipelines::wi_run_state &rs)
{
    PyObject *ret_ = wi_run_state_object::tp_new(&wi_run_state_type, NULL, NULL);
    if (!ret_)
	return NULL;

    wi_run_state_object *ret = (wi_run_state_object *) (ret_);
    ret->pbare = &rs;
    return ret_;
}

// static member function
bool wi_run_state_object::isinstance(PyObject *obj)
{
    return PyObject_IsInstance(obj, (PyObject *) &wi_run_state_type);
}


// -------------------------------------------------------------------------------------------------
//
// wi_stream wrapper class


// "Upcalling" stream whose stream_body() virtual function is implemented by python upcalls.
struct upcalling_wi_stream : public rf_pipelines::wi_stream
{
    object weakref;

    upcalling_wi_stream(PyObject *self) :
	weakref(PyWeakref_NewRef(self,NULL), false)
    { }

    virtual ~upcalling_wi_stream() { }

    // returns borrowed reference
    PyObject *get_pyobj()
    {
	PyObject *ret = PyWeakref_GetObject(weakref.ptr);
	if (!ret)
	    throw python_exception();
	if (ret == Py_None)
	    throw runtime_error("rf_pipelines: internal error: weak reference expired [should never happen]");
	return ret;
    }
    
    virtual void stream_body(rf_pipelines::wi_run_state &run_state)
    {
	PyObject *rs = wi_run_state_object::make(run_state);
	object rs_ref(rs, false);

	PyObject *ret = PyObject_CallMethod(this->get_pyobj(), (char *)"stream_body", (char *)"O", rs);
	object ret_ref(ret, false);
	
	if (rs_ref.get_refcount() > 1)
	    throw runtime_error("fatal: wi_stream.stream_body() callback kept a reference to the wi_run_state object");
    }
};


struct wi_stream_object {
    PyObject_HEAD

    // "Borrowed" reference (delete will not be called in tp_dealloc())
    rf_pipelines::wi_stream *pbare;

    // Reference held through shared_ptr.  Using a pointer-to-a-shared pointer was least awkward here.
    shared_ptr<rf_pipelines::wi_stream> *pshared;

    // forward declarations (these guys need to come after 'wi_stream_type'
    static PyObject *make(const shared_ptr<rf_pipelines::wi_stream> &p);
    static bool isinstance(PyObject *obj);
    
    static PyObject *tp_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
    {
	PyObject *self_ = type->tp_alloc(type, 0);
	wi_stream_object *self = (wi_stream_object *) self_;

	self->pbare = nullptr;
	self->pshared = nullptr;
	return self_;
    }

    // note: no tp_alloc function is necessary, since the default (PyType_GenericAlloc()) calls memset().    
    static void tp_dealloc(PyObject *self_)
    {
	wi_stream_object *self = (wi_stream_object *)self_;

	delete self->pshared;

	self->pbare = nullptr;
	self->pshared = nullptr;
	Py_TYPE(self)->tp_free((PyObject*) self);
    }

    // Get a bare pointer from a (PyObject *) which is known to be a wi_stream_object
    static inline rf_pipelines::wi_stream *get_pbare(PyObject *self)
    {
	if (!wi_stream_object::isinstance(self))
	    throw runtime_error("rf_pipelines: 'self' argument to wi_stream method was not an object of type wi_stream");

	wi_stream_object *s = (wi_stream_object *) self;

	if (!s->pbare) {
	    shared_ptr<rf_pipelines::wi_stream> p = make_shared<upcalling_wi_stream> (self);
	    s->pshared = new shared_ptr<rf_pipelines::wi_stream> (p);
	    s->pbare = p.get();
	}

	return s->pbare;
    }
    
    static PyObject *run(PyObject *self, PyObject *args, PyObject *kwds)
    {
	static const char *kwlist[] = { "transforms", "outdir", "noisy", "clobber", "return_json", NULL };

	PyObject *transforms_obj = nullptr;
	const char *outdir = ".";
	int noisy = 1;
	int clobber = 1;
	int return_json = 0;

	if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|siii", (char **)kwlist, &transforms_obj, (char **)&outdir, &noisy, &clobber, &return_json))
	    return NULL;

	rf_pipelines::wi_stream *stream = get_pbare(self);

	PyObject *transforms_iter = PyObject_GetIter(transforms_obj);
	if (!transforms_iter)
	    throw runtime_error("rf_pipelines: expected 'transforms' argument to wi_stream.run() to be a list/iterator");

	object iter_reference(transforms_iter, false);
	vector<object> item_references;

	vector<shared_ptr<rf_pipelines::wi_transform> > transform_list;
	transform_list.push_back(make_shared<exception_monitor> (stream->nt_maxwrite));

	for (;;) {
	    PyObject *item_ptr = PyIter_Next(transforms_iter);
	    if (!item_ptr)
		break;

	    item_references.push_back(object(item_ptr,false));

	    if (!wi_transform_object::isinstance(item_ptr))
		throw runtime_error("rf_pipelines: expected 'transforms' argument to wi_stream.run() to be a list/iterator of wi_transform objects");

	    transform_list.push_back(wi_transform_object::get_pshared(item_ptr));
	}	

	if (PyErr_Occurred())
	    throw python_exception();

	Json::Value json_out;
	Json::Value *json_outp = return_json ? &json_out : nullptr;
	stream->run(transform_list, outdir, json_outp, noisy, clobber);

	if (!return_json) {
	    Py_INCREF(Py_None);
	    return Py_None;
	}

	Json::FastWriter w;
	string ret = w.write(json_out);
	return Py_BuildValue("s", ret.c_str());   // Note: Py_BuildValue() copies the string
    }

    // Properties

    static PyObject *nfreq_getter(PyObject *self, void *closure)
    {
	return Py_BuildValue("i", get_pbare(self)->nfreq);
    }

    static int nfreq_setter(PyObject *self, PyObject *value, void *closure)
    {
	get_pbare(self)->nfreq = ssize_t_from_python(value);
	return 0;
    }

    static PyObject *nt_maxwrite_getter(PyObject *self, void *closure)
    {
	return Py_BuildValue("i", get_pbare(self)->nt_maxwrite);
    }

    static int nt_maxwrite_setter(PyObject *self, PyObject *value, void *closure)
    {
	get_pbare(self)->nt_maxwrite = ssize_t_from_python(value);
	return 0;
    }

    static PyObject *freq_lo_MHz_getter(PyObject *self, void *closure)
    {
	return Py_BuildValue("d", get_pbare(self)->freq_lo_MHz);
    }

    static int freq_lo_MHz_setter(PyObject *self, PyObject *value, void *closure)
    {
	get_pbare(self)->freq_lo_MHz = double_from_python(value);
	return 0;
    }

    static PyObject *freq_hi_MHz_getter(PyObject *self, void *closure)
    {
	return Py_BuildValue("d", get_pbare(self)->freq_hi_MHz);
    }

    static int freq_hi_MHz_setter(PyObject *self, PyObject *value, void *closure)
    {
	get_pbare(self)->freq_hi_MHz = double_from_python(value);
	return 0;
    }

    static PyObject *dt_sample_getter(PyObject *self, void *closure)
    {
	return Py_BuildValue("d", get_pbare(self)->dt_sample);
    }

    static int dt_sample_setter(PyObject *self, PyObject *value, void *closure)
    {
	get_pbare(self)->dt_sample = double_from_python(value);
	return 0;
    }

    static constexpr const char *dummy_docstring = 
	"wi_stream is a C++ base class, and streams written in C++ inherit from it.\n"
	"Streams written in python will inherit from rf_pipelines.py_wi_stream.\n"
	"For documentation of wi_stream and its methods, see the rf_pipelines.py_wi_stream docstring.";
};


static PyMethodDef wi_stream_methods[] = {
    { "run", (PyCFunction) tc_wrap3<wi_stream_object::run>, METH_VARARGS | METH_KEYWORDS, wi_stream_object::dummy_docstring },
    { NULL, NULL, 0, NULL }
};


static PyGetSetDef wi_stream_getseters[] = {
    { (char *)"nfreq", 
      tc_wrap_getter<wi_stream_object::nfreq_getter>, 
      tc_wrap_setter<wi_stream_object::nfreq_setter>, 
      (char *)wi_stream_object::dummy_docstring, NULL },

    { (char *)"nt_maxwrite", 
      tc_wrap_getter<wi_stream_object::nt_maxwrite_getter>, 
      tc_wrap_setter<wi_stream_object::nt_maxwrite_setter>, 
      (char *)wi_stream_object::dummy_docstring, NULL },

    { (char *)"freq_lo_MHz", 
      tc_wrap_getter<wi_stream_object::freq_lo_MHz_getter>, 
      tc_wrap_setter<wi_stream_object::freq_lo_MHz_setter>, 
      (char *)wi_stream_object::dummy_docstring, NULL },

    { (char *)"freq_hi_MHz", 
      tc_wrap_getter<wi_stream_object::freq_hi_MHz_getter>, 
      tc_wrap_setter<wi_stream_object::freq_hi_MHz_setter>, 
      (char *)wi_stream_object::dummy_docstring, NULL },

    { (char *)"dt_sample", 
      tc_wrap_getter<wi_stream_object::dt_sample_getter>, 
      tc_wrap_setter<wi_stream_object::dt_sample_setter>, 
      (char *)wi_stream_object::dummy_docstring, NULL },

    { NULL, NULL, NULL, NULL, NULL }
};


static PyTypeObject wi_stream_type = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "rf_pipelines_c.wi_stream",  /* tp_name */
    sizeof(wi_stream_object),    /* tp_basicsize */
    0,                         /* tp_itemsize */
    wi_stream_object::tp_dealloc,  /* tp_dealloc */
    0,                         /* tp_print */
    0,                         /* tp_getattr */
    0,                         /* tp_setattr */
    0,                         /* tp_reserved */
    0,                         /* tp_repr */
    0,                         /* tp_as_number */
    0,                         /* tp_as_sequence */
    0,                         /* tp_as_mapping */
    0,                         /* tp_hash  */
    0,                         /* tp_call */
    0,                         /* tp_str */
    0,                         /* tp_getattro */
    0,                         /* tp_setattro */
    0,                         /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,   /* tp_flags */
    wi_stream_object::dummy_docstring,          /* tp_doc */
    0,                         /* tp_traverse */
    0,                         /* tp_clear */
    0,                         /* tp_richcompare */
    0,                         /* tp_weaklistoffset */
    0,                         /* tp_iter */
    0,                         /* tp_iternext */
    wi_stream_methods,         /* tp_methods */
    0,                         /* tp_members */
    wi_stream_getseters,       /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    0,                         /* tp_init */
    0,                         /* tp_alloc */
    wi_stream_object::tp_new,  /* tp_new */
};


// static member function
PyObject *wi_stream_object::make(const shared_ptr<rf_pipelines::wi_stream> &ptr)
{
    if (!ptr)
	throw runtime_error("rf_pipelines: internal error: empty pointer passed to wi_stream_object::make()");

    PyObject *ret_ = wi_stream_object::tp_new(&wi_stream_type, NULL, NULL);
    if (!ret_)
	return NULL;

    wi_stream_object *ret = (wi_stream_object *) (ret_);
    ret->pbare = ptr.get();
    ret->pshared = new shared_ptr<rf_pipelines::wi_stream> (ptr);
    return ret_;
}

// FIXME I'd like to make this a static member function
static PyObject *make_temporary_stream(const rf_pipelines::wi_stream &s)
{
    PyObject *ret_ = wi_stream_object::tp_new(&wi_stream_type, NULL, NULL);
    if (!ret_)
	return NULL;

    wi_stream_object *ret = (wi_stream_object *) (ret_);
    ret->pbare = const_cast<rf_pipelines::wi_stream *> (&s);
    return ret_;
}

// static member function
bool wi_stream_object::isinstance(PyObject *obj)
{
    return PyObject_IsInstance(obj, (PyObject *) &wi_stream_type);
}


// -------------------------------------------------------------------------------------------------
//
// Library


static PyObject *make_psrfits_stream(PyObject *self, PyObject *args)
{
    const char *filename = nullptr;
    if (!PyArg_ParseTuple(args, "s", &filename))
	return NULL;

    shared_ptr<rf_pipelines::wi_stream> ret = rf_pipelines::make_psrfits_stream(filename);
    return wi_stream_object::make(ret);
}


static PyObject *make_chime_stream_from_acqdir(PyObject *self, PyObject *args)
{
    const char *filename = nullptr;
    ssize_t nt_chunk = 0;
    ssize_t noise_source_align = 0;

    if (!PyArg_ParseTuple(args, "snn", &filename, &nt_chunk, &noise_source_align))
	return NULL;

    shared_ptr<rf_pipelines::wi_stream> ret = rf_pipelines::make_chime_stream_from_acqdir(filename, nt_chunk, noise_source_align);
    return wi_stream_object::make(ret);    
}


static PyObject *make_chime_stream_from_filename(PyObject *self, PyObject *args)
{
    const char *filename = nullptr;
    ssize_t nt_chunk = 0;
    ssize_t noise_source_align = 0;

    if (!PyArg_ParseTuple(args, "snn", &filename, &nt_chunk, &noise_source_align))
	return NULL;

    shared_ptr<rf_pipelines::wi_stream> ret = rf_pipelines::make_chime_stream_from_filename(filename, nt_chunk, noise_source_align);
    return wi_stream_object::make(ret);    
}


static PyObject *make_chime_stream_from_filename_list(PyObject *self, PyObject *args)
{
    PyObject *arg_ptr = nullptr;
    ssize_t nt_chunk = 0;
    ssize_t noise_source_align = 0;

    if (!PyArg_ParseTuple(args, "Onn", &arg_ptr, &nt_chunk, &noise_source_align))
	return NULL;
    
    object arg(arg_ptr, false);

    if (PyString_Check(arg_ptr))
	throw runtime_error("rf_pipelines: expected argument to make_chime_stream_from_filename_list() to be a list/iterator of strings, not a single string");
    PyObject *iter_ptr = PyObject_GetIter(arg_ptr);
    if (!iter_ptr)
	throw runtime_error("rf_pipelines: expected argument to make_chime_stream_from_filename_list() to be a list/iterator of strings");

    object iter(iter_ptr, false);
    vector<object> item_list;
    vector<string> filename_list;
    
    for (;;) {
	PyObject *item_ptr = PyIter_Next(iter_ptr);
	if (!item_ptr)
	    break;

	item_list.push_back(object(item_ptr, false));

	char *s = PyString_AsString(item_ptr);
	if (!s)
	    throw runtime_error("rf_pipelines: expected argument to make_chime_stream_from_filename_list() to be a list/iterator of strings");

	filename_list.push_back(string(s));
    }

    if (PyErr_Occurred())
	throw python_exception();

    shared_ptr<rf_pipelines::wi_stream> ret = rf_pipelines::make_chime_stream_from_filename_list(filename_list, nt_chunk, noise_source_align);
    return wi_stream_object::make(ret);
}


static PyObject *make_gaussian_noise_stream(PyObject *self, PyObject *args)
{
    ssize_t nfreq, nt_chunk, nt_tot;
    double freq_lo_MHz, freq_hi_MHz, dt_sample, sample_rms;

    if (!PyArg_ParseTuple(args, "nnddddn", &nfreq, &nt_tot, &freq_lo_MHz, &freq_hi_MHz, &dt_sample, &sample_rms, &nt_chunk))
	return NULL;

    shared_ptr<rf_pipelines::wi_stream> ret = rf_pipelines::make_gaussian_noise_stream(nfreq, nt_tot, freq_lo_MHz, freq_hi_MHz, dt_sample, sample_rms, nt_chunk);
    return wi_stream_object::make(ret);
}


static PyObject *make_simple_detrender(PyObject *self, PyObject *args)
{
    ssize_t nt_detrend = 0;
    if (!PyArg_ParseTuple(args, "n", &nt_detrend))
	return NULL;
    
    shared_ptr<rf_pipelines::wi_transform> ret = rf_pipelines::make_simple_detrender(nt_detrend);
    return wi_transform_object::make(ret);
}


static PyObject *make_chime_file_writer(PyObject *self, PyObject *args)
{
    char *filename = nullptr;
    int clobber = false;   // "int" (not "bool") is deliberate here
    int bitshuffle = 2;
    int nt_chunk = 0;

    if (!PyArg_ParseTuple(args, "siii", &filename, &clobber, &bitshuffle, &nt_chunk))
	return NULL;

    shared_ptr<rf_pipelines::wi_transform> ret = rf_pipelines::make_chime_file_writer(filename, clobber, bitshuffle, nt_chunk);
    return wi_transform_object::make(ret);
}


static PyObject *make_bonsai_dedisperser(PyObject *self, PyObject *args)
{
    const char *config_hdf5_filename = nullptr;
    const char *output_hdf5_filename = nullptr;
    int nt_per_file = 0;
    int ibeam = 0;
    
    // FIXME there should be a way to disable core-pinning entirely

    if (!PyArg_ParseTuple(args, "ssii", &config_hdf5_filename, &output_hdf5_filename, &nt_per_file, &ibeam))
	return NULL;

    shared_ptr<rf_pipelines::wi_transform> ret = rf_pipelines::make_bonsai_dedisperser(config_hdf5_filename, output_hdf5_filename, nt_per_file, ibeam);
    return wi_transform_object::make(ret);
}


static constexpr const char *dummy_module_method_docstring = 
    "This is a C++ function in the rf_pipelines_c module.\n"
    "For documentation, see the docstring of the similarly-named python function in the rf_pipelines module.\n";


// -------------------------------------------------------------------------------------------------


static PyMethodDef module_methods[] = {
    { "make_psrfits_stream", tc_wrap2<make_psrfits_stream>, METH_VARARGS, dummy_module_method_docstring },
    { "make_chime_stream_from_acqdir", tc_wrap2<make_chime_stream_from_acqdir>, METH_VARARGS, dummy_module_method_docstring },
    { "make_chime_stream_from_filename", tc_wrap2<make_chime_stream_from_filename>, METH_VARARGS, dummy_module_method_docstring },
    { "make_chime_stream_from_filename_list", tc_wrap2<make_chime_stream_from_filename_list>, METH_VARARGS, dummy_module_method_docstring },
    { "make_gaussian_noise_stream", tc_wrap2<make_gaussian_noise_stream>, METH_VARARGS, dummy_module_method_docstring },
    { "make_simple_detrender", tc_wrap2<make_simple_detrender>, METH_VARARGS, dummy_module_method_docstring },
    { "make_chime_file_writer", tc_wrap2<make_chime_file_writer>, METH_VARARGS, dummy_module_method_docstring },
    { "make_bonsai_dedisperser", tc_wrap2<make_bonsai_dedisperser>, METH_VARARGS, dummy_module_method_docstring },
    { NULL, NULL, 0, NULL }
};


PyMODINIT_FUNC initrf_pipelines_c(void)
{
    import_array();

    if (PyType_Ready(&wi_stream_type) < 0)
        return;
    if (PyType_Ready(&wi_transform_type) < 0)
        return;
    if (PyType_Ready(&wi_run_state_type) < 0)
        return;

    PyObject *m = Py_InitModule3("rf_pipelines_c", module_methods, "rf_pipelines_c: a C++ library containing low-level rf_pipelines code");
    if (!m)
	return;

    Py_INCREF(&wi_stream_type);
    PyModule_AddObject(m, "wi_stream", (PyObject *)&wi_stream_type);

    Py_INCREF(&wi_transform_type);
    PyModule_AddObject(m, "wi_transform", (PyObject *)&wi_transform_type);

    Py_INCREF(&wi_run_state_type);
    PyModule_AddObject(m, "wi_run_state", (PyObject *)&wi_run_state_type);
}
