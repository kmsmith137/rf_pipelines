// Note: I haven't systematically documented the C++ interface to rf_pipelines,
// so the level of documentation will be hit-or-miss.  Also please note that the
// python-wrapping in rf_pipelines_c.cpp is kind of a mess which I hope to improve
// soon.  In the meantime if you want to python-wrap a C++ class, just email me
// and I'll help navigate the mess!

#include "python_extension_helpers.hpp"
#include "rf_pipelines_internals.hpp"

using namespace std;

static PyObject *make_temporary_stream(const rf_pipelines::wi_stream &s);


// -------------------------------------------------------------------------------------------------
//
// arr_2d_helper: a helper class for receiving an array from python
// arr_wi_helper: a helper class for receiving an (intensity, weights) array pair from python


struct arr_2d_helper {
    object a_ref;
    bool writeback = false;

    float *data = nullptr;
    int nfreq = 0;
    int nt = 0;
    int stride = 0;


    // Helper function for constructor
    inline bool try_array(PyArrayObject *a, bool writeback)
    {
	if (PyArray_TYPE(a) != NPY_FLOAT)
	    return false;
	if (PyArray_ITEMSIZE(a) != sizeof(float))
	    return false;
	if (PyArray_NDIM(a) != 2)
	    return false;
	if (PyArray_STRIDE(a,1) != sizeof(float))
	    return false;
	if (PyArray_STRIDE(a,0) % sizeof(float))
	    return false;
	if (writeback && ((PyArray_FLAGS(a) & NPY_ARRAY_WRITEABLE) != NPY_ARRAY_WRITEABLE))
	    return false;

	this->data = reinterpret_cast<float *> (PyArray_DATA(a));
	this->nfreq = PyArray_DIM(a, 0);
	this->nt = PyArray_DIM(a, 1);
	this->stride = PyArray_STRIDE(a,0) / sizeof(float);

	return true;
    }

    
    void set_contiguous(PyObject *obj)
    {
	// The NPY_ARRAY_FORCECAST flag allows conversion double -> float
	int requirements = NPY_ARRAY_C_CONTIGUOUS | NPY_ARRAY_ALIGNED | NPY_ARRAY_ENSUREARRAY | NPY_ARRAY_FORCECAST;

	if (writeback)
	    requirements |= NPY_ARRAY_UPDATEIFCOPY;

	// FIXME I'm a little hazy on whether PyArray_DescrFromType() increments the refcount.
	PyObject *a0 = PyArray_FromAny(obj, PyArray_DescrFromType(NPY_FLOAT), 2, 2, requirements, NULL);
	this->a_ref = object(a0, false);   // manages refcount, throws exception on NULL

	PyArrayObject *a = (PyArrayObject *) a0;

	// If either of these exception is thrown, then my understsanding of the python array API is incomplete!
	if (!try_array(a, writeback))
	    throw runtime_error("rf_pipelines internal error: try_array() failed in arr_2d_helper::set_contiguous()");
	if (stride != nt)
	    throw runtime_error("rf_pipelines internal error: unexpected stride in arr_2d_helper::set_contiguous()");
    }


    void set_contiguous()
    {
	assert(a_ref.ptr);
	set_contiguous(a_ref.ptr);
    }

    arr_2d_helper(PyObject *obj, bool writeback_)
	: writeback(writeback_)
    {
	// FIXME do I want PyArray_Check() or PyArray_CheckExact() here?

	if (PyArray_CheckExact(obj) && try_array((PyArrayObject *) obj, writeback)) {
	    // fast track: array works "as is", no need to convert or cast anything
	    this->a_ref = object(obj, true);   // increment_refcount = true
	    return;
	}

	set_contiguous(obj);
    }
};


struct arr_wi_helper {
    arr_2d_helper intensity;
    arr_2d_helper weights;

    // Same for both arrays.
    int nfreq;
    int nt;
    int stride;
    
    arr_wi_helper(PyObject *intensity_obj, PyObject *weights_obj, bool intensity_writeback, bool weights_writeback) :
	intensity(intensity_obj, intensity_writeback),
	weights(weights_obj, weights_writeback)
    {
	if ((intensity.nfreq != weights.nfreq) || (intensity.nt != weights.nt))
	    throw runtime_error("intensity and weights arrays have mismatched dimensions (nfreq, nt)");

	this->nfreq = intensity.nfreq;
	this->nt = intensity.nt;

	if (intensity.stride != weights.stride) {
	    if (intensity.stride != nt)
		intensity.set_contiguous();
	    if (weights.stride != nt)
		weights.set_contiguous();
	}

	this->stride = intensity.stride;
    }
};


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


    static PyObject *name_getter(PyObject *self, void *closure)
    {
	return Py_BuildValue("s", get_pbare(self)->name.c_str());   // Note: Py_BuildValue() copies the string
    }

    static int name_setter(PyObject *self, PyObject *value, void *closure)
    {
	get_pbare(self)->name = string_from_python(value);
	return 0;
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

    static constexpr const char *name_docstring =
	"Transform name (string).  Ends up in rf_pipelines.json.  Should be initialized in the transform constructor.";

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
    { (char *)"name",
      tc_wrap_getter<wi_transform_object::name_getter>, 
      tc_wrap_setter<wi_transform_object::name_setter>, 
      (char *)wi_transform_object::name_docstring, NULL },
      
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
    exception_monitor()
    {
	this->name = "exception_monitor";
    }

    virtual ~exception_monitor() { }

    virtual void set_stream(const rf_pipelines::wi_stream &stream) override
    {
	this->nfreq = stream.nfreq;
	this->nt_chunk = stream.nt_maxwrite;
    }

    virtual void start_substream(int isubstream, double t0) override { }

    virtual void process_chunk(double t0, double t1, float *intensity, float *weights, ssize_t stride, float *pp_intensity, float *pp_weights, ssize_t pp_stride) override
    {
	if (PyErr_Occurred() || PyErr_CheckSignals())
	    throw python_exception();
    }

    virtual void end_substream() override { }
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

    virtual void stream_start()
    {
	PyObject *p = PyObject_CallMethod(this->get_pyobj(), (char *)"stream_start", NULL);
	object ret(p, false);
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
	object default_outdir(Py_BuildValue("s","."), false);

	rf_pipelines::wi_stream *stream = get_pbare(self);
	PyObject *transforms_obj = Py_None;
	PyObject *outdir_obj = default_outdir.ptr;
	int noisy = 1;
	int clobber = 1;
	int return_json = 0;

	if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|Oiii", (char **)kwlist, &transforms_obj, &outdir_obj, &noisy, &clobber, &return_json))
	    return NULL;

	string outdir;
	if (outdir_obj != Py_None)
	    outdir = string_from_python(outdir_obj);

	PyObject *transforms_iter = PyObject_GetIter(transforms_obj);
	if (!transforms_iter)
	    throw runtime_error("rf_pipelines: expected 'transforms' argument to wi_stream.run() to be a list/iterator");

	object iter_reference(transforms_iter, false);
	vector<object> item_references;

	vector<shared_ptr<rf_pipelines::wi_transform> > transform_list;
	transform_list.push_back(make_shared<exception_monitor> ());

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

    static constexpr const char *run_docstring =
	"run(self, transform_list, outdir='.', noisy=True, clobber=True, return_json=False)\n"
	"\n"
	"This function is called to run an rf_pipeline.  Arguments:\n"
        "\n"
	"  - 'transform_list' is a list (or generator) of objects of type wi_transform (including\n"
	"     its subclass py_wi_transform).\n"
	"\n"
	"  - 'outdir' is the rf_pipelines output directory, where the rf_pipelines json file will\n"
	"     be written, in addition to other transform-specific output files such as plots\n"
	"\n"
	"  -  If 'outdir' is None or an empty string, then the json file will not be written,\n"
	"     and any transform which tries to write an output file (such as a plotter_transform)\n"
	"     will throw an exception.\n"
	"\n"
	"  -  If 'clobber' is False, then an exception will be thrown if the pipeline tries to\n"
	"     overwrite an old rf_pipelines.json file.\n"
	"\n"
	"  -  If 'return_json' is True, then the return value from run() will be the rf_pipelines\n"
	"     json output (i.e. same data which is written to rf_pipelines.json)\n"
	"\n"
	"     A kludge: eventually, the run() return value will be a json object, but for now it returns\n"
	"     the string representation, which can be converted to a json object by calling json.loads().\n";

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
    { "run", (PyCFunction) tc_wrap3<wi_stream_object::run>, METH_VARARGS | METH_KEYWORDS, wi_stream_object::run_docstring },
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


static PyObject *make_chime_network_stream(PyObject *self, PyObject *args)
{
    int udp_port, beam_id;
    
    if (!PyArg_ParseTuple(args, "ii", &udp_port, &beam_id))
	return NULL;

    shared_ptr<rf_pipelines::wi_stream> ret = rf_pipelines::make_chime_network_stream(udp_port, beam_id);
    return wi_stream_object::make(ret);
}


static PyObject *make_gaussian_noise_stream(PyObject *self, PyObject *args)
{
    ssize_t nfreq, nt_chunk, nt_tot;
    double freq_lo_MHz, freq_hi_MHz, dt_sample, sample_rms;
    int randomize_weights;  // will be converted to boolean

    if (!PyArg_ParseTuple(args, "nnddddni", &nfreq, &nt_tot, &freq_lo_MHz, &freq_hi_MHz, &dt_sample, &sample_rms, &nt_chunk, &randomize_weights))
	return NULL;

    shared_ptr<rf_pipelines::wi_stream> ret = rf_pipelines::make_gaussian_noise_stream(nfreq, nt_tot, freq_lo_MHz, freq_hi_MHz, dt_sample, sample_rms, nt_chunk, (randomize_weights != 0));
    return wi_stream_object::make(ret);
}


static PyObject *make_chime_packetizer(PyObject *self, PyObject *args)
{
    const char *dstname = nullptr;
    int nfreq_per_packet;
    int nt_per_chunk;
    int nt_per_packet;
    double wt_cutoff;
    double target_gbps;

    if (!PyArg_ParseTuple(args, "siiidd", &dstname, &nfreq_per_packet, &nt_per_chunk, &nt_per_packet, &wt_cutoff, &target_gbps))
	return NULL;

    shared_ptr<rf_pipelines::wi_transform> ret = rf_pipelines::make_chime_packetizer(dstname, nfreq_per_packet, nt_per_chunk, nt_per_packet, wt_cutoff, target_gbps);
    return wi_transform_object::make(ret);
}


static rf_pipelines::axis_type axis_type_from_python(const char *function_name, PyObject *obj)
{
    if (obj == Py_None)
	return rf_pipelines::AXIS_NONE;

    if (!PyInt_Check(obj))
	throw runtime_error(string(function_name) + ": bad 'axis' parameter");

    ssize_t ret = PyInt_AsSsize_t(obj);

    if (ret == 0)
	return rf_pipelines::AXIS_FREQ;
    if (ret == 1)
	return rf_pipelines::AXIS_TIME;

    if ((ret == -1) && PyErr_Occurred())
	throw python_exception();

    throw runtime_error(string(function_name) + ": bad 'axis' parameter");    
}


static PyObject *make_polynomial_detrender(PyObject *self, PyObject *args, PyObject *kwds)
{
    static const char *kwlist[] = { "nt_chunk", "axis", "polydeg", "epsilon", NULL };

    int nt_chunk = 0;
    PyObject *axis_ptr = Py_None;
    int polydeg = 0;
    double epsilon = 1.0e-2;   // meaningful default value

    // Note: the object pointers will be borrowed references
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "iOi|d", (char **)kwlist, &nt_chunk, &axis_ptr, &polydeg, &epsilon))
	return NULL;

    rf_pipelines::axis_type axis = axis_type_from_python("make_intensity_clipper()", axis_ptr);

    shared_ptr<rf_pipelines::wi_transform> ret = rf_pipelines::make_polynomial_detrender(nt_chunk, axis, polydeg, epsilon);
    return wi_transform_object::make(ret);
}


static PyObject *make_intensity_clipper(PyObject *self, PyObject *args, PyObject *kwds)
{
    static const char *kwlist[] = { "nt_chunk", "axis", "sigma", "niter", "iter_sigma", "Df", "Dt", "two_pass", NULL }; 

    int nt_chunk = 0;
    PyObject *axis_ptr = Py_None;
    double sigma = 0.0;
    int niter = 1;             // meaningful default value
    double iter_sigma = 0.0;   // meaningful default value
    int Df = 1;                // meaningful default value
    int Dt = 1;                // meaningful default value
    int two_pass = 0;          // meaningful default value

    // Note: the object pointers will be borrowed references
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "iOd|idiii", (char **)kwlist, &nt_chunk, &axis_ptr, &sigma, &niter, &iter_sigma, &Df, &Dt, &two_pass))
	return NULL;

    rf_pipelines::axis_type axis = axis_type_from_python("make_intensity_clipper()", axis_ptr);

    shared_ptr<rf_pipelines::wi_transform> ret = rf_pipelines::make_intensity_clipper(nt_chunk, axis, sigma, niter, iter_sigma, Df, Dt, two_pass);
    return wi_transform_object::make(ret);
}


static PyObject *make_std_dev_clipper(PyObject *self, PyObject *args, PyObject *kwds)
{
    static const char *kwlist[] = { "nt_chunk", "axis", "sigma", "Df", "Dt", "two_pass", NULL };

    int nt_chunk = 0;
    PyObject *axis_ptr = Py_None;
    double sigma = 0.0;
    int Df = 1;        // meaningful default value
    int Dt = 1;        // meaningful default value
    int two_pass = 0;  // meaningful default value

    // Note: the object pointers will be borrowed references
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "iOd|iii", (char **)kwlist, &nt_chunk, &axis_ptr, &sigma, &Df, &Dt, &two_pass))
	return NULL;

    rf_pipelines::axis_type axis = axis_type_from_python("make_std_dev_clipper()", axis_ptr);

    shared_ptr<rf_pipelines::wi_transform> ret = rf_pipelines::make_std_dev_clipper(nt_chunk, axis, sigma, Df, Dt, two_pass);
    return wi_transform_object::make(ret);
}


static PyObject *apply_polynomial_detrender(PyObject *self, PyObject *args, PyObject *kwds)
{
    static const char *kwlist[] = { "intensity", "weights", "axis", "polydeg", "epsilon", NULL };

    PyObject *intensity_obj = Py_None;
    PyObject *weights_obj = Py_None;
    PyObject *axis_ptr = Py_None;
    int polydeg = -1;
    double epsilon = 1.0e-2;   // meaningful default value

    // Note: the object pointers will be borrowed references
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "OOOi|d", (char **)kwlist, &intensity_obj, &weights_obj, &axis_ptr, &polydeg, &epsilon))
	return NULL;

    // (intensity_writeback, weights_writeback) = (true, true)
    // Note that the weights can be updated if the fit is poorly conditioned.
    arr_wi_helper wi(intensity_obj, weights_obj, true, true); 

    rf_pipelines::axis_type axis = axis_type_from_python("apply_polynomial_detrender", axis_ptr);

    rf_pipelines::apply_polynomial_detrender(wi.intensity.data, wi.weights.data, wi.nfreq, wi.nt, wi.stride, axis, polydeg, epsilon);

    Py_INCREF(Py_None);
    return Py_None;
}


static PyObject *apply_intensity_clipper(PyObject *self, PyObject *args, PyObject *kwds)
{
    static const char *kwlist[] = { "intensity", "weights", "axis", "sigma", "niter", "iter_sigma", "Df", "Dt", "two_pass", NULL };

    PyObject *intensity_obj = Py_None;
    PyObject *weights_obj = Py_None;
    PyObject *axis_ptr = Py_None;
    double sigma = 0.0;
    int niter = 1;             // meaningful default value
    double iter_sigma = 0.0;   // meaningful default value
    int Df = 1;                // meaningful default value
    int Dt = 1;                // meaningful default value
    int two_pass = 0;          // meaningful default value

    // Note: the object pointers will be borrowed references
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "OOOd|idiii", (char **)kwlist, &intensity_obj, &weights_obj, &axis_ptr, &sigma, &niter, &iter_sigma, &Df, &Dt, &two_pass))
	return NULL;

    arr_wi_helper wi(intensity_obj, weights_obj, false, true);   // (intensity_writeback, weights_writeback) = (false, true)

    rf_pipelines::axis_type axis = axis_type_from_python("apply_intensity_clipper", axis_ptr);

    rf_pipelines::apply_intensity_clipper(wi.intensity.data, wi.weights.data, wi.nfreq, wi.nt, wi.stride, axis, sigma, niter, iter_sigma, Df, Dt, two_pass);

    Py_INCREF(Py_None);
    return Py_None;    
}


static PyObject *apply_std_dev_clipper(PyObject *self, PyObject *args, PyObject *kwds)
{
    static const char *kwlist[] = { "intensity", "weights", "axis", "sigma", "Df", "Dt", "two_pass", NULL };

    PyObject *intensity_obj = Py_None;
    PyObject *weights_obj = Py_None;
    PyObject *axis_ptr = Py_None;
    double sigma = 0.0;
    int Df = 1;   // meaningful default value
    int Dt = 1;   // meaningful default value
    int two_pass = 0;  // meaningful default valuex

    // Note: the object pointers will be borrowed references
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "OOOd|iii", (char **)kwlist, &intensity_obj, &weights_obj, &axis_ptr, &sigma, &Df, &Dt, &two_pass))
	return NULL;

    arr_wi_helper wi(intensity_obj, weights_obj, false, true);   // (intensity_writeback, weights_writeback) = (false, true)

    rf_pipelines::axis_type axis = axis_type_from_python("apply_std_dev_clipper", axis_ptr);

    rf_pipelines::apply_std_dev_clipper(wi.intensity.data, wi.weights.data, wi.nfreq, wi.nt, wi.stride, axis, sigma, Df, Dt, two_pass);

    Py_INCREF(Py_None);
    return Py_None;
}


static PyObject *wi_downsample(PyObject *self, PyObject *args, PyObject *kwds)
{
    static const char *kwlist[] = { "intensity", "weights", "Df", "Dt", NULL };

    PyObject *intensity_obj = Py_None;
    PyObject *weights_obj = Py_None;
    int Df = 0;
    int Dt = 0;
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "OOii|", (char **)kwlist, &intensity_obj, &weights_obj, &Df, &Dt))
	return NULL;

    // Argument-checking is done systematically in wi_downsample(), but we need this check up-front
    // to make sure that we don't divide by zero when computing out_shape.

    if ((Df <= 0) || (Dt <= 0))
	throw runtime_error("wi_downsample(): (Df,Dt)=(" + to_string(Df) + "," + to_string(Dt) + ") is invalid");

    // (intensity_writeback, weights_writeback) = (false, false)
    arr_wi_helper wi(intensity_obj, weights_obj, false, false);

    npy_intp out_shape[2] = { (wi.nfreq / Df), (wi.nt / Dt) };

    PyObject *ds_iptr = PyArray_SimpleNew(2, out_shape, NPY_FLOAT);
    object ds_iobj(ds_iptr, false);   // manage refcount

    PyObject *ds_wptr = PyArray_SimpleNew(2, out_shape, NPY_FLOAT);
    object ds_wobj(ds_wptr, false);   // manage refcount
    
    rf_pipelines::wi_downsample((float *) PyArray_DATA((PyArrayObject *) ds_iptr),   // out_intensity
				(float *) PyArray_DATA((PyArrayObject *) ds_wptr),   // out_weights
				wi.nt / Dt, wi.intensity.data, wi.weights.data,
				wi.nfreq, wi.nt, wi.stride, Df, Dt);

    PyObject *ret = PyTuple_Pack(2, ds_iptr, ds_wptr);

    // Note: PyTuple_Pack() increments (ds_iptr, ds_wptr) refcounts when it creates the tuple.
    // When this routine exits, the (ds_iobj, ds_wobj) destructors will decrement the refcounts.
    // It follows that we don't need calls to either Py_INCREF() or Py_DECREF() here.

    return ret;
}


static PyObject *weighted_mean_and_rms(PyObject *self, PyObject *args, PyObject *kwds)
{
    static const char *kwlist[] = { "intensity", "weights", "niter", "sigma", "two_pass", NULL };

    PyObject *intensity_obj = Py_None;
    PyObject *weights_obj = Py_None;
    int niter = 1;       // meaningful default value
    double sigma = 3.0;  // meaningful default value
    int two_pass = 0;    // meaningful default value

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "OO|idi", (char **)kwlist, &intensity_obj, &weights_obj, &niter, &sigma, &two_pass))
	return NULL;

    // (intensity_writeback, weights_writeback) = (false, false)
    arr_wi_helper wi(intensity_obj, weights_obj, false, false);

    float mean, rms;
    rf_pipelines::weighted_mean_and_rms(mean, rms, wi.intensity.data, wi.weights.data, wi.nfreq, wi.nt, wi.stride, niter, sigma, two_pass);

    return Py_BuildValue("(dd)", mean, rms);
}


static PyObject *_wrms_hack_for_testing1(PyObject *self, PyObject *args, PyObject *kwds)
{
    static const char *kwlist[] = { "intensity", "weights", "niter", "sigma", "two_pass", NULL };

    PyObject *intensity_obj = Py_None;
    PyObject *weights_obj = Py_None;
    int niter = 1;       // meaningful default value
    double sigma = 3.0;  // meaningful default value
    int two_pass = 0;    // meaningful default value

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "OO|idi", (char **)kwlist, &intensity_obj, &weights_obj, &niter, &sigma, &two_pass))
	return NULL;

    // (intensity_writeback, weights_writeback) = (false, false)
    arr_wi_helper wi(intensity_obj, weights_obj, false, false);

    vector<float> mean_hint;
    rf_pipelines::_wrms_hack_for_testing1(mean_hint, wi.intensity.data, wi.weights.data, wi.nfreq, wi.nt, wi.stride, niter, sigma, two_pass);

    npy_intp nelts = mean_hint.size();
    PyObject *ret = PyArray_SimpleNew(1, &nelts, NPY_FLOAT);

    memcpy((float *) PyArray_DATA((PyArrayObject *) ret), &mean_hint[0], nelts * sizeof(float));
    return ret;
}


static PyObject *_wrms_hack_for_testing2(PyObject *self, PyObject *args, PyObject *kwds)
{
    static const char *kwlist[] = { "intensity", "weights", "mean_hint", NULL };    

    PyObject *intensity_obj = Py_None;
    PyObject *weights_obj = Py_None;
    PyObject *mean_hint_obj = Py_None;

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "OOO", (char **)kwlist, &intensity_obj, &weights_obj, &mean_hint_obj))
	return NULL;

    // (intensity_writeback, weights_writeback) = (false, false)
    arr_wi_helper wi(intensity_obj, weights_obj, false, false);

    int requirements = NPY_ARRAY_C_CONTIGUOUS | NPY_ARRAY_ALIGNED | NPY_ARRAY_ENSUREARRAY | NPY_ARRAY_FORCECAST;
    PyArrayObject *mean_hint_arr = (PyArrayObject *) PyArray_FromAny(mean_hint_obj, PyArray_DescrFromType(NPY_FLOAT), 1, 1, requirements, NULL);
    
    assert(PyArray_NDIM(mean_hint_arr) == 1);
    assert(PyArray_TYPE(mean_hint_arr) == NPY_FLOAT);
    assert(PyArray_STRIDE(mean_hint_arr,0) == sizeof(float));

    int n = PyArray_DIM(mean_hint_arr, 0);
    vector<float> mean_hint(n);
    memcpy(&mean_hint[0], (float *) PyArray_DATA(mean_hint_arr), n * sizeof(float));

    float mean, rms;
    rf_pipelines::_wrms_hack_for_testing2(mean, rms, wi.intensity.data, wi.weights.data, wi.nfreq, wi.nt, wi.stride, mean_hint);

    return Py_BuildValue("(dd)", mean, rms);
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
    const char *trigger_hdf5_filename = nullptr;
    const char *trigger_plot_stem = nullptr;
    int nt_per_file = 0;
    int ibeam = 0;
    
    if (!PyArg_ParseTuple(args, "sssii", &config_hdf5_filename, &trigger_hdf5_filename, &trigger_plot_stem, &nt_per_file, &ibeam))
	return NULL;

    shared_ptr<rf_pipelines::wi_transform> ret = rf_pipelines::make_bonsai_dedisperser(config_hdf5_filename, trigger_hdf5_filename, trigger_plot_stem, nt_per_file, ibeam);
    return wi_transform_object::make(ret);
}


// extern std::shared_ptr<wi_transform> make_badchannel_mask(const std::string &maskpath, int nt_chunk=1024);
static PyObject *make_badchannel_mask(PyObject *self, PyObject *args)
{
    const char *maskpath = nullptr;
    int nt_chunk = 0;
    
    if (!PyArg_ParseTuple(args, "si", &maskpath, &nt_chunk))
	return NULL;

    shared_ptr<rf_pipelines::wi_transform> ret = rf_pipelines::make_badchannel_mask(maskpath, nt_chunk);
    return wi_transform_object::make(ret);
}


// FIXME improve?
static constexpr const char *dummy_module_method_docstring = 
    "This is a C++ function in the rf_pipelines_c module.\n"
    "For documentation, see the docstring of the similarly-named python function in the rf_pipelines module.\n";


static constexpr const char *make_gaussian_noise_stream_docstring =
    "make_gaussian_noise_stream(nfreq, nt_tot, freq_lo_MHz, freq_hi_MHz, dt_sample, sample_rms, nt_chunk, randomize_weights)\n"
    "\n"
    "nfreq               Number of frequency channels\n"
    "nt_tot              Total number of time samples written before stream ends.\n"
    "freq_lo_MHz         Lowest frequency in band (e.g. 400 for CHIME)\n"
    "freq_hi_MHz         Highest frequency in band (e.g. 800 for CHIME)\n"
    "dt_sample           Length of a time sample in seconds\n"
    "nt_chunk            Stream block size (if zero, will default to a reasonable value)\n"
    "randomize_weights   If true, weights will be uniform random numbers (if false, all weights will be 1.0)\n";


static constexpr const char *make_polynomial_detrender_docstring =
    "make_polynomial_detrender(nt_chunk, axis, polydeg, epsilon = 1.0e-2)\n"
    "\n"
    "Detrends along the specified axis by subtracting a best-fit polynomial.\n"
    "axis=0 means 'detrend in time', axis=1 means 'detrend in frequency'.\n"
    "\n"
    "If the fit is poorly conditioned then the entire frequency channel (FIXME or the entire time sample when axis=1 ?? ) will be masked\n"
    "(by setting its weights to zero).  The threshold is controlled by the parameter\n"
    "'epsilon'.  I think that 1.0e-2 is a reasonable default here, but haven't\n"
    "experimented systematically.\n";


static constexpr const char *make_intensity_clipper_docstring =
    "make_intensity_clipper(nt_chunk, axis, sigma, niter=1, iter_sigma=0, Df=1, Dt=1)\n"
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
    "If the 'two_pass' flag is set, a more numerically stable but slightly slower algorithm will be used.\n";


static constexpr const char *make_std_dev_clipper_docstring =
    "make_std_dev_clipper(nt_chunk, axis, sigma, Df=1, Dt=1, two_pass=False)\n"
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
    "If the 'two_pass' flag is set, then a more numerically stable but slightly slower algorithm will be used.\n";


static constexpr const char *apply_polynomial_detrender_docstring =
    "apply_polynomial_detrender(intensity, weights, axis, polydeg, epsilon=1.0e-2)\n"
    "\n"
    "Detrends along the specified axis by subtracting a best-fit polynomial.\n"
    "axis=0 means 'detrend in time', axis=1 means 'detrend in frequency'.\n"
    "\n"
    "If the fit is poorly conditioned then the entire frequency channel will be masked\n"
    "(by setting its weights to zero).  The threshold is controlled by the parameter\n"
    "'epsilon'.  I think that 1.0e-2 is a reasonable default here, but haven't\n"
    "experimented systematically.\n";


static constexpr const char *apply_intensity_clipper_docstring =
    "apply_intensity_clipper(intensity, weights, axis, sigma, niter=1, iter_sigma=0.0, Df=1, Dt=1)\n"
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
    "to 'sigma', but the two thresholds need not be the same.\n";


static constexpr const char *apply_std_dev_clipper_docstring =
    "apply_std_dev_clipper(intensity, weights, axis, sigma, Df=1, Dt=1, two_pass=False)\n"
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
    "If the 'two_pass' flag is set, then a more numerically stable but slightly slower algorithm will be used.\n";


static constexpr const char *wi_downsample_docstring =
    "wi_downsample(intensity, weights, Df, Dt)\n"
    "\n"
    "Downsamples a weighted intensity array, and returns a new pair (intensity, weights).\n"
    "The downsampling factors (Df,Dt) must be powers of two.\n"
    "\n"
    "Note that the normalization of the downsampled 'weights' array differs\n"
    "(by a factor of Df*Dt) from the python version of wi_downsample().\n";


static constexpr const char *weighted_mean_and_rms_docstring =
    "weighted_mean_and_rms(intensity, weights, niter=1, sigma=3.0)\n"
    "\n"
    "Computes weighted mean/rms of a 2D intensity array.\n"
    "If the 'niter' argument is >1, then the calculation will be iterated, clipping\n"
    "outlier samples which differ from the mean by the specified number of \"sigmas\".\n"
    "If the 'two_pass' flag is set, a more numerically stable but slightly slower algorithm will be used.\n";


static constexpr const char *wrms_hack_for_testing_docstring =
    "The \"wrms_hack_for_testing\" is explained in test-cpp-python-equivalence.py";

static constexpr const char *make_badchannel_mask_docstring = 
    "make_badchannel_mask(maskpath, nt_chunk)\n"
    "\n"
    "Some day, this factory function will return a C++ implementation of the 'badchannel_mask' class.\n"
    "Right now, it is a placeholder which throws an exception if called.\n";



// -------------------------------------------------------------------------------------------------


static PyMethodDef module_methods[] = {
    { "make_psrfits_stream", tc_wrap2<make_psrfits_stream>, METH_VARARGS, dummy_module_method_docstring },
    { "make_chime_stream_from_acqdir", tc_wrap2<make_chime_stream_from_acqdir>, METH_VARARGS, dummy_module_method_docstring },
    { "make_chime_stream_from_filename", tc_wrap2<make_chime_stream_from_filename>, METH_VARARGS, dummy_module_method_docstring },
    { "make_chime_stream_from_filename_list", tc_wrap2<make_chime_stream_from_filename_list>, METH_VARARGS, dummy_module_method_docstring },
    { "make_chime_network_stream", tc_wrap2<make_chime_network_stream>, METH_VARARGS, dummy_module_method_docstring },
    { "make_gaussian_noise_stream", tc_wrap2<make_gaussian_noise_stream>, METH_VARARGS, make_gaussian_noise_stream_docstring },
    { "make_chime_packetizer", tc_wrap2<make_chime_packetizer>, METH_VARARGS, dummy_module_method_docstring },
    { "make_polynomial_detrender", (PyCFunction) tc_wrap3<make_polynomial_detrender>, METH_VARARGS, make_polynomial_detrender_docstring },
    { "make_intensity_clipper", (PyCFunction) tc_wrap3<make_intensity_clipper>, METH_VARARGS, make_intensity_clipper_docstring },
    { "make_std_dev_clipper", (PyCFunction) tc_wrap3<make_std_dev_clipper>, METH_VARARGS, make_std_dev_clipper_docstring },
    { "make_chime_file_writer", tc_wrap2<make_chime_file_writer>, METH_VARARGS, dummy_module_method_docstring },
    { "make_bonsai_dedisperser", tc_wrap2<make_bonsai_dedisperser>, METH_VARARGS, dummy_module_method_docstring },
    { "make_badchannel_mask", tc_wrap2<make_badchannel_mask>, METH_VARARGS, make_badchannel_mask_docstring },
    { "apply_polynomial_detrender", (PyCFunction) tc_wrap3<apply_polynomial_detrender>, METH_VARARGS | METH_KEYWORDS, apply_polynomial_detrender_docstring },
    { "apply_intensity_clipper", (PyCFunction) tc_wrap3<apply_intensity_clipper>, METH_VARARGS | METH_KEYWORDS, apply_intensity_clipper_docstring },
    { "apply_std_dev_clipper", (PyCFunction) tc_wrap3<apply_std_dev_clipper>, METH_VARARGS | METH_KEYWORDS, apply_std_dev_clipper_docstring },
    { "wi_downsample", (PyCFunction) tc_wrap3<wi_downsample>, METH_VARARGS | METH_KEYWORDS, wi_downsample_docstring },
    { "weighted_mean_and_rms", (PyCFunction) tc_wrap3<weighted_mean_and_rms>, METH_VARARGS | METH_KEYWORDS, weighted_mean_and_rms_docstring },
    { "_wrms_hack_for_testing1", (PyCFunction) tc_wrap3<_wrms_hack_for_testing1>, METH_VARARGS | METH_KEYWORDS, wrms_hack_for_testing_docstring },
    { "_wrms_hack_for_testing2", (PyCFunction) tc_wrap3<_wrms_hack_for_testing2>, METH_VARARGS | METH_KEYWORDS, wrms_hack_for_testing_docstring },
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
