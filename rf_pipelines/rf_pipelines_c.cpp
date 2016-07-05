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
    PyObject *get_pyobj()
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

    virtual void set_stream(const rf_pipelines::wi_stream &stream)
    {
	PyObject *sp = make_temporary_stream(stream);
	object s(sp, false);

	PyObject *retp = PyObject_CallMethod(this->get_pyobj(), (char *)"set_stream", (char *)"O", sp);
	object ret(retp, false);  // a convenient way to ensure Py_DECREF gets called, and throw an exception on failure

	if (s.get_refcount() > 1)
	    throw runtime_error("fatal: wi_transform.set_stream() callback kept a reference to the stream");
    }

    virtual void start_substream(int isubstream, double t0)
    {	
	PyObject *p = PyObject_CallMethod(this->get_pyobj(), (char *)"start_substream", (char *)"id", isubstream, t0);
	object ret(p, false);  // a convenient way to ensure Py_DECREF gets called, and throw an exception on failure
    }

    virtual void process_chunk(double t0, double t1, float *intensity, float *weights, ssize_t stride, float *pp_intensity, float *pp_weights, ssize_t pp_stride)
    {
	object np_intensity = array2d_to_python(nfreq, nt_chunk, intensity, stride);
	object np_weights = array2d_to_python(nfreq, nt_chunk, weights, stride);
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

    virtual void end_substream()
    {
	PyObject *p = PyObject_CallMethod(this->get_pyobj(), (char *)"end_substream", NULL);
	object ret(p, false);
    }
};


struct wi_transform_object {
    PyObject_HEAD

    shared_ptr<rf_pipelines::wi_transform> *pshared;

    // forward declarations (these guys need to come after 'wi_transform_type'
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
};


static PyGetSetDef wi_transform_getseters[] = {
    { (char *)"nfreq", 
      tc_wrap_getter<wi_transform_object::nfreq_getter>, 
      tc_wrap_setter<wi_transform_object::nfreq_setter>, 
      NULL, NULL },

    { (char *)"nt_chunk", 
      tc_wrap_getter<wi_transform_object::nt_chunk_getter>, 
      tc_wrap_setter<wi_transform_object::nt_chunk_setter>,
      NULL, NULL },

    { (char *)"nt_prepad", 
      tc_wrap_getter<wi_transform_object::nt_prepad_getter>, 
      tc_wrap_setter<wi_transform_object::nt_prepad_setter>,
      NULL, NULL },

    { (char *)"nt_postpad", 
      tc_wrap_getter<wi_transform_object::nt_postpad_getter>, 
      tc_wrap_setter<wi_transform_object::nt_postpad_setter>,
      NULL, NULL },

    { NULL, NULL, NULL, NULL, NULL }
};


static PyMethodDef wi_transform_methods[] = {
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
    "Transform base class (C++)",           /* tp_doc */
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
    virtual void start_substream(int isubstream, double t0) { }
    virtual void end_substream() { }

    virtual void set_stream(const rf_pipelines::wi_stream &stream)
    {
	this->nfreq = stream.nfreq;
    }

    virtual void process_chunk(double t0, double t1, float *intensity, float *weights, ssize_t stride, float *pp_intensity, float *pp_weights, ssize_t pp_stride)
    {
	if (PyErr_Occurred() || PyErr_CheckSignals())
	    throw python_exception();
    }
};


// -------------------------------------------------------------------------------------------------
//
// wi_stream wrapper class


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

	if (!s->pbare)
	    throw runtime_error("rf_pipelines: internal error: unexpected NULL pointer in wi_stream [should never happen]");

	return s->pbare;
    }
    
    static PyObject *run(PyObject *self, PyObject *arg)
    {
	rf_pipelines::wi_stream *stream = get_pbare(self);

	PyObject *iter = PyObject_GetIter(arg);
	if (!iter)
	    throw runtime_error("rf_pipelines: expected argument to wi_stream.run() to be a list/iterator");

	object iter_reference(iter, false);
	vector<object> item_references;

	vector<shared_ptr<rf_pipelines::wi_transform> > transform_list;
	transform_list.push_back(make_shared<exception_monitor> (stream->nt_maxwrite));

	for (;;) {
	    PyObject *item_ptr = PyIter_Next(iter);
	    if (!item_ptr)
		break;

	    item_references.push_back(object(item_ptr,false));

	    if (!wi_transform_object::isinstance(item_ptr))
		throw runtime_error("rf_pipelines: expected argument to wi_stream.run() to be a list/iterator of wi_transform objects");

	    transform_list.push_back(wi_transform_object::get_pshared(item_ptr));
	}	

	if (PyErr_Occurred())
	    throw python_exception();

	bool noisy=true;  // FIXME make this selectable from python
	stream->run(transform_list, noisy);

	Py_INCREF(Py_None);
	return Py_None;
    }

    // Properties

    static PyObject *nfreq_getter(PyObject *self, void *closure)
    {
	return Py_BuildValue("i", get_pbare(self)->nfreq);
    }

    static PyObject *nt_maxwrite_getter(PyObject *self, void *closure)
    {
	return Py_BuildValue("i", get_pbare(self)->nt_maxwrite);
    }

    static PyObject *freq_lo_MHz_getter(PyObject *self, void *closure)
    {
	return Py_BuildValue("d", get_pbare(self)->freq_lo_MHz);
    }

    static PyObject *freq_hi_MHz_getter(PyObject *self, void *closure)
    {
	return Py_BuildValue("d", get_pbare(self)->freq_hi_MHz);
    }

    static PyObject *dt_sample_getter(PyObject *self, void *closure)
    {
	return Py_BuildValue("d", get_pbare(self)->dt_sample);
    }
};


static PyMethodDef wi_stream_methods[] = {
    { "run", tc_wrap2<wi_stream_object::run>, METH_O, NULL },
    { NULL, NULL, 0, NULL }
};

static PyGetSetDef wi_stream_getseters[] = {
    { (char *)"nfreq", tc_wrap_getter<wi_stream_object::nfreq_getter>, NULL, NULL, NULL },
    { (char *)"nt_maxwrite", tc_wrap_getter<wi_stream_object::nt_maxwrite_getter>, NULL, NULL, NULL },
    { (char *)"freq_lo_MHz", tc_wrap_getter<wi_stream_object::freq_lo_MHz_getter>, NULL, NULL, NULL },
    { (char *)"freq_hi_MHz", tc_wrap_getter<wi_stream_object::freq_hi_MHz_getter>, NULL, NULL, NULL },
    { (char *)"dt_sample", tc_wrap_getter<wi_stream_object::dt_sample_getter>, NULL, (char *)"sample length in seconds", NULL },
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
    "Stream base class (C++)",           /* tp_doc */
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

    if (!PyArg_ParseTuple(args, "sn", &filename, &nt_chunk))
	return NULL;

    shared_ptr<rf_pipelines::wi_stream> ret = rf_pipelines::make_chime_stream_from_acqdir(filename, nt_chunk);
    return wi_stream_object::make(ret);    
}


static PyObject *make_chime_stream_from_filename(PyObject *self, PyObject *args)
{
    const char *filename = nullptr;
    ssize_t nt_chunk = 0;

    if (!PyArg_ParseTuple(args, "sn", &filename, &nt_chunk))
	return NULL;

    shared_ptr<rf_pipelines::wi_stream> ret = rf_pipelines::make_chime_stream_from_filename(filename, nt_chunk);
    return wi_stream_object::make(ret);    
}


static PyObject *make_chime_stream_from_filename_list(PyObject *self, PyObject *args)
{
    PyObject *arg_ptr = nullptr;
    ssize_t nt_chunk = 0;

    if (!PyArg_ParseTuple(args, "On", &arg_ptr, &nt_chunk))
	return NULL;
    
    object arg(arg_ptr, false);
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

    shared_ptr<rf_pipelines::wi_stream> ret = rf_pipelines::make_chime_stream_from_filename_list(filename_list, nt_chunk);
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
    ssize_t nt_chunk = 0;
    if (!PyArg_ParseTuple(args, "n", &nt_chunk))
	return NULL;
    
    shared_ptr<rf_pipelines::wi_transform> ret = rf_pipelines::make_simple_detrender(nt_chunk);
    return wi_transform_object::make(ret);
}


static PyObject *make_bonsai_dedisperser(PyObject *self, PyObject *args)
{
    const char *config_hdf5_filename = nullptr;
    const char *output_hdf5_filename = nullptr;
    int ibeam = 0;
    
    // FIXME there should be a way to disable core-pinning entirely

    if (!PyArg_ParseTuple(args, "ssi", &config_hdf5_filename, &output_hdf5_filename, &ibeam))
	return NULL;

    shared_ptr<rf_pipelines::wi_transform> ret = rf_pipelines::make_bonsai_dedisperser(config_hdf5_filename, output_hdf5_filename, ibeam);
    return wi_transform_object::make(ret);
}


// -------------------------------------------------------------------------------------------------


static PyMethodDef module_methods[] = {
    { "make_psrfits_stream", tc_wrap2<make_psrfits_stream>, METH_VARARGS, "Python interface to C++ routine" },
    { "make_chime_stream_from_acqdir", tc_wrap2<make_chime_stream_from_acqdir>, METH_VARARGS, "Python interface to C++ routine" },
    { "make_chime_stream_from_filename", tc_wrap2<make_chime_stream_from_filename>, METH_VARARGS, "Python interface to C++ routine" },
    { "make_chime_stream_from_filename_list", tc_wrap2<make_chime_stream_from_filename_list>, METH_VARARGS, "Python interface to C++ routine" },
    { "make_gaussian_noise_stream", tc_wrap2<make_gaussian_noise_stream>, METH_VARARGS, "Python interface to C++ routine" },
    { "make_simple_detrender", tc_wrap2<make_simple_detrender>, METH_VARARGS, "Python interface to C++ routine" },
    { "make_bonsai_dedisperser", tc_wrap2<make_bonsai_dedisperser>, METH_VARARGS, "Python interface to C++ routine" },
    { NULL, NULL, 0, NULL }
};


PyMODINIT_FUNC initrf_pipelines_c(void)
{
    import_array();

    if (PyType_Ready(&wi_stream_type) < 0)
        return;
    if (PyType_Ready(&wi_transform_type) < 0)
        return;

    PyObject *m = Py_InitModule3("rf_pipelines_c", module_methods, "Python interface to C++ library");
    if (!m)
	return;

    Py_INCREF(&wi_stream_type);
    PyModule_AddObject(m, "wi_stream", (PyObject *)&wi_stream_type);

    Py_INCREF(&wi_transform_type);
    PyModule_AddObject(m, "wi_transform", (PyObject *)&wi_transform_type);
}
