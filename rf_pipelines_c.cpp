#include "python_extension_helpers.hpp"
#include "rf_pipelines.hpp"

using namespace std;

static PyObject *make_dummy_stream(const rf_pipelines::wi_stream &s);


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
    static object array2d_to_python(int m, int n, const float *src, int src_stride)
    {
	if ((m <= 0) || (n <= 0) || !src)
	    throw runtime_error("rf_pipelines: array2d_to_python: internal checks failed [should never happen]");

	npy_intp dims[2] = { m, n };
	object ret(PyArray_SimpleNew(2,dims,NPY_FLOAT), false);
	PyArrayObject *arr = (PyArrayObject *) ret.ptr;
	
	if (!PyArray_Check(ret.ptr) || (PyArray_NDIM(arr) != 2) 
	    || (PyArray_DIM(arr,0) != m) || (PyArray_DIM(arr,1) != n)
	    || (PyArray_STRIDE(arr,0) != (int)(n*sizeof(float))) 
	    || (PyArray_STRIDE(arr,1) != sizeof(float)))
	    throw runtime_error("rf_pipelines: array2d_to_python: internal checks failed [should never happen]");

	float *dst = (float *) PyArray_DATA(arr);

	for (int i = 0; i < m; i++)
	    memcpy(dst + i*n, src + i*src_stride, n * sizeof(float));

	return ret;
    }

    static void array2d_from_python(const object &obj, int m, int n, float *dst, int dst_stride)
    {
	PyArrayObject *arr = (PyArrayObject *) obj.ptr;
	
	if (!PyArray_Check(obj.ptr) || (PyArray_NDIM(arr) != 2) 
	    || (PyArray_DIM(arr,0) != m) || (PyArray_DIM(arr,1) != n)
	    || (PyArray_STRIDE(arr,0) != (int)(n*sizeof(float))) 
	    || (PyArray_STRIDE(arr,1) != sizeof(float)))
	    throw runtime_error("rf_pipelines: array2d_to_python: internal checks failed [should never happen]");

	const float *src = (float *) PyArray_DATA(arr);

	for (int i = 0; i < m; i++)
	    memcpy(dst + i*dst_stride, src + i*n, n * sizeof(float));
    }

    virtual void set_stream(const rf_pipelines::wi_stream &stream)
    {
	object s(make_dummy_stream(stream), false);
	PyObject *p = PyObject_CallMethod(this->get_pyobj(), (char *)"set_stream", (char *)"O", s.ptr);
	object ret(p, false);
    }

    virtual void start_substream(double t0)
    {	
	PyObject *p = PyObject_CallMethod(this->get_pyobj(), (char *)"start_substream", (char *)"d", t0);
	object ret(p, false);
    }

    virtual void process_chunk(float *intensity, float *weights, int stride, float *pp_intensity, float *pp_weights, int pp_stride)
    {
	object np_intensity = array2d_to_python(nfreq, nt_chunk, intensity, stride);
	object np_weights = array2d_to_python(nfreq, nt_chunk, weights, stride);
	object np_pp_intensity;
	object np_pp_weights;

	if (nt_prepad > 0) {
	    np_pp_intensity = array2d_to_python(nfreq, nt_prepad, pp_intensity, pp_stride);
	    np_pp_weights = array2d_to_python(nfreq, nt_prepad, pp_weights, pp_stride);
	}

	PyObject *p = PyObject_CallMethod(this->get_pyobj(), (char *)"process_chunk", (char *)"OOOO", 
					  np_intensity.ptr, np_weights.ptr, np_pp_intensity.ptr, np_pp_weights.ptr);
	
	object ret(p, false);

	array2d_from_python(np_intensity, nfreq, nt_chunk, intensity, stride);
	array2d_from_python(np_weights, nfreq, nt_chunk, weights, stride);
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

    // Note: no tp_alloc() function is necessary, since the default (PyType_GenericAlloc()) calls memset().    
    static void tp_dealloc(PyObject *self_)
    {
	wi_transform_object *self = (wi_transform_object *)self_;

	delete self->pshared;
	self->pshared = nullptr;
	Py_TYPE(self)->tp_free((PyObject*) self);
    }

    // We do need tp_new() in order to define subclasses from Python
    static PyObject *tp_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
    {
	return type->tp_alloc(type, 0);
    }

    // Helper for get_pshared(): get a pointer-to-shared-ptr from a (PyObject *) which is known to be a wi_transform object.
    // This is the where the upcalling transform gets created, if it doesn't already exist.
    static inline shared_ptr<rf_pipelines::wi_transform> get_pshared(PyObject *obj)
    {
	wi_transform_object *t = (wi_transform_object *) obj;

	if (t->pshared)
	    return *(t->pshared);

	shared_ptr<rf_pipelines::wi_transform> p = make_shared<upcalling_wi_transform> (obj);
	t->pshared = new shared_ptr<rf_pipelines::wi_transform> (p);
	return p;
    }

    // Get a bare pointer from a (PyObject *) which is known to be a wi_transform_object
    static inline rf_pipelines::wi_transform *get_pbare(PyObject *obj)
    {
	shared_ptr<rf_pipelines::wi_transform> p = get_pshared(obj);
	
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
	get_pbare(self)->nfreq = int_from_python(value);
	return 0;
    }

    static PyObject *nt_chunk_getter(PyObject *self, void *closure)
    {
	return Py_BuildValue("i", get_pbare(self)->nt_chunk);
    }

    static int nt_chunk_setter(PyObject *self, PyObject *value, void *closure)
    {
	get_pbare(self)->nt_chunk = int_from_python(value);
	return 0;
    }

    static PyObject *nt_prepad_getter(PyObject *self, void *closure)
    {
	return Py_BuildValue("i", get_pbare(self)->nt_prepad);
    }

    static int nt_prepad_setter(PyObject *self, PyObject *value, void *closure)
    {
	get_pbare(self)->nt_prepad = int_from_python(value);
	return 0;
    }

    static PyObject *nt_postpad_getter(PyObject *self, void *closure)
    {
	return Py_BuildValue("i", get_pbare(self)->nt_postpad);
    }

    static int nt_postpad_setter(PyObject *self, PyObject *value, void *closure)
    {
	get_pbare(self)->nt_postpad = int_from_python(value);
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
    "MY AWESOME TRANSFORM",           /* tp_doc */
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


// make new python object (for factory functions)
static PyObject *make_wi_transform(const shared_ptr<rf_pipelines::wi_transform> &ptr)
{
    if (!ptr)
	throw runtime_error("rf_pipelines: internal error: empty pointer passed to make_wi_transform()");

    wi_transform_object *ret = PyObject_New(wi_transform_object, &wi_transform_type);
    if (!ret)
	return NULL;

    ret->pshared = new shared_ptr<rf_pipelines::wi_transform> (ptr);
    return (PyObject *) ret;
}


// -------------------------------------------------------------------------------------------------
//
// wi_stream wrapper class


struct wi_stream_object {
    PyObject_HEAD

    // "Borrowed" reference (delete will not be called in tp_dealloc())
    rf_pipelines::wi_stream *pbare;

    // Reference held through shared_ptr.  Using a pointer-to-a-shared pointer was least awkward here.
    shared_ptr<rf_pipelines::wi_stream> *pshared;

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
    static inline rf_pipelines::wi_stream *get_pbare(PyObject *obj)
    {
	wi_stream_object *s = (wi_stream_object *) obj;

	if (!s->pbare)
	    throw runtime_error("rf_pipelines: internal error: unexpected NULL pointer in wi_stream [should never happen]");

	return s->pbare;
    }
    
    static PyObject *run(PyObject *self, PyObject *arg)
    {
	rf_pipelines::wi_stream *stream = get_pbare(self);

	PyObject *iter = PyObject_GetIter(arg);
	if (!iter) {
	    PyErr_SetString(PyExc_RuntimeError, "rf_pipelines: expected argument to wi_stream.run() to be a list/iterator of wi_transform objects");	    
	    return NULL;
	}

	object iter_reference(iter, false);
	vector<object> item_references;
	vector<shared_ptr<rf_pipelines::wi_transform> > transform_list;

	for (;;) {
	    PyObject *item_ptr = PyIter_Next(iter);
	    if (!item_ptr)
		break;

	    item_references.push_back(object(item_ptr,false));

	    if (!PyObject_IsInstance(item_ptr, (PyObject *) &wi_transform_type))
		throw runtime_error("rf_pipelines: expected argument to wi_stream.run() to be a list/iterator of wi_transform objects");

	    transform_list.push_back(wi_transform_object::get_pshared(item_ptr));
	}	

	if (PyErr_Occurred())
	    throw python_exception();

	stream->run(transform_list);

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
    "MY AWESOME STREAM",           /* tp_doc */
    0,                         /* tp_traverse */
    0,                         /* tp_clear */
    0,                         /* tp_richcompare */
    0,                         /* tp_weaklistoffset */
    0,                         /* tp_iter */
    0,                         /* tp_iternext */
    wi_stream_methods,         /* tp_methods */
    0,                         /* tp_members */
    wi_stream_getseters,       /* tp_getset */
};


// make new python object (for factory functions)
static PyObject *make_wi_stream(const shared_ptr<rf_pipelines::wi_stream> &ptr)
{
    if (!ptr)
	throw runtime_error("rf_pipelines: internal error: empty pointer passed to make_wi_stream()");

    wi_stream_object *ret = PyObject_New(wi_stream_object, &wi_stream_type);
    if (!ret)
	return NULL;

    ret->pbare = ptr.get();
    ret->pshared = new shared_ptr<rf_pipelines::wi_stream> (ptr);
    return (PyObject *) ret;
}


struct dummy_stream : public rf_pipelines::wi_stream
{
    dummy_stream(const rf_pipelines::wi_stream &s)
    {
	this->nfreq = s.nfreq;
	this->nt_maxwrite = s.nt_maxwrite;
	this->freq_lo_MHz = s.freq_lo_MHz;
	this->freq_hi_MHz = s.freq_hi_MHz;
	this->dt_sample = s.dt_sample;
    }
    
    virtual void stream_body(rf_pipelines::wi_run_state &run_state)
    {
	throw runtime_error("rf_pipelines: attempt to call run() on 'dummy' stream used as argument to stream_start()");
    }
};

static PyObject *make_dummy_stream(const rf_pipelines::wi_stream &s)
{
    return make_wi_stream(make_shared<dummy_stream> (s));
}


// -------------------------------------------------------------------------------------------------
//
// Library


static PyObject *make_psrfits_stream(PyObject *self, PyObject *args)
{
    const char *filename = nullptr;
    if (!PyArg_ParseTuple(args, "s", &filename))
	return NULL;

    return make_wi_stream(rf_pipelines::make_psrfits_stream(filename));
}


static PyObject *make_simple_detrender(PyObject *self, PyObject *args)
{
    int nt_chunk = 0;
    if (!PyArg_ParseTuple(args, "i", &nt_chunk))
	return NULL;

    return make_wi_transform(rf_pipelines::make_simple_detrender(nt_chunk));
}


// -------------------------------------------------------------------------------------------------


static PyMethodDef module_methods[] = {
    { "make_psrfits_stream", tc_wrap2<make_psrfits_stream>, METH_VARARGS, "XXX" },
    { "make_simple_detrender", tc_wrap2<make_simple_detrender>, METH_VARARGS, "XXX" },
    { NULL, NULL, 0, NULL }
};


PyMODINIT_FUNC initrf_pipelines_c(void)
{
    import_array();

    if (PyType_Ready(&wi_stream_type) < 0)
        return;
    if (PyType_Ready(&wi_transform_type) < 0)
        return;

    PyObject *m = Py_InitModule3("rf_pipelines_c", module_methods, "XXX");
    if (!m)
	return;

    Py_INCREF(&wi_stream_type);
    PyModule_AddObject(m, "wi_stream", (PyObject *)&wi_stream_type);

    Py_INCREF(&wi_transform_type);
    PyModule_AddObject(m, "wi_transform", (PyObject *)&wi_transform_type);
}
