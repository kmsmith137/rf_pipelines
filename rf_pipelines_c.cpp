#include <Python.h>
#include "rf_pipelines.hpp"

using namespace std;


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

	if (self->pshared) {
	    delete self->pshared;
	    self->pshared = NULL;
	}

	Py_TYPE(self)->tp_free((PyObject*) self);
    }

    // Get a bare pointer from a (PyObject *) which is known to be a wi_stream_object
    static inline rf_pipelines::wi_stream *get_pbare(PyObject *obj)
    {
	wi_stream_object *s = (wi_stream_object *) obj;

	if (s->pbare)
	    return s->pbare;

	PyErr_SetString(PyExc_RuntimeError, "rf_pipelines: internal error: unexpected NULL pointer");
	return NULL;	
    }

    static PyObject *nfreq_getter(PyObject *self, void *closure)
    {
	rf_pipelines::wi_stream *p = wi_stream_object::get_pbare(self);
	return p ? Py_BuildValue("i", p->nfreq) : NULL;
    }

    static PyObject *nt_maxwrite_getter(PyObject *self, void *closure)
    {
	rf_pipelines::wi_stream *p = wi_stream_object::get_pbare(self);
	return p ? Py_BuildValue("i", p->nt_maxwrite) : NULL;
    }

    static PyObject *freq_lo_MHz_getter(PyObject *self, void *closure)
    {
	rf_pipelines::wi_stream *p = wi_stream_object::get_pbare(self);
	return p ? Py_BuildValue("d", p->freq_lo_MHz) : NULL;
    }

    static PyObject *freq_hi_MHz_getter(PyObject *self, void *closure)
    {
	rf_pipelines::wi_stream *p = wi_stream_object::get_pbare(self);
	return p ? Py_BuildValue("d", p->freq_hi_MHz) : NULL;
    }

    static PyObject *dt_sample_getter(PyObject *self, void *closure)
    {
	rf_pipelines::wi_stream *p = wi_stream_object::get_pbare(self);
	return p ? Py_BuildValue("d", p->dt_sample) : NULL;
    }
};


static PyGetSetDef wi_stream_getseters[] = {
    { (char *)"nfreq", wi_stream_object::nfreq_getter, NULL, NULL, NULL },
    { (char *)"nt_maxwrite", wi_stream_object::nt_maxwrite_getter, NULL, NULL, NULL },
    { (char *)"freq_lo_MHz", wi_stream_object::freq_lo_MHz_getter, NULL, NULL, NULL },
    { (char *)"freq_hi_MHz", wi_stream_object::freq_hi_MHz_getter, NULL, NULL, NULL },
    { (char *)"dt_sample", wi_stream_object::dt_sample_getter, NULL, (char *)"sample length in seconds", NULL },
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
    0,                         /* tp_methods */
    0,                         /* tp_members */
    wi_stream_getseters,       /* tp_getset */
};


// make new python object (for factory functions)
// can return NULL if something goes wrong
static PyObject *make_wi_stream(const shared_ptr<rf_pipelines::wi_stream> &ptr)
{
    if (!ptr) {
	PyErr_SetString(PyExc_RuntimeError, "rf_pipelines: internal error: empty pointer passed to make_wi_stream()");
	return NULL;
    }

    wi_stream_object *ret = PyObject_New(wi_stream_object, &wi_stream_type);
    if (!ret)
	return NULL;

    try {
	ret->pbare = ptr.get();
	ret->pshared = new shared_ptr<rf_pipelines::wi_stream> (ptr);
    }
    catch (...) {
	Py_XDECREF(ret);
	return PyErr_NoMemory();
    }

    return (PyObject *) ret;
}


// -------------------------------------------------------------------------------------------------
//
// Wrap wi_transform


// "Upcalling" transform whose virtual functions are implemented by python upcalls.
struct upcalling_wi_transform : public rf_pipelines::wi_transform
{
    virtual void set_stream(const rf_pipelines::wi_stream &stream)
    {
	throw runtime_error("oops no set_stream");
    }

    virtual void start_substream(double t0)
    {
	throw runtime_error("oops no start_substream");
    }

    virtual void process_chunk(float *intensity, float *weight, int stride, float *pp_intensity, float *pp_weight, int pp_stride)
    {
	throw runtime_error("oops no process_chunk");
    }

    virtual void end_substream()
    {
	throw runtime_error("oops no end_substream");
    }
};


struct wi_transform_object {
    PyObject_HEAD

    shared_ptr<rf_pipelines::wi_transform> *pshared;

    // note: no tp_alloc() function is necessary, since the default (PyType_GenericAlloc()) calls memset().    
    static void tp_dealloc(PyObject *self_)
    {
	wi_transform_object *self = (wi_transform_object *)self_;

	if (self->pshared) {
	    delete self->pshared;
	    self->pshared = NULL;
	}

	Py_TYPE(self)->tp_free((PyObject*) self);
    }

    // we do need tp_new() in order to define subclasses from Python
    static PyObject *tp_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
    {
	return type->tp_alloc(type, 0);
    }

    // Get a shared_ptr from a (PyObject *) which is known to be a wi_transform_object
    static inline shared_ptr<rf_pipelines::wi_transform> get_pshared(PyObject *obj)
    {
	wi_transform_object *t = (wi_transform_object *) obj;

	if (!t->pshared) {
	    try {
		shared_ptr<rf_pipelines::wi_transform> p = make_shared<upcalling_wi_transform> ();
		t->pshared = new shared_ptr<rf_pipelines::wi_transform> (p);
	    } 
	    catch (...) {
		PyErr_NoMemory();
		return shared_ptr<rf_pipelines::wi_transform> ();
	    }
	}

	shared_ptr<rf_pipelines::wi_transform> ret = *(t->pshared);
	if (!ret)
	    PyErr_SetString(PyExc_RuntimeError, "rf_pipelines: internal error: unexpected NULL pointer");

	return ret;
    }

    // Get a bare pointer from a (PyObject *) which is known to be a wi_transform_object
    static inline rf_pipelines::wi_transform *get_pbare(PyObject *obj)
    {
	shared_ptr<rf_pipelines::wi_transform> p = get_pshared(obj);
	return p ? p.get() : NULL;
    }

    static PyObject *nfreq_getter(PyObject *self, void *closure)
    {
	rf_pipelines::wi_transform *p = get_pbare(self);
	return p ? Py_BuildValue("i", p->nfreq) : NULL;
    }

    static int nfreq_setter(PyObject *self, PyObject *value, void *closure)
    {
	long n = PyInt_AsLong(value);
	if ((n == -1) && PyErr_Occurred())
	    return -1;

	rf_pipelines::wi_transform *p = get_pbare(self);
	if (!p)
	    return -1;

	p->nfreq = n;
	return 0;
    }


    static PyObject *nt_chunk_getter(PyObject *self, void *closure)
    {
	rf_pipelines::wi_transform *p = get_pbare(self);
	return p ? Py_BuildValue("i", p->nt_chunk) : NULL;
    }

    static int nt_chunk_setter(PyObject *self, PyObject *value, void *closure)
    {
	long n = PyInt_AsLong(value);
	if ((n == -1) && PyErr_Occurred())
	    return -1;

	rf_pipelines::wi_transform *p = get_pbare(self);
	if (!p)
	    return -1;

	p->nt_chunk = n;
	return 0;
    }


    static PyObject *nt_prepad_getter(PyObject *self, void *closure)
    {
	rf_pipelines::wi_transform *p = get_pbare(self);
	return p ? Py_BuildValue("i", p->nt_prepad) : NULL;
    }

    static int nt_prepad_setter(PyObject *self, PyObject *value, void *closure)
    {
	long n = PyInt_AsLong(value);
	if ((n == -1) && PyErr_Occurred())
	    return -1;

	rf_pipelines::wi_transform *p = get_pbare(self);
	if (!p)
	    return -1;

	p->nt_prepad = n;
	return 0;
    }


    static PyObject *nt_postpad_getter(PyObject *self, void *closure)
    {
	rf_pipelines::wi_transform *p = get_pbare(self);
	return p ? Py_BuildValue("i", p->nt_postpad) : NULL;
    }

    static int nt_postpad_setter(PyObject *self, PyObject *value, void *closure)
    {
	long n = PyInt_AsLong(value);
	if ((n == -1) && PyErr_Occurred())
	    return -1;

	rf_pipelines::wi_transform *p = get_pbare(self);
	if (!p)
	    return -1;

	p->nt_postpad = n;
	return 0;
    }
};


static PyGetSetDef wi_transform_getseters[] = {
    { (char *)"nfreq", wi_transform_object::nfreq_getter, wi_transform_object::nfreq_setter, NULL, NULL },
    { (char *)"nt_chunk", wi_transform_object::nt_chunk_getter, wi_transform_object::nt_chunk_setter, NULL, NULL },
    { NULL, NULL, NULL, NULL, NULL }
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
    0,                         /* tp_methods */
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
// can return NULL if something goes wrong
static PyObject *make_wi_transform(const shared_ptr<rf_pipelines::wi_transform> &ptr)
{
    if (!ptr) {
	PyErr_SetString(PyExc_RuntimeError, "rf_pipelines: internal error: empty pointer passed to make_wi_transform()");
	return NULL;
    }

    wi_transform_object *ret = PyObject_New(wi_transform_object, &wi_transform_type);
    if (!ret)
	return NULL;

    try {
	ret->pshared = new shared_ptr<rf_pipelines::wi_transform> (ptr);
    }
    catch (...) {
	Py_XDECREF(ret);
	return PyErr_NoMemory();
    }

    return (PyObject *) ret;
}


// -------------------------------------------------------------------------------------------------
//
// Library


static PyObject *make_psrfits_stream(PyObject *self, PyObject *args)
{
    const char *filename = nullptr;
    if (!PyArg_ParseTuple(args, "s", &filename))
	return NULL;

    try {
	shared_ptr<rf_pipelines::wi_stream> ret = rf_pipelines::make_psrfits_stream(filename);
	return make_wi_stream(ret);
    }
    catch (std::exception &e) {
	PyErr_SetString(PyExc_RuntimeError, e.what());
	return NULL;
    }
    catch (...) {
	return NULL;
    }
}


static PyObject *make_simple_detrender(PyObject *self, PyObject *args)
{
    int nt_chunk = 0;
    if (!PyArg_ParseTuple(args, "i", &nt_chunk))
	return NULL;

    try {
	shared_ptr<rf_pipelines::wi_transform> ret = rf_pipelines::make_simple_detrender(nt_chunk);
	return make_wi_transform(ret);
    }
    catch (std::exception &e) {
	PyErr_SetString(PyExc_RuntimeError, e.what());
	return NULL;
    }
    catch (...) {
	return NULL;
    }
}


// -------------------------------------------------------------------------------------------------


static PyMethodDef module_methods[] = {
    { "make_psrfits_stream", make_psrfits_stream, METH_VARARGS, "XXX" },
    { "make_simple_detrender", make_simple_detrender, METH_VARARGS, "XXX" },
    { NULL, NULL, 0, NULL }
};


PyMODINIT_FUNC initrf_pipelines_c(void)
{
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
