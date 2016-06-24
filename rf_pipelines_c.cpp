// g++ -std=c++11 -fPIC -Wall -I/usr/include/python2.7 -shared -L. -L/home/kmsmith/lib -o rf_pipelines_c.so rf_pipelines_c.cpp -lrf_pipelines -lpsrfits_utils -lcfitsio

#include <Python.h>
#include "rf_pipelines.hpp"


using namespace std;


// -------------------------------------------------------------------------------------------------
//
// Wrap wi_stream


struct wi_stream_object {
    PyObject_HEAD

    // we use a pointer-to-a-shared-pointer, to avoid figuring out how to handle the shared_ptr in tp_alloc()/tp_dealloc().
    shared_ptr<rf_pipelines::wi_stream> *p;


    // note: no tp_alloc function is necessary, since the default (PyType_GenericAlloc()) calls memset().    
    static void tp_dealloc(wi_stream_object *self)
    {
	if (self->p) {
	    delete self->p;
	    self->p = NULL;
	}
	Py_TYPE(self)->tp_free((PyObject*) self);
    }


    // get bare pointer to underlying wi_stream
    // can return NULL if something goes wrong
    inline rf_pipelines::wi_stream *get()
    {
	if (!p) {
	    PyErr_SetString(PyExc_RuntimeError, "rf_pipelines: internal error: unexpected NULL pointer");
	    return NULL;
	}
	
	shared_ptr<rf_pipelines::wi_stream> &sp = *p;
	if (!sp) {
	    PyErr_SetString(PyExc_RuntimeError, "rf_pipelines: internal error: unexpected NULL pointer");
	    return NULL;
	}

	return sp.get();
    }


    static PyObject *nfreq_getter(wi_stream_object *self, void *closure)
    {
	rf_pipelines::wi_stream *s = self->get();
	return Py_BuildValue("i", s->nfreq);
    }

    static PyObject *nt_maxwrite_getter(wi_stream_object *self, void *closure)
    {
	rf_pipelines::wi_stream *s = self->get();
	return Py_BuildValue("i", s->nt_maxwrite);
    }

    static PyObject *freq_lo_MHz_getter(wi_stream_object *self, void *closure)
    {
	rf_pipelines::wi_stream *s = self->get();
	return Py_BuildValue("d", s->freq_lo_MHz);
    }

    static PyObject *freq_hi_MHz_getter(wi_stream_object *self, void *closure)
    {
	rf_pipelines::wi_stream *s = self->get();
	return Py_BuildValue("d", s->freq_hi_MHz);
    }

    static PyObject *dt_sample_getter(wi_stream_object *self, void *closure)
    {
	rf_pipelines::wi_stream *s = self->get();
	return Py_BuildValue("d", s->dt_sample);
    }
};


static PyGetSetDef wi_stream_getseters[] = {
    { (char *)"nfreq", (getter) &wi_stream_object::nfreq_getter, NULL, NULL, NULL },
    { (char *)"nt_maxwrite", (getter) &wi_stream_object::nt_maxwrite_getter, NULL, NULL, NULL },
    { (char *)"freq_lo_MHz", (getter) &wi_stream_object::freq_lo_MHz_getter, NULL, NULL, NULL },
    { (char *)"freq_hi_MHz", (getter) &wi_stream_object::freq_hi_MHz_getter, NULL, NULL, NULL },
    { (char *)"dt_sample", (getter) &wi_stream_object::dt_sample_getter, NULL, (char *)"sample length in seconds", NULL },
    { NULL, NULL, NULL, NULL, NULL }
};


static PyTypeObject wi_stream_type = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "rf_pipelines_c.wi_stream",  /* tp_name */
    sizeof(wi_stream_object),    /* tp_basicsize */
    0,                         /* tp_itemsize */
    (destructor) wi_stream_object::tp_dealloc,  /* tp_dealloc */
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
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    0,                         /* tp_init */
    0,                         /* tp_alloc */
    0,                         /* tp_new */
};


// make new python object (for factory functions)
// can return NULL if something goes wrong
static PyObject *make_wi_stream(const shared_ptr<rf_pipelines::wi_stream> &ptr)
{
    if (!ptr) {
	PyErr_SetString(PyExc_RuntimeError, "rf_pipelines: internal error: unexpected NULL pointer");
	return NULL;
    }

    wi_stream_object *ret = PyObject_New(wi_stream_object, &wi_stream_type);
    if (!ret)
	return NULL;

    try {
	ret->p = new shared_ptr<rf_pipelines::wi_stream> (ptr);
    }
    catch (...) {
	Py_XDECREF(ret);
	return PyErr_NoMemory();  // always returns NULL
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


// -------------------------------------------------------------------------------------------------


static PyMethodDef module_methods[] = {
    { "make_psrfits_stream", make_psrfits_stream, METH_VARARGS, "XXX" },
    { NULL, NULL, 0, NULL }
};


extern "C" PyMODINIT_FUNC initrf_pipelines_c(void)
{
    if (PyType_Ready(&wi_stream_type) < 0)
        return;

    PyObject *m = Py_InitModule3("rf_pipelines_c", module_methods, "XXX");

    if (!m)
	return;

    Py_INCREF(&wi_stream_type);
    PyModule_AddObject(m, "wi_stream", (PyObject *)&wi_stream_type);
}
