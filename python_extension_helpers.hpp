#ifndef _PYTHON_EXTENSION_HELPERS_HPP
#define _PYTHON_EXTENSION_HELPERS_HPP

#include <Python.h>
#include <numpy/arrayobject.h>

#include <cstdlib>
#include <stdexcept>


// -------------------------------------------------------------------------------------------------
//
// Exception handling
//


// The following exception is thrown in C++ code if the python interpreter is found to be in its error state.


struct python_exception : std::exception
{
    const char *msg;   // NULL means "out of memory"
    bool free_flag;
    
    // a python_exception should be thrown with no 
    python_exception() : msg(NULL), free_flag(false)
    {
	PyObject *ptype = PyErr_Occurred();

	if (!ptype) {
	    this->msg = "python_exception thrown but PyErr_Occurred() returned false, not sure what's going on here";
	    return;
	}

	// special handling of out-of-memory case, to avoid any new memory-allocating operations
	if (PyErr_GivenExceptionMatches(ptype, PyExc_MemoryError))
	    return;

	ptype = NULL;
	PyObject *pvalue = NULL;
	PyObject *ptraceback = NULL;

	PyErr_Fetch(&ptype, &pvalue, &ptraceback);
	Py_XDECREF(ptype);
	Py_XDECREF(ptraceback);

	if (!pvalue) {
	    this->msg = "python_exception thrown but PyErr_Fetch() failed, not sure what's going on here";
	    return;
	}

	char *s = PyString_AsString(pvalue);
	Py_XDECREF(pvalue);

	if (!s) {
	    this->msg = "python_exception thrown but PyString_AsString() failed, not sure what's going on here";
	    return;
	}
	
	s = strdup(s);
	if (!s)
	    return;

	this->msg = s;
	this->free_flag = true;
    }

    ~python_exception()
    {
	if (msg && free_flag)
	    free((void *)msg);
	msg = nullptr;
	free_flag = false;
    }
    
    // make noncopyable but implement move constructor
    python_exception(const python_exception &) = delete;
    python_exception& operator=(const python_exception &) = delete;

    python_exception(python_exception &&x) noexcept
    {
	this->msg = x.msg;
	this->free_flag = x.free_flag;
	
	// If this message ever appears, then I don't understand C++11 as well as I thought..
	x.msg = "python exception thrown, but text was destroyed by C++11 move constructor [should never happen]";
	x.free_flag = false;
    }

    virtual char const *what() const noexcept
    {
	return msg;
    }

    void ensure_python_error_set()
    {
	if (PyErr_Occurred())
	    return;   // python interpreter is already in its error state, no need to do anything
	
	if (msg)
	    PyErr_SetString(PyExc_RuntimeError, this->msg);   // Note: PyErr_SetString() will copy the string
	else
	    PyErr_NoMemory();
    }
};


//
// Sometimes we need to map a C++ exception to a Python exception.  We do this by wrapping the
// Python C-API function with one of the following templates ("tc_wrap" = try-catch wrap)
//
//   tc_wrapN<f>: wraps a function of the form (for N=1,2,3)
//     PyObject *f(PyObject *arg1, PyObject *arg2, ..., PyObject *argN)
//
//   tc_wrap_getter<f>: wraps a "getter" function of the form 
//     PyObject *f(PyObject *self, void *closure)
//
//   tc_wrap_setter<f>: wraps a "setter" function of the form
//     int f(PyObject *self, PyObject *value, void *closure)
//
// FIXME: hmm there must be some C++ magic which avoids cut-and-paste here.  I didn't spend much time 
// trying to figure it out, decided to return to it in a future cleanup iteration.
//


template<PyObject* (*f)(PyObject *)>
inline PyObject *tc_wrap1(PyObject *arg)
{
    try {
	return f(arg);
    } catch (python_exception &e) {
	e.ensure_python_error_set();
	return NULL;
    } catch (std::exception &e) {
	PyErr_SetString(PyExc_RuntimeError, e.what());
	return NULL;
    } catch (...) {
	PyErr_SetString(PyExc_RuntimeError, "unknown C++ exception was thrown");
	return NULL;
    }    
}

template<PyObject* (*f)(PyObject *, PyObject *)>
inline PyObject *tc_wrap2(PyObject *arg1, PyObject *arg2)
{
    try {
	return f(arg1, arg2);
    } catch (python_exception &e) {
	e.ensure_python_error_set();
	return NULL;
    } catch (std::exception &e) {
	PyErr_SetString(PyExc_RuntimeError, e.what());
	return NULL;
    } catch (...) {
	PyErr_SetString(PyExc_RuntimeError, "unknown C++ exception was thrown");
	return NULL;
    }    
}

template<PyObject* (*f)(PyObject *, PyObject *, PyObject *)>
inline PyObject *tc_wrap3(PyObject *arg1, PyObject *arg2, PyObject *arg3)
{
    try {
	return f(arg1, arg2, arg3);
    } catch (python_exception &e) {
	e.ensure_python_error_set();
	return NULL;
    } catch (std::exception &e) {
	PyErr_SetString(PyExc_RuntimeError, e.what());
	return NULL;
    } catch (...) {
	PyErr_SetString(PyExc_RuntimeError, "unknown C++ exception was thrown");
	return NULL;
    }    
}

template<PyObject* (*f)(PyObject *, void *)>
inline PyObject *tc_wrap_getter(PyObject *self, void *closure)
{
    try {
	return f(self, closure);
    } catch (python_exception &e) {
	e.ensure_python_error_set();
	return NULL;
    } catch (std::exception &e) {
	PyErr_SetString(PyExc_RuntimeError, e.what());
	return NULL;
    } catch (...) {
	PyErr_SetString(PyExc_RuntimeError, "unknown C++ exception was thrown");
	return NULL;
    }
}

template<int (*f)(PyObject *, PyObject *, void *)>
inline int tc_wrap_setter(PyObject *self, PyObject *value, void *closure)
{
    try {
	return f(self, value, closure);
    } catch (python_exception &e) {
	e.ensure_python_error_set();
	return -1;
    } catch (std::exception &e) {
	PyErr_SetString(PyExc_RuntimeError, e.what());
	return -1;
    } catch (...) {
	PyErr_SetString(PyExc_RuntimeError, "unknown C++ exception was thrown");
	return -1;
    }
}


// -------------------------------------------------------------------------------------------------
//
// A simple 'object' class to help get reference counting correct (this is especially tricky in
// paths which can throw exceptions).  Another use for the object class is to streamline exception
// handling, since the object class automatically throws an exception if a NULL pointer arises.
//
// FIXME does it make sense to implement a C++11 move constructor here?


struct object {
    // can never be NULL
    PyObject *ptr;

    object(PyObject *p, bool increment_refcount)
    {
	this->ptr = p;

	if (!p)
	    throw python_exception();
	if (!increment_refcount)
	    Py_INCREF(p);
    }

    // if an 'object' is default-constructed, it points to Py_None
    object()
    {
	Py_INCREF(Py_None);
	this->ptr = Py_None;
    }

    ~object()
    {
	Py_XDECREF(ptr);
	ptr = NULL;
    }

    object(const object &x) : object(x.ptr,true) { }

    object &operator=(const object &x)
    {
	// this ordering of calls handles the self-assignment case correctly
	Py_XINCREF(x.ptr);
	Py_XDECREF(this->ptr);
	this->ptr = x.ptr;
	return *this;
    }
};


// -------------------------------------------------------------------------------------------------
//
// Misc


inline int int_from_python(PyObject *obj)
{
    int n = PyInt_AsLong(obj);
    if ((n == -1) && PyErr_Occurred())
	throw python_exception();
    return n;
}


#endif  // _PYTHON_EXTENSION_HELPERS_HPP
