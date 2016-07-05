// Note: I haven't systematically documented the C++ interface to rf_pipelines,
// so the level of documentation will be hit-or-miss!

#ifndef _PYTHON_EXTENSION_HELPERS_HPP
#define _PYTHON_EXTENSION_HELPERS_HPP

#include <Python.h>
#include <frameobject.h>
#include <numpy/arrayobject.h>

#include <cstdlib>
#include <iostream>
#include <stdexcept>


// -------------------------------------------------------------------------------------------------
//
// Here is something that seems like it should be part of the Python C-API but isn't!  
// A standalone function
//
//    const char *get_exception_text(int &free_flag)
//
// which returns a C-style string containing a textual description of the current exception.  If the 
// returned string was malloc()-ed, then 'free_flag' will be set, to tell the caller that free() is needed.
//
// We use C strings instead of std::string, so that we are guaranteed never to throw a C++ exception
// if the process is out of memory.  (In this case we would return "Out of Memory" for the exception
// text.)  
//
//   PyErr_Display()           [ pythonrun.c ]
//      PyTraceBack_Print()    [ traceback.c ]
//      parse_syntax_error()   [ pythonrun.c ]


//
// get_traceback_text(): helper for get_exception_text()
//
// Reference: PyTraceBack_Print() in Python/traceback.c
//
// FIXME room for improvement here
//   - It gets confused by callbacks, but works well enough to give the filename and line number
//   - There is some signal handling stuff in PyTraceBack_Print() which I didn't try to emulate but probably should
//   - I didn't dig up the actual line of text in the source file (by emalating _Py_DisplaySourceLine() in traceback.c), 
//     I just took whatever information was available in the traceback object.
//
static const char *get_traceback_text(PyObject *ptraceback, int &free_flag)
{
    free_flag = 0;

    if (!ptraceback)
	return "No traceback available [traceback object was NULL, this is probably some sort of internal error]\n";
    if (!PyTraceBack_Check(ptraceback))
	return "No traceback available [PyTraceBack_Check() returned false, this is probably some sort of internal error]\n";

    static const int tb_maxcount = 10;
    static const int tb_maxline = 1000;
    static const int tb_nheader = 100;
    static const int nalloc = tb_nheader + (tb_maxcount * tb_maxline) + 1;
    
    char *ret = (char *)malloc(nalloc);
    if (!ret)
	return "No traceback available [out of memory]\n";
	
    free_flag = 1;
    memset(ret, 0, nalloc);
    int istr = snprintf(ret, tb_nheader, "Traceback (most recent call last):\n");
    int itb = 0;
    
    for (PyTracebackObject *tb = (PyTracebackObject *)ptraceback; tb; tb = tb->tb_next) {
	if (itb >= tb_maxcount)
	    break;

	const char *filename = PyString_AsString(tb->tb_frame->f_code->co_filename);
	const char *code = PyString_AsString(tb->tb_frame->f_code->co_name);
	int lineno = tb->tb_lineno;

	if (!filename)
	    filename = "[filename unavailable]";
	if (!code)
	    filename = "[code unavailable]";

	istr += snprintf(ret+istr, tb_maxline, "  File %s, line %d\n    %s\n", filename, lineno, code);
	itb++;
    }

    return ret;
}


//
// get_exception_type_text(): helper for get_exception_text()
//
// Reference: PyErr_Display() in Python/pythonrun.c
//
// FIXME room for improvement here
//   - In PyErr_Display(), there is a magic print_file_and_line() method which I didn't look into
//   - In PyErr_Display(), there is module lookup logic which I didn't look into
//
static const char *get_exception_type_text(PyObject *ptype, int &free_flag)
{
    free_flag = 0;

    if (!ptype)
       return "[unknown exception type]";
    if (!PyExceptionClass_Check(ptype))
       return "[unknown exception type]";	

    char *s = PyExceptionClass_Name(ptype);
    if (!s)
       return "[unknown exception type]";

    char *dot = strrchr(s, '.');
    if (dot != NULL)
	s = dot+1;
    
    const char *ret = strdup(s);
    if (!ret)
	return "MemoryError";

    free_flag = 1;
    return ret;
}


//
// get_exception_type_text(): helper for get_exception_text()
//
// Reference: PyErr_Display() in Python/pythonrun.c
//
// FIXME room for improvement here
//   - In PyErr_Display(), there is a PyFile_WriteObject(..., Py_PRINT_RAW) path
//
static const char *get_exception_value_text(PyObject *pvalue, int &free_flag)
{
    free_flag = 0;

    if (!pvalue)
	return "[exception value unavailable]";

    PyObject *s = PyObject_Str(pvalue);
    if (!s)
	return "[exception value unavailable]";

    const char *s2 = PyString_AsString(s);
    if (!s2) {
	Py_XDECREF(s);
	return "[unavailable exception value]";
    }
	
    char *ret = strdup(s2);
    if (!ret) {
	Py_XDECREF(s);
	return "[out of memory]";
    }

    free_flag = 1;
    return ret;
}


// Using the three helper functions above, here is get_exception_text()
static const char *get_exception_text(int &free_flag)
{
    free_flag = 0;

    PyObject *ptype = PyErr_Occurred();  // returns borrowed reference
    
    if (!ptype)
	return "Unknown exception occurred [maybe internal error? get_exception_text() was called, but PyErr_Occurred() returned NULL]";
    if (PyErr_GivenExceptionMatches(ptype, PyExc_MemoryError))
	return "Out of memory";  // special handling of out-of-memory case, to avoid any malloc()-ing

    ptype = NULL;
    PyObject *pvalue = NULL;
    PyObject *ptraceback = NULL;

    // Caller of PyErr_Fetch() owns all 3 references.  Calling PyErr_Fetch() clears the interpreter error state!
    PyErr_Fetch(&ptype, &pvalue, &ptraceback);

    int ptype_free_flag = 0;
    int pvalue_free_flag = 0;
    int ptraceback_free_flag = 0;

    const char *ptype_text = get_exception_type_text(ptype, ptype_free_flag);
    const char *pvalue_text = get_exception_value_text(pvalue, pvalue_free_flag);
    const char *ptraceback_text = get_traceback_text(ptraceback, ptraceback_free_flag);

    // Undoes PyErr_Fetch(), by setting the interpreter error state and returning the 3 references.
    PyErr_Restore(ptype, pvalue, ptraceback);

    int nalloc = strlen(ptype_text) + strlen(pvalue_text) + strlen(ptraceback_text) + 10;
    char *ret = (char *) malloc(nalloc);

    if (ret) {
	snprintf(ret, nalloc, "%s%s: %s", ptraceback_text, ptype_text, pvalue_text);
	free_flag = 1;
    }
    else
	ret = (char *) "Out of memory";

    if (ptype_free_flag)
	free((void *) ptype_text);
    if (pvalue_free_flag)
	free((void *) pvalue_text);
    if (ptraceback_free_flag)
	free((void *) ptraceback_text);

    return ret;
}


// -------------------------------------------------------------------------------------------------
//
// Exception handling



// The following exception is thrown in C++ code if the python interpreter is found to be in its error state.

struct python_exception : std::exception
{
    const char *msg;   // NULL means "out of memory"
    int free_flag;
    
    // a python_exception 
    python_exception() : msg(NULL), free_flag(false)
    {
	this->msg = get_exception_text(this->free_flag);
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

	if (msg)  // FIXME is there a better way to handle this case?
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
	if (increment_refcount)
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

    ssize_t get_refcount() const { return Py_REFCNT(ptr); }

    // FIXME not currently using this, since there are corner cases I don't understand (e.g. control-C can cause callback to 
    // keep the reference, maybe it's in the traceback stack?)
    void die_unless_refcount1(const char *msg)
    {
	if (this->get_refcount() > 1) {
	    std::cerr << msg << "\nThis is such a grievous sin that this process will be terminated now, without even throwing an exception!\n";
	    exit(1);
	}
    }
};


// -------------------------------------------------------------------------------------------------
//
// Misc


inline ssize_t ssize_t_from_python(PyObject *obj)
{
    ssize_t n = PyInt_AsSsize_t(obj);
    if ((n == -1) && PyErr_Occurred())
	throw python_exception();
    return n;
}


#endif  // _PYTHON_EXTENSION_HELPERS_HPP
