#ifndef _RF_PIPELINES_CYTHON_HPP
#define _RF_PIPELINES_CYTHON_HPP

#include <Python.h>
#include "rf_pipelines_internals.hpp"


using namespace std;


// -------------------------------------------------------------------------------------------------
//
// Transform and stream objects


struct _wi_transform {
    shared_ptr<rf_pipelines::wi_transform> p;

    _wi_transform(const shared_ptr<rf_pipelines::wi_transform> &p_) : p(p_)
    {
	rf_assert(p_);
    }
};


struct _wi_stream {
    shared_ptr<rf_pipelines::wi_stream> p;
    vector<shared_ptr<rf_pipelines::wi_transform> > tbuf;

    _wi_stream(const shared_ptr<rf_pipelines::wi_stream> &p_) : p(p_)
    {
	rf_assert(p_);
    }

    int get_nfreq() { return p->nfreq; }
    int get_nt_maxwrite() { return p->nt_maxwrite; }
    double get_freq_lo_MHz() { return p->freq_lo_MHz; }
    double get_freq_hi_MHz() { return p->freq_hi_MHz; }
    double get_dt_sample() { return p->dt_sample; }    

    void clear_transforms() { tbuf.clear(); }

    void add_transform(_wi_transform *t)
    {
	rf_assert(t != nullptr);
	tbuf.push_back(t->p);
    }
    
    void run() { p->run(tbuf); }
};


// ------------------------------------------------------------------------------------------------
//
// Library


inline _wi_stream *_make_psrfits_stream(string filename)
{
    return new _wi_stream(rf_pipelines::make_psrfits_stream(filename));
}


inline _wi_transform *_make_simple_detrender(int nt_chunk)
{
    return new _wi_transform(rf_pipelines::make_simple_detrender(nt_chunk));
}


// -------------------------------------------------------------------------------------------------
//
// upcalling_transform: a wi_transform whose virtual functions call "up" from C++ to Python


struct upcalling_transform : public rf_pipelines::wi_transform {
    PyObject *weakref;

    upcalling_transform(PyObject *pyobj)
    {
	if (!pyobj)
	    throw runtime_error("rf_pipelines: NULL pointer in upcalling_transform constructor?!");
	
	this->weakref = PyWeakref_NewRef(pyobj, NULL);
	if (!weakref)
	    throw runtime_error("rf_pipelines: PyWeakref_NewRef() failed?!");
    }

    virtual ~upcalling_transform()
    { 
	if (weakref) {
	    Py_XDECREF(weakref);
	    weakref = NULL;
	}
    }

    // Note: returns a borrowed (not new) reference
    PyObject *get_pyobj()
    {
	PyObject *ret = PyWeakref_GetObject(this->weakref);
	if (!ret)
	    throw runtime_error("rf_pipelines: PyWeakref_GetObject() failed?!");    
	if (ret == Py_None)
	    throw runtime_error("rf_pipelines: weak reference expired?! [should never happen]");
	return ret;
    }

    // Returns a new reference
    PyObject *getattr(const char *attr_name)
    {
	PyObject *ret = PyObject_GetAttrString(this->get_pyobj(), attr_name);
	if (!ret)
	    throw runtime_error("rf_pipelines: wi_transform object has no attribute named '" + string(attr_name) + "'");
	return ret;
    }

    virtual void set_stream(const rf_pipelines::wi_stream &stream)
    {
	PyObject *f = this->getattr("set_stream");
	cerr << "XXX yay!\n";	
	Py_XDECREF(f);
	throw runtime_error("XXX set_stream() blowing up now\n");
    }

    virtual void start_substream(double t0)
    {
	throw runtime_error("XXX start_substream() blowing up now\n");
    }

    virtual void process_chunk(float *intensity, float *weight, int stride, float *pp_intensity, float *pp_weight, int pp_stride)
    {
	throw runtime_error("XXX process_chunk() blowing up now\n");
    }

    virtual void end_substream()
    {	
	throw runtime_error("XXX end_substream() blowing up now\n");
    }
};

inline _wi_transform *_make_upcalling_transform(PyObject *pyobj)
{
    return new _wi_transform(make_shared<upcalling_transform> (pyobj));
}


#endif  // _RF_PIPELINES_CYTHON
