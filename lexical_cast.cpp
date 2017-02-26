// This file defines
//
//   template<typename T> T lexical_cast(const std::string &x);
//
// which converts a string to type T.  Currently the following types T are supported:
//
//   string  (trivial)
//   long
//   int
//   double
//   float
//   uint16_t
//
// Note that a more general implementation of lexical_cast is already defined in boost, and considered 
// for inclusion in C++ TR2, but it's not in C++11, so we're forced to define it here!

#include "rf_pipelines_internals.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif


template<> const char *typestr<string>()    { return "string"; }
template<> const char *typestr<long>()      { return "long"; }
template<> const char *typestr<int>()       { return "int"; }
template<> const char *typestr<double>()    { return "double"; }
template<> const char *typestr<float>()     { return "float"; }
template<> const char *typestr<uint16_t>()  { return "uint16_t"; }
template<> const char *typestr<bool>()      { return "bool"; }


// trivial case: convert string -> string
template<> bool lexical_cast(const string &x, string &ret) { ret = x; return true;}


inline bool is_all_spaces(const char *s)
{
    if (!s)
	throw runtime_error("fatal: NULL pointer passed to is_all_spaces()");

    for (;;) {
	if (!*s)
	    return true;
	if (!isspace(*s))
	    return false;
	s++;
    }
}


template<> bool lexical_cast(const string &x, long &ret)
{ 
    const char *ptr = x.c_str();
    char *endptr = NULL;

    ret = strtol(ptr, &endptr, 10);
    return (endptr != ptr) && (ret != LONG_MIN) && (ret != LONG_MAX) && is_all_spaces(endptr);
}


template<> bool lexical_cast(const string &x, int &ret)
{
    long retl;

    if (!lexical_cast(x, retl))
	return false;
    if ((sizeof(int) != sizeof(long)) && ((retl < INT_MIN) || (retl > INT_MAX)))
	return false;

    ret = retl;
    return true;
}


template<> bool lexical_cast(const string &x, uint16_t &ret)
{
    long retl;

    if (!lexical_cast(x, retl))
	return false;
    if ((retl < 0) || (retl > 65535))
	return false;

    ret = retl;
    return true;
}


template<> bool lexical_cast(const string &x, double &ret)
{ 
    const char *ptr = x.c_str();
    char *endptr = NULL;

    ret = strtod(ptr, &endptr);
    return (endptr != ptr) && (ret != -HUGE_VAL) && (ret != HUGE_VAL) && is_all_spaces(endptr);
}


template<> bool lexical_cast(const string &x, float &ret)
{ 
    const char *ptr = x.c_str();
    char *endptr = NULL;

    ret = strtof(ptr, &endptr);
    return (endptr != ptr) && (ret != -HUGE_VALF) && (ret != HUGE_VALF) && is_all_spaces(endptr);
}

template<> bool lexical_cast(const string &x, bool &ret)
{
    const char *ptr = x.c_str();

    if (!strcasecmp(ptr,"t") || !strcasecmp(ptr,"true"))
	ret = true;
    else if (!strcasecmp(ptr,"f") || !strcasecmp(ptr,"false"))
	ret = false;
    else
	return false;

    return true;
}


// -------------------------------------------------------------------------------------------------
//
// Unit test


template<typename T> 
static void check_convert(const string &x, T y)
{
    T ret;
    if (!lexical_cast(x, ret))
	throw runtime_error("test_lexical_cast(): didn't successfully convert");
    if (fabs(double(ret) - double(y)) > 1.0e-5)
	throw runtime_error("test_lexical_cast(): didn't correctly convert");
}


template<typename T> 
static void check_convert_fails(const string &x)
{
    T ret;
    if (lexical_cast(x, ret))
	throw runtime_error("test_lexical_cast(): conversion succeeded where it was expected to fail");
}


void test_lexical_cast()
{
    check_convert<int>("0", 0);
    check_convert<int>("-0", 0);
    check_convert<int>("12", 12);
    check_convert<int>("-123", -123);
    check_convert<int>(" \t 1234  \n\t", 1234);

    check_convert_fails<int>("");
    check_convert_fails<int>("  ");
    check_convert_fails<int>("oops");
    check_convert_fails<int>(" oops ");
    check_convert_fails<int>("1234abc");
    check_convert_fails<int>("1234 abc");
    check_convert_fails<int>("0.1");

    check_convert<uint16_t>("0", 0);
    check_convert<uint16_t> ("65535", 65535);
    check_convert_fails<uint16_t> ("-1");
    check_convert_fails<uint16_t> ("65536");

    check_convert<double>("1.23", 1.23);
    check_convert<double>("-1.23e-5", -1.23e-5);
    check_convert<double>("-5", -5.0);
    check_convert<double>(".23", 0.23);
    check_convert<double>("-.034e3", -0.034e3);
    check_convert<double>("  0.03e20  ", 0.03e20);

    check_convert_fails<double>("");
    check_convert_fails<double>("  ");
    check_convert_fails<double>("oops");
    check_convert_fails<double>(" oops ");
    check_convert_fails<double>("5x");
    check_convert_fails<double>("-1.3e20x");
    
    check_convert<bool>("t",true);
    check_convert<bool>("TRUE",true);
    check_convert<bool>("F",false);
    check_convert<bool>("False",false);
    check_convert_fails<bool>("False2");

    cerr << "test_lexical_cast(): success\n";
}


}  // namespace rf_pipelines
