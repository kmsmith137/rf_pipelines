// FIXME exception text could be improved

#include <cstring>
#include <sstream>
#include <fstream>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <unistd.h>

#include "rf_pipelines_internals.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
};  // pacify emacs c-mode!
#endif


string json_stringify(const Json::Value &j)
{
    Json::StyledWriter w;
    return w.write(j);
}


Json::Value json_read(const std::string &filename, bool noisy)
{
    ifstream f(filename);

    if (f.fail())
	throw runtime_error(filename + ": couldn't open file");

    Json::Value j;

    try {
	f >> j;
    } catch (...) {
	throw runtime_error(filename + ": couldn't parse json file");
    }

    if (noisy)
	cout << "read " << filename << endl;

    return j;
}


void json_write(const string &filename, const Json::Value &j, bool noisy)
{
    ofstream f(filename);
	
    if (f.fail())
	throw runtime_error(filename + ": couldn't open file for writing");
    
    f << json_stringify(j);
    
    if (noisy)
	cout << "wrote " << filename << endl;
}


static const Json::Value &get_member(const Json::Value &j, const string &k)
{
    if (!j.isObject())
	throw runtime_error("rf_pipelines: json value was not an Object as expected");
    if (!j.isMember(k))
	throw runtime_error("rf_pipelines: json field '" + k + "' was expected but not found");
    return j[k];
}


string string_from_json(const Json::Value &j, const string &k)
{
    const Json::Value &v = get_member(j, k);

    if (!v.isString())
	throw runtime_error("rf_pipelines: json field '" + k + "' was not a string as expected");

    return v.asString();
}


int int_from_json(const Json::Value &j, const string &k)
{
    const Json::Value &v = get_member(j, k);

    if (!v.isIntegral())
	throw runtime_error("rf_pipelines: json field '" + k + "' was not an integer as expected");

    return v.asInt();
}


bool bool_from_json(const Json::Value &j, const string &k)
{
    const Json::Value &v = get_member(j, k);

    if (!v.isBool())
	throw runtime_error("rf_pipelines: json field '" + k + "' was not boolean as expected");

    return v.asBool();
}


ssize_t ssize_t_from_json(const Json::Value &j, const string &k)
{
    const Json::Value &v = get_member(j, k);

    if (!v.isIntegral())
	throw runtime_error("rf_pipelines: json field '" + k + "' was not an integer as expected");

    return v.asInt64();
}


uint64_t uint64_t_from_json(const Json::Value &j, const string &k)
{
    const Json::Value &v = get_member(j, k);

    if (!v.isIntegral())
	throw runtime_error("rf_pipelines: json field '" + k + "' was not an integer as expected");

    return v.asUInt64();
}


double double_from_json(const Json::Value &j, const string &k)
{
    const Json::Value &v = get_member(j, k);

    if (!v.isDouble())
	throw runtime_error("rf_pipelines: json field '" + k + "' was not a floating-point number as expected");

    return v.asDouble();
}


rf_kernels::axis_type axis_type_from_json(const Json::Value &j, const string &k)
{
    string s = string_from_json(j, k);
    return rf_kernels::axis_type_from_string(s, "json field");
}


Json::Value array_from_json(const Json::Value &j, const string &k)
{
    const Json::Value &v = get_member(j, k);

    if (!v.isArray())
	throw runtime_error("rf_pipelines: json field '" + k + "' was not an array as expected");

    return v;
}


void add_json_object(Json::Value &dst, const Json::Value &src)
{
    if (src.isNull())
	return;

    if (!src.isObject())
	throw runtime_error("rf_pipelines internal error: 'src' argument to add_json_object() was not an Object as expected");

    for (const auto &key: src.getMemberNames())
	dst[key] = src[key];
}


}  // namespace rf_pipelines
