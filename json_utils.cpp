// FIXME exception text could be improved

#include <cstring>
#include <sstream>
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

    if (!v.isInt())
	throw runtime_error("rf_pipelines: json field '" + k + "' was not an int as expected");

    return v.asInt();
}


bool bool_from_json(const Json::Value &j, const string &k)
{
    const Json::Value &v = get_member(j, k);

    if (!v.isBool())
	throw runtime_error("rf_pipelines: json field '" + k + "' was not boolean as expected");

    return v.asBool();
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


}  // namespace rf_pipelines
