#include <fstream>
#include <sstream>
#include "rf_pipelines_internals.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif


int wi_transform::add_plot_group(const string &name, int nt_per_pix, int ny)
{
    // Note: argument checking is done in plot_group constructor
    this->plot_groups.push_back(make_shared<plot_group>(name, nt_per_pix, ny));
    return plot_groups.size()-1;
}


string wi_transform::add_plot(const string &basename, int64_t it0, int nt, int nx, int ny, int group_id)
{
    if (plot_groups.size() == 0)
	throw runtime_error("wi_transform::add_plot(): no plot groups defined, maybe you forgot to call add_plot_group()?");

    if ((group_id < 0) || (group_id >= (int)plot_groups.size()))
	throw runtime_error("wi_transform::add_plot(): bad group_id specified");

    shared_ptr<plot_group> g = plot_groups[group_id];

    if (nt != g->nt_per_pix * nx)
	throw runtime_error("wi_transform::add_plot(): requirement (nt == nx*nt_per_pix) failed");
    if (ny != g->ny)
	throw runtime_error("wi_transform::add_plot(): ny doesn't match value specified in add_plot_group()");

    if (g->is_empty) {
	g->is_empty = false;
	g->curr_it0 = it0;
    }
    else if (it0 != g->curr_it1)
	throw runtime_error("wi_transform::add_plot(): plot time ranges are not contiguous");

    string filename = this->add_file(basename);

    Json::Value file;
    file["filename"] = basename;
    file["it0"] = Json::Value::Int64(it0);
    file["nx"] = nx;

    g->curr_it1 = it0 + nt;
    g->files.append(file);
    return filename;
}


string wi_transform::add_file(const string &basename)
{
    if (!outdir_manager)
	throw runtime_error("rf_pipelines: internal error: no outdir_manager in wi_transform::add_file()");
    if (outdir_manager->outdir.size() == 0)
	throw runtime_error("rf_pipelines: transform '" + this->name + "' attempted to write output file, but outdir=None was specified in the stream constructor");

    return this->outdir_manager->add_file(basename);
}


static void merge_json(Json::Value &dst, const Json::Value &src)
{
    if (src.isNull())
	return;

    if (!src.isObject())
	throw runtime_error("expected json value to be of 'object' type");

    for (const auto &key : src.getMemberNames())
        dst[key] = src[key];
}


// default implementation of virtual 
void wi_transform::_get_json(Json::Value &dst) const
{
    merge_json(dst, this->json_persistent);
    merge_json(dst, this->json_per_stream);
    merge_json(dst, this->json_per_substream);
}


// default implementation of virtual 
void wi_transform::_clear_json(bool substream_only)
{
    this->json_per_substream.clear();
    
    if (!substream_only)
	this->json_per_stream.clear();
}


// Default implementation of virtual
Json::Value wi_transform::serialize_to_json() const
{
    throw runtime_error("rf_pipelines: serialization-to-json not implemeted for transform '" + this->name + "'");
}


// -------------------------------------------------------------------------------------------------
//
// More json serialization/deserialization


// Helper for deserialize_transform_from_json()
static string _get_string(const Json::Value &x, const string &k)
{
    if (!x.isMember(k))
	throw runtime_error("rf_pipelines: deserialize_transform_from_json(): member '" + k + "' was expected but not found");
    
    const Json::Value &v = x[k];
    if (!v.isString())
	throw runtime_error("rf_pipelines: deserialize_transform_from_json(): member '" + k + "' was not a string as expected");

    return v.asString();
}

// Helper for deserialize_transform_from_json()
static int _get_int(const Json::Value &x, const string &k)
{
    if (!x.isMember(k))
	throw runtime_error("rf_pipelines: deserialize_transform_from_json(): member '" + k + "' was expected but not found");
    
    const Json::Value &v = x[k];
    if (!v.isInt())
	throw runtime_error("rf_pipelines: deserialize_transform_from_json(): member '" + k + "' was not an int as expected");

    return v.asInt();
}

// Helper for deserialize_transform_from_json()
static double _get_double(const Json::Value &x, const string &k)
{
    if (!x.isMember(k))
	throw runtime_error("rf_pipelines: deserialize_transform_from_json(): member '" + k + "' was expected but not found");
    
    const Json::Value &v = x[k];
    if (!v.isDouble())
	throw runtime_error("rf_pipelines: deserialize_transform_from_json(): member '" + k + "' was not a floating-point number as expected");

    return v.asDouble();
}

// Helper for deserialize_transform_from_json()
static axis_type _get_axis(const Json::Value &x, const string &k)
{
    string s = _get_string(x, k);

    if (s == "AXIS_FREQ")
	return AXIS_FREQ;
    if (s == "AXIS_TIME")
	return AXIS_TIME;
    if (s == "AXIS_NONE")
	return AXIS_NONE;

    throw runtime_error("rf_pipelines: deserialize_transform_from_json(): member '" + k + "' was not an axis_type as expected");
}


shared_ptr<wi_transform> deserialize_transform_from_json(const Json::Value &x)
{
    if (!x.isObject())
	throw runtime_error("rf_pipelines: deserialize_transform_from_json(): argument is not a json object as expected");

    string transform_name = _get_string(x, "transform_name");

    if (transform_name == "polynomial_detrender") {
	return make_polynomial_detrender(_get_int(x, "nt_chunk"),
					 _get_axis(x, "axis"),
					 _get_int(x, "polydeg"),
					 _get_double(x, "epsilon"));
    }

    throw runtime_error("rf_pipelines::deserialize_transform_from_json(): transform_name='" + transform_name + "' not recognized");
}


}  // namespace rf_pipelines
