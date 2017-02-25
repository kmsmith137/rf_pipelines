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


}  // namespace rf_pipelines
