#include <fstream>
#include <sstream>
#include "rf_pipelines_internals.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif

    
outdir_manager::outdir_manager(const string &outdir_, bool clobber_ok_) :
    outdir(outdir_),
    clobber_ok(clobber_ok_)
{
    if (outdir.size() == 0)
	return;

    // Ensure trailing slash
    if (outdir.back() != '/')
	outdir = outdir + "/";

    makedirs(outdir);
	
    if (clobber_ok)
	return;

    for (const string &s: listdir(outdir)) {
	if (is_json_basename(s))
	    throw runtime_error("Directory \"" + outdir + "\" contains stray files of the form rf_pipeline_NN.json, and 'clobber' flag wasn't set");
    }
}


// Returns the full pathname
string outdir_manager::add_file(const string &basename)
{
    if (outdir.size() == 0)
	throw runtime_error("rf_pipelines: transform attempted to write output file, but outdir=None was specified in the stream constructor");

    bool is_new = basename_set.insert(basename).second;
    string ret = outdir + basename;

    if (!is_new)
	throw runtime_error("rf_pipelines: output file '" + ret + "' was written twice in same pipeline run");

    return ret;
}


void outdir_manager::write_per_substream_json_file(int isubstream, const Json::Value &data, bool noisy)
{
    if (isubstream < 0)
	throw runtime_error("outdir_manager::write_json_file(): 'isubstream' arg was negative");

    stringstream ss;
    ss << this->outdir << "rf_pipeline_" << isubstream << ".json";

    string filename = ss.str();
    ofstream f(filename);

    if (f.fail())
	throw runtime_error("couldn't open output file " + filename);

    Json::StyledWriter w;
    f << w.write(data);

    if (noisy)
	cerr << ("wrote " + filename + "\n");
}


// Static member function
bool outdir_manager::is_json_basename(const string &basename)
{
    if (!startswith(basename, "rf_pipeline_"))
	return false;
    if (!endswith(basename, ".json"))
	return false;
    
    const char *s = basename.c_str();
    int len = strlen(s);
    
    if (len <= 14)
	return false;
    
    for (int i = 9; i < len-5; i++)
	if (!isdigit(s[i]))
	    return false;

    return true;
}
    

}  // namespace rf_pipelines
