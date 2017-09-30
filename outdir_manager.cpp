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

    // By convention, the toplevel pipeline json file is named 'rf_pipeline_0.json'.
    string toplevel_filename = outdir + "rf_pipeline_0.json";
    
    if (!clobber_ok && file_exists(toplevel_filename))
	throw runtime_error("Directory \"" + outdir + "\" contains rf_pipeline_0.json, and 'clobber' flag wasn't set");
}


// Returns the full pathname
string outdir_manager::add_file(const string &basename)
{
    if (outdir.size() == 0)
	throw runtime_error("rf_pipelines: transform attempted to write output file, but outdir=None was specified in the stream constructor");

    // The return value of unordered_set::insert() is a std::pair whose 
    // second element is 'true' if the inserted value doesn't already exist.
    bool is_new = basenames.insert(basename).second;
    string ret = outdir + basename;

    if (!is_new)
	throw runtime_error("rf_pipelines: output file '" + ret + "' was written twice in same pipeline run");

    return ret;
}
    

}  // namespace rf_pipelines
