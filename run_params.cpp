#include "rf_pipelines_internals.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}  // emacs pacifier
#endif


template<typename T>
static void _mismatch_helper(string &s, const char *k, const T &v1, const T &v2)
{
    if (v1 == v2)
	return;
    if (s.size() > 0)
	s.append(", ");
    s.append(k);
}


string run_params::mismatch(const run_params &p) const
{
    string ret;
    
    _mismatch_helper(ret, "outdir", outdir, p.outdir);
    _mismatch_helper(ret, "clobber", clobber, p.clobber);
    _mismatch_helper(ret, "img_nzoom", img_nzoom, p.img_nzoom);
    _mismatch_helper(ret, "img_nds", img_nds, p.img_nds);
    _mismatch_helper(ret, "img_nx", img_nx, p.img_nx);
    _mismatch_helper(ret, "verbosity", verbosity, p.verbosity);
    _mismatch_helper(ret, "debug", debug, p.debug);

    return ret;
}


void run_params::check() const
{
    if (extra_attrs.isObject())
	throw runtime_error("rf_pipelines: run_params::extra_attrs must be a Json::Object");

    if ((img_nzoom < 1) || (img_nzoom > 10))
	throw runtime_error("rf_pipelines: expected img_nzoom(=" + to_string(img_nzoom) + ") to be between 1 and 10");
    if ((img_nds <= 0) || !is_power_of_two(img_nds))
	throw runtime_error("rf_pipelines: img_nds(=" + to_string(img_nds) + ") must be a power of two");

    // Note: the upsampling logic in zoomable_tileset_state::_emit_plot() assumes img_nx is even.
    if ((img_nx <= 0) || (img_nx % 2))
	throw runtime_error("rf_pipelines: img_nx(=" + to_string(img_nx) + ") must be positive and even");
}


}  // namespace rf_pipelines
