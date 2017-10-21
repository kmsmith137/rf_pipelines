#include "rf_pipelines_internals.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}  // emacs pacifier
#endif


void run_params::check() const
{
    if ((img_nzoom < 1) || (img_nzoom > 10))
	throw runtime_error("expected img_nzoom(=" + to_string(img_nzoom) + ") to be between 1 and 10");
    if ((img_nds <= 0) || !is_power_of_two(img_nds))
	throw runtime_error("img_nds(=" + to_string(img_nds) + ") must be a power of two");

    // Note: the upsampleing logic in zoomable_tileset_state::_emit_plot() assumes img_nx is even.
    if ((img_nx <= 0) || (img_nx % 2))
	throw runtime_error("img_nx(=" + to_string(img_nx) + ") must be positive and even");
}


}  // namespace rf_pipelines
