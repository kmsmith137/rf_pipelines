#include "rf_pipelines_internals.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif


wi_transform::wi_transform() :
    nt_chunk(0), nt_prepad(0), nt_postpad(0)
{ }


wi_transform::wi_transform(int nt_chunk_, int nt_prepad_, int nt_postpad_) :
    nt_chunk(nt_chunk_), nt_prepad(nt_prepad_), nt_postpad(nt_postpad_)
{
    this->check_invariants();
}


void wi_transform::check_invariants() const
{
    if (nt_chunk <= 0)
	throw runtime_error("wi_transform: nt_chunk is non-positive or uninitialized");
    if (nt_prepad < 0)
	throw runtime_error("wi_transform: nt_prepad is negative");
    if (nt_postpad < 0)
	throw runtime_error("wi_transform: nt_postpad is negative");
}


}  // namespace rf_pipelines
