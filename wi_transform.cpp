#include "rf_pipelines_internals.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif


wi_transform::wi_transform()
{
    this->reset_nt();
}


wi_transform::wi_transform(int nt_chunk_, int nt_prepad_, int nt_postpad_) :
    wi_transform()
{
    this->set_nt(nt_chunk_, nt_prepad_, nt_postpad_);
}


void wi_transform::set_nt(int nt_chunk_, int nt_prepad_, int nt_postpad_)
{
    this->nt_chunk = nt_chunk_;
    this->nt_prepad = nt_prepad_;
    this->nt_postpad = nt_postpad_;

    if (nt_chunk <= 0)
	throw runtime_error("wi_transform::set_nt(): expected nt_chunk > 0");
    if (nt_prepad < 0)
	throw runtime_error("wi_transform::set_nt(): expected nt_prepad >= 0");
    if (nt_postpad < 0)
	throw runtime_error("wi_transform::set_nt(): expected nt_postpad >= 0");    
}


void wi_transform::reset_nt()
{
    this->nt_chunk = 0;
    this->nt_prepad = 0;
    this->nt_postpad = 0;
}


}  // namespace rf_pipelines
