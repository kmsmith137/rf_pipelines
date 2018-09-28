#include "rf_pipelines_internals.hpp"
#include "rf_pipelines_inventory.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif




std::vector<rf_pipelines::mask_measurements>
mask_measurements_ringbuf::get_all_measurements() const {
    std::vector<rf_pipelines::mask_measurements> rtn;
    throw runtime_error("unimplemented!");
    return rtn;
}

std::unordered_map<std::string, float> 
mask_measurements_ringbuf::get_stats(float period) const {
    std::unordered_map<std::string, float> rtn;
    throw runtime_error("unimplemented!");
    return rtn;
}



}  // namespace rf_pipelines
