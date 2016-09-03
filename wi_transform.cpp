#include <fstream>
#include <sstream>
#include "rf_pipelines_internals.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif

// Not much in this source file right now, but I think it will grow with time!

string wi_transform::add_file(const string &basename)
{
    return this->outdir_manager->add_file(basename);
}
    

}  // namespace rf_pipelines
