#include "pyclops.hpp"

using namespace pyclops;

PyMODINIT_FUNC initpyclops(void)
{
    import_array();

    extension_module m("pyclops", "Some hacks for writing hybrid C++/python code");
    m.finalize();
}

