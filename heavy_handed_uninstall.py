#!/usr/bin/env python

import build_helpers

# If called recursively in superbuild, a global persistent HeavyHandedUninstaller will be returned.
u = build_helpers.get_global_heavy_handed_uninstaller()

u.uninstall_headers('rf_pipelines*.hpp')
u.uninstall_headers('rf_pipelines/*.hpp')
u.uninstall_headers('rf_pipelines/')
u.uninstall_libraries('librf_pipelines*')

u.uninstall_executables('rfp-*')

u.uninstall_python_package('rf_pipelines/streams')
u.uninstall_python_package('rf_pipelines/transforms')
u.uninstall_python_package('rf_pipelines/retirement_home')
u.uninstall_python_package('rf_pipelines/streams')
u.uninstall_pyfiles('rf_pipelines*.py')
u.uninstall_pyfiles('rf_pipelines*.pyc')
u.uninstall_pyfiles('rf_pipelines*.so')

# FIXME currently hacking in pyclops as a git subtree!
# This doesn't seem like the right thing to do, but I'm planning to phase out pyclops soon anyway.
u.uninstall_headers('pyclops/*.hpp')
u.uninstall_headers('pyclops.hpp')
u.uninstall_headers('pyclops/')
u.uninstall_libraries('libpyclops*')
u.uninstall_pyfiles('pyclops.so')

# Files which no longer exist, but did exist in previous versions of rf_pipelines
u.uninstall_headers('chime_packetizer.hpp')
u.uninstall_headers('chime_file_stream_base.hpp')
u.uninstall_headers('reverter.hpp')

# If called recursively in superbuild, run() will not be called here.
if __name__ == '__main__':
    u.run()
