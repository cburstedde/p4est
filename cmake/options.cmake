option(enable_p6est "build p6est" on)
option(enable_p8est "build p8est" on)

option(vtk "VTK interface" on)

option(sc_external "force build of libsc")

set(CMAKE_EXPORT_COMPILE_COMMANDS on)


# --- default install directory under build/local
# users can specify like "cmake -B build -DCMAKE_INSTALL_PREFIX=~/mydir"
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  # will not take effect without FORCE
  set(CMAKE_INSTALL_PREFIX "${PROJECT_BINARY_DIR}/local" CACHE PATH "Install top-level directory" FORCE)
endif()
