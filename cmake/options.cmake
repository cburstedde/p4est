option(enable_p6est "build p6est" on)
option(enable_p8est "build p8est" on)

option(vtk "VTK interface" on)

option(mpi "use MPI library" off)
option(openmp "use OpenMP" off)

option(sc_external "force build of libsc")

set(CMAKE_EXPORT_COMPILE_COMMANDS on)
set(CMAKE_TLS_VERIFY on)

if(dev)

else()
  set_directory_properties(PROPERTIES EP_UPDATE_DISCONNECTED true)
endif()

# --- default install directory under build/local
# users can specify like "cmake -B build -DCMAKE_INSTALL_PREFIX=~/mydir"
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  # will not take effect without FORCE
  set(CMAKE_INSTALL_PREFIX "${PROJECT_BINARY_DIR}/local" CACHE PATH "Install top-level directory" FORCE)
endif()

# --- auto-ignore build directory
if(NOT EXISTS ${PROJECT_BINARY_DIR}/.gitignore)
  file(WRITE ${PROJECT_BINARY_DIR}/.gitignore "*")
endif()
