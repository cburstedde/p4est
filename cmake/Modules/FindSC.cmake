# FindSC.cmake
# ---------------
#
# Find libsc library
#
# Result Variables
# ----------------
#
# This module defines the following variables::
#
#   SC_FOUND
#   SC_INCLUDE_DIRS   - include directories for libsc
#   SC_LIBRARIES      - link against this library to use libsc
#
# Imported Targets
# ^^^^^^^^^^^^^^^^
#   SC::SC


find_path (SC_INCLUDE_DIR
  NAMES sc.h
  DOC "libsc header")

find_library (SC_LIBRARY
  NAMES sc
  DOC "libsc library")


if(SC_INCLUDE_DIR AND SC_LIBRARY)
  set(SC_mpi_ok true)

  # check if libsc was configured in compatible way
  include(CheckSymbolExists)
  set(CMAKE_REQUIRED_FLAGS)
  set(CMAKE_REQUIRED_INCLUDES)
  set(CMAKE_REQUIRED_LIBRARIES)

  # libsc and current project must both be compiled with/without MPI
  check_symbol_exists(SC_ENABLE_MPI ${SC_INCLUDE_DIR}/sc_config.h SC_has_mpi)
  check_symbol_exists(SC_ENABLE_MPIIO ${SC_INCLUDE_DIR}/sc_config.h SC_has_mpi_io)

  if(MPI_C_FOUND)
    # a sign the current project is using MPI
    if(NOT (SC_has_mpi AND SC_has_mpi_io))
      set(SC_mpi_ok false)
    endif()
  else()
    if(SC_has_mpi OR SC_has_mpi_io)
      set(SC_mpi_ok false)
    endif()
  endif()

endif()


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args (SC
  REQUIRED_VARS SC_LIBRARY SC_INCLUDE_DIR SC_mpi_ok)

if(SC_FOUND)

set(SC_INCLUDE_DIRS ${SC_INCLUDE_DIR})
set(SC_LIBRARIES ${SC_LIBRARY})

if(NOT TARGET SC::SC)
  add_library(SC::SC INTERFACE IMPORTED)

  set_target_properties(SC::SC PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${SC_INCLUDE_DIR}"
    INTERFACE_LINK_LIBRARIES "${SC_LIBRARY}")
endif()

endif(SC_FOUND)

mark_as_advanced(SC_INCLUDE_DIR SC_LIBRARY)
