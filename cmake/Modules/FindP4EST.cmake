# FindP4EST.cmake
# ---------------
#
# Find p4est library
#
# Result Variables
# ----------------
#
# This module defines the following variables::
#
#   P4EST_FOUND
#   P4EST_INCLUDE_DIRS   - include directories for p4est
#   P4EST_LIBRARIES      - link against this library to use p4est
#
# Imported Targets
# ^^^^^^^^^^^^^^^^
#   P4EST::P4EST


find_path (P4EST_INCLUDE_DIR
  NAMES p4est.h
  DOC "p4est header")

find_library (P4EST_LIBRARY
  NAMES p4est
  DOC "p4est library")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args (P4EST
  REQUIRED_VARS P4EST_LIBRARY P4EST_INCLUDE_DIR
  )

if(P4EST_FOUND)

set(P4EST_INCLUDE_DIRS ${P4EST_INCLUDE_DIR})
set(P4EST_LIBRARIES ${P4EST_LIBRARY})

if(NOT TARGET P4EST::P4EST)
    add_library(P4EST::P4EST INTERFACE IMPORTED)
    set_property(TARGET P4EST::P4EST PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${P4EST_INCLUDE_DIR}")
    set_property(TARGET P4EST::P4EST PROPERTY INTERFACE_LINK_LIBRARIES "${P4EST_LIBRARY}")
endif()

endif()

mark_as_advanced(P4EST_INCLUDE_DIR P4EST_LIBRARY)
