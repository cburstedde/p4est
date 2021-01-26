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

if(P4EST_LIBRARY AND P4EST_INCLUDE_DIR)
  set(P4EST_P4EST_FOUND true)
endif()

if(P8EST IN_LIST P4EST_FIND_COMPONENTS)
  find_library (P8EST_LIBRARY
    NAMES p8est
    DOC "p8est library")

  if(P8EST_LIBRARY)
    set(P4EST_P8EST_FOUND true)
  endif()
endif()

if(P6EST IN_LIST P4EST_FIND_COMPONENTS)
  find_library (P6EST_LIBRARY
    NAMES p6est
    DOC "p6est library")

  if(P6EST_LIBRARY)
    set(P4EST_P6EST_FOUND true)
  endif()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args (P4EST
  REQUIRED_VARS P4EST_LIBRARY P4EST_INCLUDE_DIR SC_LIBRARY SC_INCLUDE_DIR
  HANDLE_COMPONENTS)

if(P4EST_FOUND)

set(P4EST_INCLUDE_DIRS ${P4EST_INCLUDE_DIR})
set(P4EST_LIBRARIES ${P4EST_LIBRARY})

if(NOT TARGET P4EST::P4EST)
    add_library(P4EST::P4EST INTERFACE IMPORTED)
    set_target_properties(P4EST::P4EST PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${P4EST_INCLUDE_DIR}"
        INTERFACE_LINK_LIBRARIES "${P4EST_LIBRARY}"
    )
endif()

endif()

mark_as_advanced(P4EST_INCLUDE_DIR P4EST_LIBRARY)
