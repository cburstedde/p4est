# Findsc.cmake
# ---------------
#
# Find sc library
#
# Result Variables
# ----------------
#
# This module defines the following variables::
#
#   sc_FOUND
#   sc_INCLUDE_DIRS   - include directories for sc
#   sc_LIBRARIES      - link against this library to use sc
#
# Imported Targets
# ^^^^^^^^^^^^^^^^
#   SC::SC


find_path (sc_INCLUDE_DIR
  NAMES sc.h
  DOC "sc header")

find_library (sc_LIBRARY
  NAMES sc
  DOC "sc library")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args (sc
  REQUIRED_VARS sc_LIBRARY sc_INCLUDE_DIR)

if (sc_FOUND)
  set(sc_INCLUDE_DIRS ${sc_INCLUDE_DIR})
endif()

if(sc_FOUND AND NOT TARGET SC::SC)
    add_library(SC::SC INTERFACE IMPORTED)

    set_target_properties(SC::SC PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${sc_INCLUDE_DIRS}"
        INTERFACE_LINK_LIBRARIES "${sc_LIBRARY}"
    )
endif()

mark_as_advanced (sc_INCLUDE_DIR sc_LIBRARY)
