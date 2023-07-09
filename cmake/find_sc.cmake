# The purpose of this cmake macro is to detect libsc from environnment
# and make target SC::SC available
# We are using 3 types of detections:
# - find_package in CONFIG mode
# - find_package in MODULE mode
# - pkg-config (usefull for detecting libsc built with autotools)

# First try to detect libsc by searching for SCConfig.cmake in CMAKE_PREFIX_PATH
find_package(SC CONFIG QUIET)
if (SC_FOUND)
  message(STATUS "libsc found via SCConfig.cmake found in CMAKE_PREFIX_PATH")
endif()

# Second try to find libsc by searching for FindSC.cmake in CMAKE_MODULE_PATH
if (NOT SC_FOUND)
  find_package(SC QUIET)
  if (SC_FOUND)
    message(STATUS "libsc found via FindSC.cmake found in CMAKE_MODULE_PATH")
  endif()
endif()

# Third and at last resort, try pkg-config
if (NOT SC_FOUND)

  # make sure cmake macro pkg_check_modules is available
  find_package(PkgConfig)

  if (PKG_CONFIG_FOUND)
    pkg_check_modules(P4EST_NEEDS_LIBSC QUIET IMPORTED_TARGET libsc)

    if (P4EST_NEEDS_LIBSC_FOUND)
      message(STATUS "libsc found via pkg-config in PKG_CONFIG_PATH")

      # we can't make an alias library here, because it is read only, and thus can't be used
      # later in top-level CMakeLists.txt to interface with required dependencies
      #add_library(SC::SC ALIAS PkgConfig::P4EST_NEEDS_LIBSC)

      # so just re-create target SC::SC from scratch, by copying properties from pkgconfig
      get_target_property(SC_INTERFACE_INCLUDE_DIRECTORIES PkgConfig::P4EST_NEEDS_LIBSC INTERFACE_INCLUDE_DIRECTORIES)
      get_target_property(SC_INTERFACE_LINK_LIBRARIES PkgConfig::P4EST_NEEDS_LIBSC INTERFACE_LINK_LIBRARIES)
      add_library(SC::SC INTERFACE IMPORTED)
      set_property(TARGET SC::SC PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${SC_INTERFACE_INCLUDE_DIRECTORIES}")
      set_property(TARGET SC::SC PROPERTY INTERFACE_LINK_LIBRARIES "${SC_INTERFACE_LINK_LIBRARIES}")

      # does libsc has mpi support ?
      execute_process(COMMAND ${PKG_CONFIG_EXECUTABLE} --variable have_mpi libsc
        OUTPUT_VARIABLE SC_HAVE_MPI)

      # does libsc has json support ?
      execute_process(COMMAND ${PKG_CONFIG_EXECUTABLE} --variable have_json libsc
        OUTPUT_VARIABLE SC_HAVE_JSON)

      set(SC_FOUND True)
      set(SC_FOUND_VIA_PKGCONFIG True)
    endif()
  else()
    message(NOTICE "pkg-config executable is not available.")
  endif()

endif()

if(SC_HAVE_JSON OR SC_have_json)
  include(cmake/jansson.cmake)
  if(NOT jansson_FOUND)
    message(FATAL_ERROR "libsc requires jansson.")
  endif()
endif()

if (mpi AND NOT SC_HAVE_MPI)
  # notice find_package(MPI) already in top-level CMakeLists.txt, don't need to do it again here.

  # this is actually a fatal error, but could be made a warning (TODO)
  message(WARNING "You enabled MPI at p4est level, but libsc was not compiled with MPI support.")
endif()
