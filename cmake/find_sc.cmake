# The purpose of this cmake macro is to detect libsc from environnment
# and make target SC:SC available
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
      add_library(SC::SC ALIAS PkgConfig::P4EST_NEEDS_LIBSC)
      set(SC_FOUND True)
      set(SC_FOUND_VIA_PKGCONFIG True)
    endif()
  else()
    message(NOTICE "pkg-config executable is not available.")
  endif()

endif()

if(NOT SC_FOUND)
  message(STATUS "libsc was not found in environment, we will try to build it.")
endif()
