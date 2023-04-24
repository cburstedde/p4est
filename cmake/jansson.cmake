# optional json library :
# 1. try cmake find_package in config mode
# 2. try pkg-config

find_package(jansson CONFIG)

if(jansson_FOUND)

  message(STATUS "jansson library found via find_package")

else()

  # make sure cmake macro pkg_check_modules is available
  find_package(PkgConfig)

  if (PKG_CONFIG_FOUND)
    pkg_check_modules(LIBSC_JANSSON QUIET IMPORTED_TARGET jansson)

    if (LIBSC_JANSSON_FOUND)
      message(STATUS "jansson library found via pkg-config")
      add_library(jansson::jansson ALIAS PkgConfig::LIBSC_JANSSON)
      set(jansson_FOUND 1)
    else()
      set(jansson_FOUND 0)
    endif()
  else()
    message(NOTICE "pkg-config executable is not available.")
    set(jansson_FOUND 0)
  endif()

endif()

if( NOT jansson_FOUND )
  message(NOTICE "libjansson was not found")
endif()

if (jansson_FOUND)
  set(SC_HAVE_JSON 1)
else()
  set(SC_HAVE_JSON 0)
endif()
