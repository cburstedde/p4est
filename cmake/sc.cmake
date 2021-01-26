# provides imported target SC::SC

set(sc_external true CACHE BOOL "build sc library" FORCE)

if(NOT EXISTS sc/CMakeLists.txt)
  find_package(Git REQUIRED)

  execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init
    RESULT_VARIABLE _err)

  if(NOT _err EQUAL 0)
    message(SEND_ERROR "Could not get Git submodule for libsc.")
  endif()
endif()

add_subdirectory(sc)
