include(FetchContent)

set(sc_external true CACHE BOOL "build sc library" FORCE)

# provides imported target SC::SC

FetchContent_Declare(SC
  GIT_REPOSITORY https://github.com/scivision/libsc.git
  GIT_TAG prev3-develop
  CMAKE_ARGS -Dmpi=${mpi})

if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.14)
  FetchContent_MakeAvailable(SC)
elseif(NOT sc_POPULATED)
  FetchContent_Populate(SC)
  add_subdirectory(${sc_SOURCE_DIR} ${sc_BINARY_DIR})
endif()
