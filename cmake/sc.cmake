include(FetchContent)

set(sc_external true CACHE BOOL "build sc library" FORCE)

# provides imported target SC::SC

FetchContent_Declare(SC
  GIT_REPOSITORY https://github.com/scivision/libsc.git
  GIT_TAG prev3-develop
  CMAKE_ARGS -Dmpi=${mpi})
FetchContent_MakeAvailable(SC)
