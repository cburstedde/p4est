# provides imported target SC::SC
include(ExternalProject)
include(${CMAKE_CURRENT_LIST_DIR}/GitSubmodule.cmake)

set(sc_external true CACHE BOOL "build sc library" FORCE)

git_submodule("${PROJECT_SOURCE_DIR}/sc")

# --- libsc externalProject
# this keeps libsc scope totally separate from p4est, which avoids
# tricky to diagnose behaviors

if(NOT SC_ROOT)
  set(SC_ROOT ${CMAKE_INSTALL_PREFIX})
endif()

if(BUILD_SHARED_LIBS)
  set(SC_LIBRARIES ${SC_ROOT}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}sc${CMAKE_SHARED_LIBRARY_SUFFIX})
else()
  set(SC_LIBRARIES ${SC_ROOT}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}sc${CMAKE_STATIC_LIBRARY_SUFFIX})
endif()

set(SC_INCLUDE_DIRS ${SC_ROOT}/include)

ExternalProject_Add(SC
SOURCE_DIR ${PROJECT_SOURCE_DIR}/sc
CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${SC_ROOT} -Dmpi:BOOL=${mpi} -Dopenmp:BOOL=${openmp}
BUILD_BYPRODUCTS ${SC_LIBRARIES}
)

# --- imported target

file(MAKE_DIRECTORY ${SC_INCLUDE_DIRS})
# avoid race condition

# this GLOBAL is required to be visible via other
# project's FetchContent of this project
add_library(SC::SC INTERFACE IMPORTED GLOBAL)
target_include_directories(SC::SC INTERFACE "${SC_INCLUDE_DIRS}")
target_link_libraries(SC::SC INTERFACE "${SC_LIBRARIES}")

add_dependencies(SC::SC SC)
