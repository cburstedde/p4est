cmake_minimum_required(VERSION 3.15)
project(p4est
  LANGUAGES C
  VERSION 0.0.0)

include(CTest)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/Modules)

include(cmake/compilers.cmake)
include(cmake/options.cmake)

if(mpi)
  find_package(MPI COMPONENTS C REQUIRED)
endif()

find_package(ZLIB)

# --- libsc

if(NOT sc_external)
  find_package(sc)
endif()
if(sc_external OR NOT sc_FOUND)
  include(cmake/sc.cmake)
endif()

# --- p4est

if(Ep4est)
add_library(p4est)
target_include_directories(p4est
  PRIVATE src ${PROJECT_BINARY_DIR}/include
  INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:include>)
target_link_libraries(p4est PRIVATE SC::SC)
if(MPI_FOUND)
  target_link_libraries(p4est PRIVATE MPI::MPI_C)
endif()
endif(Ep4est)

#--- p6est
if(Ep6est AND Ep4est AND Ep8est)
# TODO: check that p6est linked properly
add_library(p6est)
target_include_directories(p6est
  PRIVATE src ${PROJECT_BINARY_DIR}/include
  INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:include)
target_link_libraries(p6est PRIVATE SC::SC)
if(MPI_FOUND)
  target_link_libraries(p6est PRIVATE MPI::MPI_C)
endif()
endif()

#--- p8est
if(Ep8est)
add_library(p8est)
target_include_directories(p8est
  PRIVATE src ${PROJECT_BINARY_DIR}/include
  INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:include)
target_link_libraries(p8est PRIVATE SC::SC)
if(MPI_FOUND)
  target_link_libraries(p8est PRIVATE MPI::MPI_C)
endif()
endif(Ep8est)

add_subdirectory(src)

if(BUILD_TESTING)
  add_subdirectory(test)
endif()
