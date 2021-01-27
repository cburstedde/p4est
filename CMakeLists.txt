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
# p4est is always needed.
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

# imported target, for use from FetchContent
add_library(P4EST::P4EST INTERFACE IMPORTED GLOBAL)
target_link_libraries(P4EST::P4EST INTERFACE p4est)

#--- p8est
if(Ep8est)
add_library(p8est)
target_include_directories(p8est
  PRIVATE src ${PROJECT_BINARY_DIR}/include
  INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:include)
target_link_libraries(p8est PRIVATE p4est SC::SC)
if(MPI_FOUND)
  target_link_libraries(p8est PRIVATE MPI::MPI_C)
endif()

# imported target, for use from FetchContent
add_library(P4EST::P8EST INTERFACE IMPORTED GLOBAL)
target_link_libraries(P4EST::P8EST INTERFACE p8est)

endif(Ep8est)

#--- p6est
if(Ep6est AND Ep8est)
add_library(p6est)
target_include_directories(p6est
  PRIVATE src ${PROJECT_BINARY_DIR}/include
  INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:include)
target_link_libraries(p6est PRIVATE p8est p4est SC::SC)
if(MPI_FOUND)
  target_link_libraries(p6est PRIVATE MPI::MPI_C)
endif()

# imported target, for use from FetchContent
add_library(P4EST::P6EST INTERFACE IMPORTED GLOBAL)
target_link_libraries(P4EST::P6EST INTERFACE p6est)

endif()

add_subdirectory(src)

if(BUILD_TESTING)
  add_subdirectory(test)
endif()
