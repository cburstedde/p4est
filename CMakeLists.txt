cmake_minimum_required(VERSION 3.13...3.20)

include(cmake/git.cmake)

project(P4EST
LANGUAGES C
HOMEPAGE_URL https://github.com/cburstedde/p4est
DESCRIPTION "The 'p4est' forest-of-octrees library"
VERSION ${PROJECT_VERSION})

message(STATUS "${PROJECT_NAME} ${PROJECT_VERSION}")

include(CTest)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/Modules)

# --- user options

include(cmake/options.cmake)


# --- external libs

if(mpi)
  find_package(MPI COMPONENTS C REQUIRED)
  if(NOT MPIEXEC_EXECUTABLE)
    message(FATAL_ERROR "MPIEXEC is missing, MPI will not work properly")
  endif()
endif()

find_package(ZLIB)

# --- libsc

if(NOT sc_external)
  find_package(SC)
endif()
if(NOT TARGET SC::SC)
  include(cmake/sc.cmake)
endif()

# --- configure p4est

include(cmake/config.cmake)
include(cmake/compilers.cmake)

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

install(DIRECTORY ${PROJECT_SOURCE_DIR}/src/ ${PROJECT_BINARY_DIR}/include/
  DESTINATION include
  FILES_MATCHING PATTERN "*.h")

install(TARGETS p4est
  EXPORT ${PROJECT_NAME}Targets
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)

#--- p8est

if(Ep8est)
add_library(p8est)
target_include_directories(p8est
  PRIVATE src ${PROJECT_BINARY_DIR}/include
  INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:include>)
target_link_libraries(p8est PRIVATE p4est SC::SC)
if(MPI_FOUND)
  target_link_libraries(p8est PRIVATE MPI::MPI_C)
endif()

# imported target, for use from FetchContent
add_library(P4EST::P8EST INTERFACE IMPORTED GLOBAL)
target_link_libraries(P4EST::P8EST INTERFACE p8est)

install(TARGETS p8est
  EXPORT ${PROJECT_NAME}Targets
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)

endif(Ep8est)

#--- p6est

if(Ep6est AND Ep8est)
add_library(p6est)
target_include_directories(p6est
  PRIVATE src ${PROJECT_BINARY_DIR}/include
  INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:include>)
target_link_libraries(p6est PRIVATE p8est p4est SC::SC)
if(MPI_FOUND)
  target_link_libraries(p6est PRIVATE MPI::MPI_C)
endif()

# imported target, for use from FetchContent
add_library(P4EST::P6EST INTERFACE IMPORTED GLOBAL)
target_link_libraries(P4EST::P6EST INTERFACE p6est)

install(TARGETS p6est
  EXPORT ${PROJECT_NAME}Targets
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)

endif()

add_subdirectory(src)

# --- optional test and install

if(BUILD_TESTING)
  add_subdirectory(test)
endif(BUILD_TESTING)

# --- packaging

include(cmake/install.cmake)
