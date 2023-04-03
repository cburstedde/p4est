cmake_minimum_required(VERSION 3.13...3.24)

include(cmake/git.cmake)

project(P4EST
LANGUAGES C
HOMEPAGE_URL https://www.p4est.org/
DESCRIPTION "p4est manages a collection---a forest---of octrees in parallel."
VERSION ${PROJECT_VERSION})

include(CTest)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/Modules)

# --- user options

include(cmake/options.cmake)
message(STATUS "p4est ${PROJECT_VERSION} "
               "install prefix: ${CMAKE_INSTALL_PREFIX}")

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

target_link_libraries(SC::SC INTERFACE
$<$<BOOL:${MPI_C_FOUND}>:MPI::MPI_C>
$<$<BOOL:${ZLIB_FOUND}>:ZLIB::ZLIB>
$<$<BOOL:${P4EST_NEED_M}>:m>
)
# --- p4est
# p4est is always needed.

add_library(p4est)
set_property(TARGET p4est PROPERTY EXPORT_NAME P4EST)
target_include_directories(p4est PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:include>)
target_link_libraries(p4est PUBLIC SC::SC)

# imported target, for use from FetchContent
add_library(P4EST::P4EST INTERFACE IMPORTED GLOBAL)
target_link_libraries(P4EST::P4EST INTERFACE p4est)

install(DIRECTORY ${PROJECT_SOURCE_DIR}/src/ ${PROJECT_BINARY_DIR}/include/
  DESTINATION include
  FILES_MATCHING PATTERN "*.h")

install(TARGETS p4est
EXPORT ${PROJECT_NAME}-targets
ARCHIVE DESTINATION lib
LIBRARY DESTINATION lib
RUNTIME DESTINATION bin
)

#--- p8est

if(enable_p8est)
  add_library(p8est OBJECT)
  target_include_directories(p8est PRIVATE src ${PROJECT_BINARY_DIR}/include)
  target_link_libraries(p8est PRIVATE SC::SC)
  target_sources(p4est PRIVATE $<TARGET_OBJECTS:p8est>)
endif(enable_p8est)

#--- p6est

if(enable_p6est AND enable_p8est)
  add_library(p6est OBJECT)
  target_include_directories(p6est PRIVATE src ${PROJECT_BINARY_DIR}/include)
  target_link_libraries(p6est PRIVATE SC::SC)
  target_sources(p4est PRIVATE $<TARGET_OBJECTS:p6est>)
endif()

add_subdirectory(src)

# --- optional test and install

if(BUILD_TESTING)
  add_subdirectory(test)
endif(BUILD_TESTING)

# --- packaging

install(FILES
${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindSC.cmake
${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindP4EST.cmake
DESTINATION cmake)

include(cmake/pkgconf.cmake)
include(cmake/install.cmake)

include(FeatureSummary)
add_feature_info(MPI mpi "MPI features of ${PROJECT_NAME}")
add_feature_info(OpenMP openmp "OpenMP features of ${PROJECT_NAME}")
add_feature_info(P6EST enable_p6est "2D-3D p6est")
add_feature_info(P8EST enable_p8est "3D p8est")
feature_summary(WHAT ENABLED_FEATURES DISABLED_FEATURES)
