cmake_minimum_required(VERSION 3.13...3.20)

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
set(CMAKE_FIND_PACKAGE_PREFER_CONFIG true)

if(NOT sc_external)
  find_package(SC)
endif()
if(NOT TARGET SC::SC)
  include(cmake/sc.cmake)
endif()

if(ZLIB_FOUND)
  target_link_libraries(SC::SC INTERFACE ZLIB::ZLIB)
endif()

# --- configure p4est

include(cmake/config.cmake)
include(cmake/compilers.cmake)

# --- p4est
# p4est is always needed.

add_library(p4est)
set_target_properties(p4est PROPERTIES EXPORT_NAME P4EST)
target_include_directories(p4est
  PRIVATE src ${PROJECT_BINARY_DIR}/include
  INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:include>)
target_link_libraries(p4est PUBLIC SC::SC)
if(MPI_FOUND)
  target_link_libraries(p4est PUBLIC MPI::MPI_C)
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

if(enable_p8est)
add_library(p8est)
set_target_properties(p8est PROPERTIES EXPORT_NAME P8EST)
target_include_directories(p8est
  PRIVATE src ${PROJECT_BINARY_DIR}/include
  INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:include>)
target_link_libraries(p8est PUBLIC p4est)

# imported target, for use from FetchContent
add_library(P4EST::P8EST INTERFACE IMPORTED GLOBAL)
target_link_libraries(P4EST::P8EST INTERFACE p8est)

install(TARGETS p8est
  EXPORT ${PROJECT_NAME}Targets
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)

endif(enable_p8est)

#--- p6est

if(enable_p6est AND enable_p8est)
add_library(p6est)
set_target_properties(p6est PROPERTIES EXPORT_NAME P6EST)
target_include_directories(p6est
  PRIVATE src ${PROJECT_BINARY_DIR}/include
  INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:include>)
target_link_libraries(p6est PUBLIC p8est p4est)

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

install(FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/FindSC.cmake
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/FindP4EST.cmake
  DESTINATION lib/cmake/Modules)

include(cmake/pkgconf.cmake)
include(cmake/install.cmake)

include(FeatureSummary)
add_feature_info(MPI mpi "MPI features of ${PROJECT_NAME}")
add_feature_info(OpenMP openmp "OpenMP features of ${PROJECT_NAME}")
add_feature_info(P6EST enable_p6est "2D-3D p6est")
add_feature_info(P8EST enable_p8est "3D p8est")
feature_summary(WHAT ENABLED_FEATURES DISABLED_FEATURES)
