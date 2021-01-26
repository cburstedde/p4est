cmake_minimum_required(VERSION 3.13)
project(p4estExamples LANGUAGES C)

include(CTest)
include(CheckIncludeFile)
include(CheckSymbolExists)

# --- check system capabilities
check_include_file(arpa/inet.h P4EST_HAVE_ARPA_INET_H)
check_include_file(netinet/in.h P4EST_HAVE_NETINET_IN_H)
if(WIN32 AND NOT P4EST_HAVE_ARPA_INET_H AND NOT P4EST_HAVE_NETINET_IN_H)
  check_include_file(Winsock2.h P4EST_HAVE_WINSOCK2_H)
  set(WINSOCK_LIBRARIES wsock32 ws2_32) # Iphlpapi
endif()

if(NOT (P4EST_HAVE_ARPA_INET_H OR P4EST_HAVE_NETINET_IN_H OR P4EST_HAVE_WINSOCK2_H))
  return()
endif()

check_include_file(unistd.h P4EST_HAVE_UNISTD_H)
if(P4EST_HAVE_UNISTD_H)
  check_include_file(getopt.h P4EST_HAVE_GETOPT_H)
endif()


list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../cmake/Modules)

# --- configure external libraries
find_package(MPI COMPONENTS C)
find_package(ZLIB REQUIRED)

find_package(P4EST CONFIG)

if(P4EST_FOUND)
  # include in this order

  if(EXISTS ${P4EST_DIR}/../SC/SCTargets.cmake)
    include(${P4EST_DIR}/../SC/SCTargets.cmake)

    add_library(SC::SC INTERFACE IMPORTED)
    target_link_libraries(SC::SC INTERFACE SC::sc)
  endif()

  include(${P4EST_DIR}/P4ESTTargets.cmake)
  add_library(P4EST::P4EST INTERFACE IMPORTED)
  target_link_libraries(P4EST::P4EST INTERFACE P4EST::p4est)

  # see if p8est or p6est available
  if(TARGET P4EST::p6est)
    add_library(P4EST::P6EST INTERFACE IMPORTED)
    target_link_libraries(P4EST::P6EST INTERFACE P4EST::p6est)
  endif()

  if(TARGET P4EST::p8est)
    add_library(P4EST::P8EST INTERFACE IMPORTED)
    target_link_libraries(P4EST::P8EST INTERFACE P4EST::p8est)
  endif()

else()
  find_package(P4EST REQUIRED)
  find_package(SC REQUIRED)
endif()


cmake_host_system_information(RESULT Ncpu QUERY NUMBER_OF_PHYSICAL_CORES)
if(Ncpu LESS 2)
  ProcessorCount(n)
  if(n GREATER Ncpu)
    set(Ncpu ${n})
  endif()
endif()


function(p4est_example name files arg1 arg2)

add_executable(${name} ${files})
target_link_libraries(${name} PRIVATE P4EST::P4EST SC::SC)
if(MPI_FOUND)
  target_link_libraries(${name} PRIVATE MPI::MPI_C)
  add_test(NAME p4est:example:${name} COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${Ncpu} $<TARGET_FILE:${name}> ${arg1} ${arg2})
else()
  add_test(NAME p4est:example:${name} COMMAND $<TARGET_FILE:${name}> ${arg1} ${arg2})
endif()
set_tests_properties(p4est:example:${name} PROPERTIES TIMEOUT 10)

endfunction(p4est_example)


p4est_example(mesh2 mesh/mesh2.c unit 4)
# p4est_example(points2 points/points2.c unit 4)  # missing .pts input files
p4est_example(simple2 simple/simple2.c unit 4)

if(P4EST_HAVE_GETOPT_H)

foreach(n bricks timings)
  p4est_example(${n}2 timings/${n}2.c "" "")
endforeach()

p4est_example(loadconn2 timings/loadconn2.c "-l 4" "")

set(n spheres)
p4est_example(${n}2 "${n}/${n}2.c;${n}/p4est_${n}.c" "" "")

set(n particles)
p4est_example(${n}2 ${n}/${n}2.c "" "")

endif(P4EST_HAVE_GETOPT_H)

foreach(i 1 3 4)
  set(n p4est_step${i})
  p4est_example(${n} steps/${n}.c "" "")
endforeach()

p4est_example(p4est_step2 steps/p4est_step2.c ${CMAKE_CURRENT_SOURCE_DIR}/steps/hole_2d_cubit.inp 4)
