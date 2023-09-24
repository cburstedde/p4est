configure_file(${PROJECT_SOURCE_DIR}/cmake/CTestCustom.cmake ${CMAKE_BINARY_DIR})
option(mpi "use MPI" off)

# --- find external libraries
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../cmake)

# --- get system capabilities

check_symbol_exists(sqrt math.h P4EST_NONEED_M)
if(NOT P4EST_NONEED_M)
  set(CMAKE_REQUIRED_LIBRARIES m)
  check_symbol_exists(sqrt math.h P4EST_NEED_M)
endif()

# --- get system capabilities
check_include_file(arpa/inet.h P4EST_HAVE_ARPA_INET_H)
check_include_file(netinet/in.h P4EST_HAVE_NETINET_IN_H)
if(WIN32 AND NOT P4EST_HAVE_ARPA_INET_H AND NOT P4EST_HAVE_NETINET_IN_H)
  check_include_file(Winsock2.h P4EST_HAVE_WINSOCK2_H)
  target_link_libraries(SC::SC INTERFACE wsock32 ws2_32) # Iphlpapi
endif()

if(NOT (P4EST_HAVE_ARPA_INET_H OR P4EST_HAVE_NETINET_IN_H OR P4EST_HAVE_WINSOCK2_H))
  message(FATAL_ERROR "A networking library was not found.")
endif()

check_include_file(unistd.h P4EST_HAVE_UNISTD_H)
if(P4EST_HAVE_UNISTD_H)
  check_include_file(getopt.h P4EST_HAVE_GETOPT_H)
endif()

cmake_host_system_information(RESULT Ncpu QUERY NUMBER_OF_PHYSICAL_CORES)
if(Ncpu LESS 2)
  include(ProcessorCount)
  ProcessorCount(n)
  if(n GREATER Ncpu)
    set(Ncpu ${n})
  endif()
endif()


# --- helper functions
# it is not intended to run examples as tests

function(p4est_example name files arg1 arg2)

add_executable(${name} ${files})
target_link_libraries(${name} PRIVATE P4EST::P4EST SC::SC)

if(mpi)
  add_test(NAME ${name} COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${Ncpu} $<TARGET_FILE:${name}> ${arg1} ${arg2})
else()
  add_test(NAME ${name} COMMAND ${name} ${arg1} ${arg2})
endif()
set_property(TEST ${name} PROPERTY LABELS p8est)

endfunction(p4est_example)


function(p8est_example name files arg1 arg2)

add_executable(${name} ${files})
target_link_libraries(${name} PRIVATE P4EST::P4EST SC::SC)

if(mpi)
  add_test(NAME ${name} COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${Ncpu} $<TARGET_FILE:${name}> ${arg1} ${arg2})
else()
  add_test(NAME ${name} COMMAND ${name} ${arg1} ${arg2})
endif()
set_property(TEST ${name} PROPERTY LABELS p8est)

endfunction(p8est_example)


# --- setup examples

set(n balance_seeds)
p4est_example(${n}2 balance/${n}2.c "" "")
if(enable_p8est)
  p8est_example(${n}3 balance/${n}3.c "" "")
endif()

set(n mesh)
p4est_example(${n}2 ${n}/${n}2.c unit 4)
if(enable_p8est)
  p8est_example(${n}3 ${n}/${n}3.c unit 4)
endif()

p4est_example(periodicity3 mesh/periodicity3.c "" "")
# p4est_example(points2 points/points2.c unit 4)  # missing .pts input files

p4est_example(simple2 simple/simple2.c unit 4)
if(enable_p8est)
  p8est_example(simple3 simple/simple3.c unit 4)
endif()


if(P4EST_HAVE_GETOPT_H)

foreach(n bricks timings)
  p4est_example(${n}2 timings/${n}2.c "" "")
  if(enable_p8est)
    p8est_example(${n}3 timings/${n}3.c "" "")
  endif()
endforeach()

# optparse not recognizing -l --level
# p4est_example(loadconn2 timings/loadconn2.c "-l 4" "")

set(n spheres)
p4est_example(${n}2 "${n}/${n}2.c;${n}/p4est_${n}.c" "" "")
if(enable_p8est)
  p8est_example(${n}3 "${n}/${n}3.c;${n}/p8est_${n}.c" "" "")
endif()

set(n particles)
p4est_example(${n}2 ${n}/${n}2.c "" "")
if(enable_p8est)
  p8est_example(${n}3 ${n}/${n}3.c "" "")
endif()

endif(P4EST_HAVE_GETOPT_H)

foreach(i 1 3 4)
  set(n p4est_step${i})
  p4est_example(${n} steps/${n}.c "" "")

  if(enable_p8est)
    set(n p8est_step${i})
    p8est_example(${n} steps/${n}.c "" "")
  endif()
endforeach()

# optparse not recognizing -l --level
# p4est_example(p4est_step2 steps/p4est_step2.c ${CMAKE_CURRENT_SOURCE_DIR}/steps/hole_2d_cubit.inp 4)
# set_property(TEST p4est:p4est_step2 PROPERTY TIMEOUT 300)

get_property(test_names DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY TESTS)

if(mpi)
  set_property(TEST ${test_names} PROPERTY RESOURCE_LOCK cpu_mpi)
endif()

set_property(TEST ${test_names} PROPERTY TIMEOUT 60)

if(WIN32 AND CMAKE_VERSION VERSION_GREATER_EQUAL 3.22)
  set_property(TEST ${test_names} PROPERTY
  ENVIRONMENT_MODIFICATION "PATH=path_list_append:${P4EST_INCLUDE_DIRS}/../bin;PATH=path_list_append:${SC_INCLUDE_DIRS}/../bin"
  )
endif()
