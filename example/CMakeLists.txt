option(mpi "use MPI" off)

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
# for everything that uses sc_options.h
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

function(p4est_example name files dir arg1 arg2)

add_executable(${name} ${files})
target_link_libraries(${name} PRIVATE P4EST::P4EST SC::SC)
set_target_properties(${name}
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/example/${dir}"
  LABELS p4est
)
endfunction(p4est_example)


function(p8est_example name files dir arg1 arg2)

add_executable(${name} ${files})
target_link_libraries(${name} PRIVATE P4EST::P4EST SC::SC)
set_target_properties(${name}
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/example/${dir}"
  LABELS p8est
)

endfunction(p8est_example)


# --- setup examples
if(P4EST_HAVE_GETOPT_H)

set(n particles)
p4est_example(${n}2 ${n}/${n}2.c ${n} "" "")
if(enable_p8est)
  p8est_example(${n}3 ${n}/${n}3.c ${n} "" "")
endif()

set(n p4est_step3)
p4est_example(${n} steps/${n}.c "steps" "" "")
if(enable_p8est)
  set(n p8est_step3)
  p8est_example(${n} steps/${n}.c "steps" "" "")
endif()

set(n spheres)
p4est_example(${n}2 "${n}/${n}2.c;${n}/p4est_${n}.c" ${n} "" "")
if(enable_p8est)
  p8est_example(${n}3 "${n}/${n}3.c;${n}/p8est_${n}.c" ${n} "" "")
endif()

foreach(n bricks timings loadconn)
  p4est_example(${n}2 timings/${n}2.c "timings" "" "")
  if(enable_p8est)
    p8est_example(${n}3 timings/${n}3.c ${n} "timings" "" "")
  endif()
endforeach()

p8est_example(tsearch3 timings/tsearch3.c "timings" "" "")

endif(P4EST_HAVE_GETOPT_H)

set(n balance_seeds)
p4est_example(${n}2 balance/${n}2.c "balance" "" "")
if(enable_p8est)
  p8est_example(${n}3 balance/${n}3.c "balance" "" "")
endif()

foreach(n mesh points simple)
  p4est_example(${n}2 ${n}/${n}2.c ${n} unit 4)
  if(enable_p8est)
    p8est_example(${n}3 ${n}/${n}3.c ${n} unit 4)
  endif()
endforeach()
p4est_example(periodicity3 mesh/periodicity3.c "mesh" "" "")

p4est_example(count_quadrants2 search/count_quadrants2.c "search" "" "")
if(enable_p8est)
  p8est_example(count_quadrants3 search/count_quadrants3.c "search" "" "")
endif()

foreach(i 1 2 4)
  set(n p4est_step${i})
  p4est_example(${n} steps/${n}.c  "steps" "" "")

  if(enable_p8est)
    set(n p8est_step${i})
    p8est_example(${n} steps/${n}.c "steps" "" "")
  endif()
endforeach()

foreach(n read_conn write_conn)
  p4est_example(${n}2 tetgen/${n}2.c "tetgen" "" "")
  if(enable_p8est)
    p8est_example(${n}3 tetgen/${n}3.c "tetgen" "" "")
  endif()
endforeach()
if(enable_p8est)
  p8est_example(read_tetgen tetgen/read_tetgen.c "tetgen" "" "")
endif()

get_property(target_names DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY NAME)

if(WIN32 AND CMAKE_VERSION VERSION_GREATER_EQUAL 3.22)
  set_property(TARGET ${target_names} PROPERTY
  ENVIRONMENT_MODIFICATION "PATH=path_list_append:${P4EST_INCLUDE_DIRS}/../bin;PATH=path_list_append:${SC_INCLUDE_DIRS}/../bin"
  )
endif()
