# --- get system capabilities
if(WIN32 AND NOT P4EST_HAVE_ARPA_INET_H AND NOT P4EST_HAVE_NETINET_IN_H)
  target_link_libraries(SC::SC INTERFACE wsock32 ws2_32) # Iphlpapi
endif()

if(NOT (P4EST_HAVE_ARPA_INET_H OR P4EST_HAVE_NETINET_IN_H OR P4EST_HAVE_WINSOCK2_H))
  message(FATAL_ERROR "A networking library was not found.")
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

function(p4est_example name files dir)

add_executable(${name} ${files})
target_link_libraries(${name} PRIVATE P4EST::P4EST SC::SC)
set_target_properties(${name}
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/example/${dir}"
  LABELS p4est
)
endfunction(p4est_example)

function(p8est_example name files dir)

add_executable(${name} ${files})
target_link_libraries(${name} PRIVATE P4EST::P4EST SC::SC)
set_target_properties(${name}
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/example/${dir}"
  LABELS p8est
)
endfunction(p8est_example)

function(p4est_copy_resource dir res_file)
configure_file(
  ${PROJECT_SOURCE_DIR}/example/${dir}/${res_file}
  ${CMAKE_BINARY_DIR}/example/${dir}/${res_file}
  COPYONLY
)
endfunction()

# --- setup examples

p4est_example(mesh2 mesh/mesh2.c unit 4)
if(P4EST_ENABLE_P8EST)
  p4est_example(mesh3 mesh/mesh3.c unit 4)
  p4est_example(periodicity3 mesh/periodicity3.c unit 4)
  configure_file(
    ${CMAKE_SOURCE_DIR}/mesh/conndebug.p8c
    ${CMAKE_BINARY_DIR}/conndebug.p8c
    COPYONLY)
endif()

# p4est_example(points2 points/points2.c unit 4)  # missing .pts input files

p4est_example(simple2 simple/simple2.c unit 4)
if(P4EST_ENABLE_P8EST)
  p8est_example(simple3 simple/simple3.c unit 4)
endif()


if(P4EST_HAVE_GETOPT_H)

set(n particles)
p4est_example(${n}2 ${n}/${n}2.c ${n})
if(enable_p8est)
  p8est_example(${n}3 ${n}/${n}3.c ${n})
endif()
p4est_copy_resource(${n} separt.pl)

set(n p4est_step3)
p4est_example(${n} steps/${n}.c "steps")
if(enable_p8est)
  set(n p8est_step3)
  p8est_example(${n} steps/${n}.c "steps")
endif()

set(n spheres)
p4est_example(${n}2 "${n}/${n}2.c;${n}/p4est_${n}.c" ${n})
if(enable_p8est)
  p8est_example(${n}3 "${n}/${n}3.c;${n}/p8est_${n}.c" ${n})
endif()

foreach(n bricks timings loadconn)
  p4est_example(${n}2 timings/${n}2.c "timings")
  if(enable_p8est)
    p8est_example(${n}3 timings/${n}3.c ${n} "timings")
  endif()
endforeach()
foreach(n timana tsrana)
  foreach(r awk sh)
    p4est_copy_resource(timings ${n}.${r})
  endforeach()
endforeach()
p4est_copy_resource(timings perfscript.sh)

p8est_example(tsearch3 timings/tsearch3.c "timings")

endif(P4EST_HAVE_GETOPT_H)

set(n balance_seeds)
p4est_example(${n}2 balance/${n}2.c "balance")
if(enable_p8est)
  p8est_example(${n}3 balance/${n}3.c "balance")
endif()

foreach(n mesh points simple)
  p4est_example(${n}2 ${n}/${n}2.c ${n})
  if(enable_p8est)
    p8est_example(${n}3 ${n}/${n}3.c ${n})
  endif()
endforeach()
p4est_copy_resource(mesh conndebug.p8c)
p4est_example(periodicity3 mesh/periodicity3.c "mesh")

p4est_example(count_quadrants2 search/count_quadrants2.c "search")
if(enable_p8est)
  p8est_example(count_quadrants3 search/count_quadrants3.c "search")
endif()

set(s steps)
foreach(i 1 2 4)
  set(n p4est_step${i})
  p4est_example(${n} ${s}/${n}.c ${s})

  if(enable_p8est)
    set(n p8est_step${i})
    p8est_example(${n} ${s}/${n}.c ${s})
  endif()
endforeach()

foreach(n cubit.inp cubit.jou gmsh.geo gmsh.inp)
  p4est_copy_resource(${s} hole_2d_${n})
endforeach()
if(enable_p8est)
  foreach(n cubit.inp cubit.jou gmsh.geo gmsh.inp)
    p4est_copy_resource(${s} hole_3d_${n})
  endforeach()
  endif()

set(t tetgen)
foreach(n read_conn write_conn)
  p4est_example(${n}2 ${t}/${n}2.c ${t})
  if(enable_p8est)
    p8est_example(${n}3 ${t}/${n}3.c ${t})
  endif()
endforeach()
if(enable_p8est)
  p8est_example(read_${t} ${t}/read_${t}.c ${t})
  foreach(r ele node)
    p4est_copy_resource(${t} p8est_box_${t}.${r})
  endforeach()
endif()

get_property(target_names DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY NAME)

if(WIN32 AND CMAKE_VERSION VERSION_GREATER_EQUAL 3.22)
  set_property(TARGET ${target_names} PROPERTY
  ENVIRONMENT_MODIFICATION "PATH=path_list_append:${P4EST_INCLUDE_DIRS}/../bin;PATH=path_list_append:${SC_INCLUDE_DIRS}/../bin"
  )
endif()
