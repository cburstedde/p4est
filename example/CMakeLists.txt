cmake_minimum_required(VERSION 3.22)
project(p4estExamples LANGUAGES C)

enable_testing()

include(CheckSymbolExists)
include(CheckIncludeFile)

if(PROJECT_IS_TOP_LEVEL)
  find_package(P4EST CONFIG REQUIRED)
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
# it is not intended to run examples as tests in the top-level project
if(NOT PROJECT_IS_TOP_LEVEL)
  return()
endif()


function(p4est_example name files dir)

add_executable(${name} ${files})
target_link_libraries(${name} PRIVATE P4EST::P4EST)
set_target_properties(${name}
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/example/${dir}"
  LABELS p4est
)

endfunction(p4est_example)


function(p8est_example name files dir)

add_executable(${name} ${files})
target_link_libraries(${name} PRIVATE P4EST::P4EST)
set_target_properties(${name}
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/example/${dir}"
  LABELS p8est
)

endfunction(p8est_example)


function(p4est_copy_resource dir res_file)
configure_file(
  ${PROJECT_SOURCE_DIR}/${dir}/${res_file}
  ${PROJECT_BINARY_DIR}/example/${dir}/${res_file}
  COPYONLY
)
endfunction()

# --- setup examples

p4est_example(points2 points/points2.c "point")
if(P4EST_ENABLE_P8EST)
  p8est_example(points3 points/points3.c "point")
endif()

p4est_example(generate_points2 points/generate_points2.c "point")
if(P4EST_ENABLE_P8EST)
  p8est_example(generate_points3 points/generate_points3.c "point")
endif()


if(P4EST_HAVE_GETOPT_H)

set(n particles)
p4est_example(${n}2 ${n}/${n}2.c ${n})
if(P4EST_ENABLE_P8EST)
  p8est_example(${n}3 ${n}/${n}3.c ${n})
endif()
p4est_copy_resource(${n} separt.pl)

set(n spheres)
p4est_example(${n}2 "${n}/${n}2.c;${n}/p4est_${n}.c" ${n})
if(P4EST_ENABLE_P8EST)
  p8est_example(${n}3 "${n}/${n}3.c;${n}/p8est_${n}.c" ${n})
endif()

foreach(n IN ITEMS bricks timings loadconn)
  p4est_example(${n}2 timings/${n}2.c "timings")
  if(P4EST_ENABLE_P8EST)
    p8est_example(${n}3 timings/${n}3.c ${n} "timings")
  endif()
endforeach()
foreach(n IN ITEMS timana tsrana)
  foreach(r IN ITEMS awk sh)
    p4est_copy_resource(timings ${n}.${r})
  endforeach()
endforeach()
p4est_copy_resource(timings perfscript.sh)

p8est_example(tsearch3 timings/tsearch3.c "timings")

endif(P4EST_HAVE_GETOPT_H)

set(n balance_seeds)
p4est_example(${n}2 balance/${n}2.c "balance")
if(P4EST_ENABLE_P8EST)
  p8est_example(${n}3 balance/${n}3.c "balance")
endif()

foreach(n IN ITEMS mesh simple)
  p4est_example(${n}2 ${n}/${n}2.c ${n})
  if(P4EST_ENABLE_P8EST)
    p8est_example(${n}3 ${n}/${n}3.c ${n})
  endif()
endforeach()
p4est_copy_resource(mesh conndebug.p8c)
p4est_example(periodicity3 mesh/periodicity3.c "mesh")

p4est_example(count_quadrants2 search/count_quadrants2.c "search")
if(P4EST_ENABLE_P8EST)
  p8est_example(count_quadrants3 search/count_quadrants3.c "search")
endif()

set(s steps)
foreach(i RANGE 1 5)
  set(n p4est_step${i})
  p4est_example(${n} ${s}/${n}.c ${s})

  if(P4EST_ENABLE_P8EST)
    set(n p8est_step${i})
    p8est_example(${n} ${s}/${n}.c ${s})
  endif()
endforeach()

foreach(n IN ITEMS cubit.inp cubit.jou gmsh.geo gmsh.inp)
  p4est_copy_resource(${s} hole_2d_${n})
endforeach()
if(P4EST_ENABLE_P8EST)
  foreach(n IN ITEMS cubit.inp cubit.jou gmsh.geo gmsh.inp)
    p4est_copy_resource(${s} hole_3d_${n})
  endforeach()
  endif()

set(t tetgen)
foreach(n IN ITEMS read_conn write_conn)
  p4est_example(${n}2 ${t}/${n}2.c ${t})
  if(P4EST_ENABLE_P8EST)
    p8est_example(${n}3 ${t}/${n}3.c ${t})
  endif()
endforeach()
if(P4EST_ENABLE_P8EST)
  p8est_example(read_${t} ${t}/read_${t}.c ${t})
  foreach(r IN ITEMS ele node)
    p4est_copy_resource(${t} p8est_box_${t}.${r})
  endforeach()
endif()
