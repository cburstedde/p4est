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
function(p4est_example name files dir)

add_executable(${name} ${files})
target_link_libraries(${name} PRIVATE P4EST::P4EST)
set_target_properties(${name}
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/${dir}"
  LABELS p4est
)
endfunction(p4est_example)


function(p8est_example name files dir)

add_executable(${name} ${files})
target_link_libraries(${name} PRIVATE P4EST::P4EST)
set_target_properties(${name}
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/${dir}"
  LABELS p8est
)
endfunction(p8est_example)


function(p4est_p8est_example name dir)

p4est_example(${name}2 ${dir}/${name}2.c ${dir})

if(enable_p8est)
  p8est_example(${name}3 ${dir}/${name}3.c ${dir})
endif()
endfunction(p4est_p8est_example)

function(p4est_copy_resource dir res_file)

configure_file(
  ${PROJECT_SOURCE_DIR}/${dir}/${res_file}
  ${PROJECT_BINARY_DIR}/${dir}/${res_file}
  COPYONLY
)
endfunction()

# --- setup examples

p4est_p8est_example(points points)
p4est_p8est_example(generate_points points)

if(P4EST_HAVE_GETOPT_H)

p4est_p8est_example(particles particles)
p4est_copy_resource(particles separt.pl)

p4est_p8est_example(bricks timings)
p4est_p8est_example(timings timings)
p4est_p8est_example(loadconn timings)
foreach(n IN ITEMS timana.awk timana.sh tsrana.awk tsrana.sh perfscript.sh)
  p4est_copy_resource(timings ${n})
endforeach()

if(enable_p8est)
  p8est_example(tsearch3 timings/tsearch3.c "timings")
endif()

endif(P4EST_HAVE_GETOPT_H)

p4est_p8est_example(balance_seeds balance)
p4est_p8est_example(balance_corner balance)
p4est_p8est_example(mesh mesh)
p4est_p8est_example(simple simple)
p4est_copy_resource(mesh conndebug.p8c)

if(enable_p8est)
  p4est_example(periodicity3 mesh/periodicity3.c "mesh")
endif()

p4est_p8est_example(count_quadrants search)

p4est_p8est_example(read_conn tetgen)
p4est_p8est_example(write_conn tetgen)
if(enable_p8est)
  p8est_example(read_tetgen tetgen/read_tetgen.c "tetgen")
  p4est_copy_resource(tetgen p8est_box_tetgen.ele)
  p4est_copy_resource(tetgen p8est_box_tetgen.node)
endif()

p4est_example(spheres2 "spheres/spheres2.c;spheres/p4est_spheres.c" "spheres")
if(enable_p8est)
  p8est_example(spheres3 "spheres/spheres3.c;spheres/p8est_spheres.c" "spheres")
endif()

foreach(i RANGE 1 5)
  p4est_example(p4est_step${i} steps/p4est_step${i}.c "steps")
  if(enable_p8est)
    p8est_example(p8est_step${i} steps/p8est_step${i}.c "steps")
  endif()
endforeach()

foreach(n IN ITEMS cubit.inp cubit.jou gmsh.geo gmsh.inp)
  p4est_copy_resource(steps hole_2d_${n})
  if(enable_p8est)
    p4est_copy_resource(steps hole_3d_${n})
  endif()
endforeach()
