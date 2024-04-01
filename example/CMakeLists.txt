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


function(p4est_example_test name dir)

p4est_example(${name}2 ${dir}/${name}2.c ${dir})

if(P4EST_ENABLE_P8EST)
  p8est_example(${name}3 ${dir}/${name}3.c ${dir})
endif()
endfunction(p4est_example_test)

function(p4est_copy_resource dir res_file)

configure_file(
  ${PROJECT_SOURCE_DIR}/${dir}/${res_file}
  ${PROJECT_BINARY_DIR}/example/${dir}/${res_file}
  COPYONLY
)
endfunction()

# --- setup examples

p4est_example_test(points points)
p4est_example_test(generate_points points)

if(P4EST_HAVE_GETOPT_H)

p4est_example_test(particles particles)
p4est_copy_resource(particles separt.pl)

set(n spheres)
p4est_example(${n}2 "${n}/${n}2.c;${n}/p4est_${n}.c" ${n})
if(P4EST_ENABLE_P8EST)
  p8est_example(${n}3 "${n}/${n}3.c;${n}/p8est_${n}.c" ${n})
endif()

p4est_example_test(bricks timings)
p4est_example_test(timings timings)
p4est_example_test(loadconn timings)

foreach(n IN ITEMS timana.awk timana.sh tsrana.awk tsrana.sh perfscript.sh)
  p4est_copy_resource(timings ${n})
endforeach()

p8est_example(tsearch3 timings/tsearch3.c "timings")

endif(P4EST_HAVE_GETOPT_H)

p4est_example_test(balance_seeds balance)
p4est_example_test(mesh mesh)
p4est_example_test(simple simple)

p4est_copy_resource(mesh conndebug.p8c)

p4est_example(periodicity3 mesh/periodicity3.c "mesh")

p4est_example_test(count_quadrants search)

foreach(i RANGE 1 5)
  p4est_example(p4est_step${i} steps/p4est_step${i}.c "steps")
  if(P4EST_ENABLE_P8EST)
    p8est_example(p8est_step${i} steps/p8est_step${i}.c "steps")
  endif()
endforeach()

foreach(n IN ITEMS cubit.inp cubit.jou gmsh.geo gmsh.inp)
  p4est_copy_resource(steps hole_2d_${n})
  if(P4EST_ENABLE_P8EST)
    p4est_copy_resource(steps hole_3d_${n})
  endif()
endforeach()

p4est_example_test(read_conn tetgen)
p4est_example_test(write_conn tetgen)
if(P4EST_ENABLE_P8EST)
  p8est_example(read_tetgen tetgen/read_tetgen.c "tetgen")
  p4est_copy_resource(tetgen p8est_box_tetgen.ele)
  p4est_copy_resource(tetgen p8est_box_tetgen.node)
endif()
