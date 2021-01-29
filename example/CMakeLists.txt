
if(NOT (P4EST_HAVE_ARPA_INET_H OR P4EST_HAVE_NETINET_IN_H OR P4EST_HAVE_WINSOCK2_H))
  return()
endif()


function(p4est_example name files)

add_executable(${name} ${files})
target_link_libraries(${name} PRIVATE P4EST::P4EST SC::SC)
if(MPI_FOUND)
  target_link_libraries(${name} PRIVATE MPI::MPI_C)
  add_test(NAME p4est:example:${name} COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${Ncpu} $<TARGET_FILE:${name}>)
else()
  add_test(NAME p4est:example:${name} COMMAND $<TARGET_FILE:${name}>)
endif()
set_tests_properties(p4est:example:${name} PROPERTIES TIMEOUT 10)

endfunction(p4est_example)


foreach(n mesh points simple)
  p4est_example(${n}2 ${n}/${n}2.c)
endforeach()

if(P4EST_HAVE_GETOPT_H)

foreach(n bricks loadconn timings)
  p4est_example(${n}2 timings/${n}2.c)
endforeach()

set(n spheres)
p4est_example(${n}2 "${n}/${n}2.c;${n}/p4est_${n}.c")

set(n particles)
p4est_example(${n}2 ${n}/${n}2.c)

endif()

foreach(i 1 2 3 4)
  set(n p4est_step${i})
  p4est_example(${n} steps/${n}.c)
endforeach()
