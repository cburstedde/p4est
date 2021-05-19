if(P4EST_HAVE_STDLIB_H)
  check_symbol_exists(random stdlib.h P4EST_HAVE_RANDOM)
  check_symbol_exists(srandom stdlib.h P4EST_HAVE_SRANDOM)
endif()

set(test test_comm test_hash test_order test_complete_subtree)

list(APPEND tests test_conn_transformation2 test_brick2 test_join2 test_conn_reduce2)
if(P4EST_HAVE_ARPA_INET_H OR P4EST_HAVE_NETINET_IN_H OR P4EST_HAVE_WINSOCK2_H)
  # htonl
  list(APPEND tests test_balance2 test_partition_corr2 test_coarsen2 test_balance_type2 test_lnodes2 test_plex2 test_connrefine2 test_search2 test_subcomm2 test_replace2 test_ghost2 test_iterate2 test_nodes2 test_partition2 test_quadrants2 test_valid2 test_conn_complete2 test_wrap2)

  if(P4EST_HAVE_GETOPT_H)
    list(APPEND tests test_load2 test_loadsave2)
  endif()

endif()
if(P4EST_HAVE_RANDOM AND P4EST_HAVE_SRANDOM)
  list(APPEND tests test_balance_seeds2)
endif()


if(enable_p8est)
  list(APPEND tests test_conn_transformation3 test_brick3 test_join3 test_conn_reduce3)
  if(P4EST_HAVE_ARPA_INET_H OR P4EST_HAVE_NETINET_IN_H OR P4EST_HAVE_WINSOCK2_H)
    # htonl
    list(APPEND tests test_balance3 test_partition_corr3 test_coarsen3 test_balance_type3 test_lnodes3 test_plex3 test_connrefine3 test_subcomm3 test_replace3 test_ghost3 test_iterate3 test_nodes3 test_partition3 test_quadrants3 test_valid3 test_conn_complete3 test_wrap3)
  endif()

  if(P4EST_HAVE_GETOPT_H)
    list(APPEND tests test_load3)
  endif()
endif()


if(enable_p6est AND enable_p8est AND P4EST_HAVE_GETOPT_H)
if(P4EST_HAVE_ARPA_INET_H OR P4EST_HAVE_NETINET_IN_H OR P4EST_HAVE_WINSOCK2_H)
  # htonl
  add_executable(p6est_test_all test_all6.c)
  target_link_libraries(p6est_test_all PRIVATE p4est)

  if(P4EST_HAVE_WINSOCK2_H)
    target_link_libraries(p6est_test_all PRIVATE ${WINSOCK_LIBRARIES})
  endif()

  if(MPI_FOUND)
    add_test(NAME p6est:test_all COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${Ncpu} $<TARGET_FILE:p6est_test_all>)
  else()
    add_test(NAME p6est:test_all COMMAND $<TARGET_FILE:p6est_test_all>)
  endif()
  set_tests_properties(p6est:test_all PROPERTIES
    TIMEOUT 60
    RESOURCE_LOCK cpu_mpi)
endif()
endif()


foreach(t ${tests})

add_executable(${t} ${t}.c)
target_link_libraries(${t} PRIVATE p4est)

if(MPI_FOUND)
  add_test(NAME p4est:${t} COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${Ncpu} $<TARGET_FILE:${t}>)
else()
  add_test(NAME p4est:${t} COMMAND $<TARGET_FILE:${t}>)
endif()
set_tests_properties(p4est:${t} PROPERTIES
  TIMEOUT 60
  RESOURCE_LOCK cpu_mpi)

endforeach()

get_property(_tests DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY TESTS)

if(p4est:test_loadsave2 IN_LIST _tests)
  set_tests_properties(p4est:test_loadsave2 PROPERTIES TIMEOUT 600)
endif()
