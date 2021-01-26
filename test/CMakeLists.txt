if(P4EST_HAVE_STDLIB_H)
  check_symbol_exists(random stdlib.h P4EST_HAVE_RANDOM)
  check_symbol_exists(srandom stdlib.h P4EST_HAVE_SRANDOM)
endif()

set(test test_comm test_hash)

if(Ep4est)
 list(APPEND tests  test_conn_transformation2 test_brick2 test_join2 test_conn_reduce2)
    if(P4EST_HAVE_ARPA_INET_H OR P4EST_HAVE_NETINET_IN_H)
      # missing htonl
      list(APPEND tests test_balance2 test_partition_corr2 test_coarsen2 test_balance_type2 test_lnodes2 test_plex2 test_connrefine2 test_subcomm2 test_replace2 test_search2 test_load2 test_order test_ghost2 test_iterate2 test_loadsave2 test_nodes2 test_partition2 test_quadrants2 test_complete_subtree test_valid2 test_conn_complete2 test_wrap2)
    endif()
    if(P4EST_HAVE_RANDOM AND P4EST_HAVE_SRANDOM)
      list(APPEND tests test_balance_seeds2)
    endif()
endif()



foreach(t ${tests})

add_executable(${t} ${t}.c)

target_link_libraries(${t} PRIVATE p4est SC::SC)
if(ZLIB_FOUND)
  target_link_libraries(${t} PRIVATE ZLIB::ZLIB)
endif()
if(MPI_FOUND)
  target_link_libraries(${t} PRIVATE MPI::MPI_C)
endif()
if(NOT WIN32)
  target_link_libraries(${t} PRIVATE m)
endif()

add_test(NAME p4est:${t} COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} $<TARGET_FILE:${t}>)

endforeach()
