include(CTest)

if(P4EST_HAVE_STDLIB_H)
  check_symbol_exists(random stdlib.h P4EST_HAVE_RANDOM)
  check_symbol_exists(srandom stdlib.h P4EST_HAVE_SRANDOM)
endif()

set(p4est_tests test_comm test_hash test_order test_complete_subtree
test_conn_transformation2 test_brick2 test_join2 test_conn_reduce2
)

list(APPEND tests test_conn_transformation2 test_brick2 test_join2 test_conn_reduce2 test_version)
if(P4EST_HAVE_ARPA_INET_H OR P4EST_HAVE_NETINET_IN_H OR P4EST_HAVE_WINSOCK2_H)
  # htonl
  list(APPEND p4est_tests test_balance2 test_partition_corr2 test_coarsen2 test_balance_type2 test_lnodes2 test_plex2 test_connrefine2 test_search2 test_subcomm2 test_replace2 test_ghost2 test_iterate2 test_nodes2 test_partition2 test_quadrants2 test_valid2 test_conn_complete2 test_wrap2)

  if(P4EST_HAVE_GETOPT_H)
    list(APPEND p4est_tests test_load2 test_loadsave2)
  endif()
endif()
if(P4EST_HAVE_RANDOM AND P4EST_HAVE_SRANDOM)
  list(APPEND p4est_tests test_balance_seeds2)
endif()


if(enable_p8est)
  set(p8est_tests test_conn_transformation3 test_brick3 test_join3 test_conn_reduce3 test_mesh_corners3)
  if(P4EST_HAVE_ARPA_INET_H OR P4EST_HAVE_NETINET_IN_H OR P4EST_HAVE_WINSOCK2_H)
    # htonl
    list(APPEND p8est_tests test_balance3 test_partition_corr3 test_coarsen3 test_balance_type3 test_lnodes3 test_plex3 test_connrefine3 test_subcomm3 test_replace3 test_ghost3 test_iterate3 test_nodes3 test_partition3 test_quadrants3 test_valid3 test_conn_complete3 test_wrap3)
  endif()

  if(P4EST_HAVE_GETOPT_H)
    list(APPEND p8est_tests test_load3)
  endif()
endif()


if(enable_p6est AND enable_p8est AND P4EST_HAVE_GETOPT_H)
  if(P4EST_HAVE_ARPA_INET_H OR P4EST_HAVE_NETINET_IN_H OR P4EST_HAVE_WINSOCK2_H)
    set(p6est_tests test_all6)
  endif()
endif()


foreach(t IN LISTS p4est_tests p6est_tests p8est_tests)

add_executable(${t} ${t}.c)

target_link_libraries(${t} PRIVATE
p4est
$<$<BOOL:${P4EST_HAVE_WINSOCK2_H}>:${WINSOCK_LIBRARIES}>
)

if(P4EST_ENABLE_MPI)
  add_test(NAME ${t} COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_NUMPROC_MAX} $<TARGET_FILE:${t}>)
else()
  add_test(NAME ${t} COMMAND ${t})
endif()

endforeach()

set_property(TEST ${p4est_tests} ${p6est_tests} ${p8est_tests} PROPERTY TIMEOUT 60)

if(P4EST_ENABLE_MPI)
  set_property(TEST ${p4est_tests} ${p6est_tests} ${p8est_tests} PROPERTY RESOURCE_LOCK cpu_mpi)
endif()

set_property(TEST ${p4est_tests} PROPERTY LABELS p4est)
set_property(TEST ${p6est_tests} PROPERTY LABELS p6est)
set_property(TEST ${p8est_tests} PROPERTY LABELS p8est)

message(DEBUG "p4est tests: ${p4est_tests}
p6est tests: ${p6est_tests}
p8est tests: ${p8est_tests}
")


if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.22)
  if(WIN32)
    set_property(TEST ${p4est_tests} ${p6est_tests} ${p8est_tests} PROPERTY
    ENVIRONMENT_MODIFICATION "PATH=path_list_append:$<TARGET_FILE_DIR:p4est>"
    )
  elseif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    set_property(TEST ${p4est_tests} ${p6est_tests} ${p8est_tests} PROPERTY
      ENVIRONMENT_MODIFICATION "LD_LIBRARY_PATH=path_list_append:$<TARGET_FILE_DIR:p4est>"
    )
  elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    set_property(TEST ${p4est_tests} ${p6est_tests} ${p8est_tests} PROPERTY
      ENVIRONMENT_MODIFICATION "DYLD_LIBRARY_PATH=path_list_append:$<TARGET_FILE_DIR:p4est>"
    )
  endif()
endif()

if(test_loadsave2 IN_LIST p4est_tests)
  set_property(TEST test_loadsave2 PROPERTY TIMEOUT 600)
endif()
