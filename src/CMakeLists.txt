include(CheckIncludeFile)
include(CheckSymbolExists)

if(Ep4est)
target_sources(p4est PRIVATE p4est_base.c p4est_connectivity.c p4est.c p4est_bits.c p4est_search.c p4est_algorithms.c p4est_communication.c p4est_ghost.c p4est_nodes.c p4est_vtk.c p4est_points.c p4est_geometry.c p4est_iterate.c p4est_lnodes.c p4est_mesh.c p4est_balance.c p4est_io.c p4est_connrefine.c p4est_wrap.c p4est_plex.c p4est_empty.c)
endif(Ep4est)

if(Ep8est)
target_sources(p8est PRIVATE p8est_connectivity.c p8est.c p8est_bits.c p8est_search.c p8est_algorithms.c p8est_communication.c p8est_ghost.c p8est_nodes.c p8est_vtk.c p8est_points.c p8est_geometry.c p8est_iterate.c p8est_lnodes.c p8est_mesh.c p8est_tets_hexes.c p8est_balance.c p8est_io.c p8est_connrefine.c p8est_wrap.c p8est_plex.c p8est_empty.c)
endif(Ep8est)

if(Ep6est AND Ep4est AND Ep8est)
target_sources(p6est PRIVATE p6est.c p6est_ghost.c p6est_lnodes.c p6est_profile.c p6est_vtk.c p6est_communication.c p6est_empty.c)
endif()
# --- generate p4est_config.h

set(CMAKE_REQUIRED_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_BINARY_DIR}/include)
if(NOT WIN32)
  set(CMAKE_REQUIRED_LIBRARIES m)
endif()
if(ZLIB_FOUND)
  set(CMAKE_REQUIRED_LIBRARIES ZLIB::ZLIB ${CMAKE_REQUIRED_LIBRARIES})
endif()

if(MPI_FOUND)
  set(CMAKE_REQUIRED_LIBRARIES MPI::MPI_C ${CMAKE_REQUIRED_LIBRARIES})
  set(P4EST_CC \"${MPI_C_COMPILER}\")
  set(P4EST_CPP ${MPI_C_COMPILER})
  set(P4EST_CXX \"${MPI_CXX_COMPILER}\")
  SET(P4EST_F77 \"${MPI_Fortran_COMPILER}\")
else()
  set(P4EST_CC \"${CMAKE_C_COMPILER}\")
  set(P4EST_CPP ${CMAKE_C_COMPILER})
  set(P4EST_CXX \"${CMAKE_CXX_COMPILER}\")
  SET(P4EST_F77 \"${CMAKE_Fortran_COMPILER}\")
endif()

string(APPEND P4EST_CPP " -E")
set(P4EST_CPP \"${P4EST_CPP}\")

set(P4EST_CFLAGS "${CMAKE_C_FLAGS} ${MPI_C_COMPILE_OPTIONS}")
set(P4EST_CFLAGS \"${P4EST_CFLAGS}\")

set(P4EST_CPPFLAGS \"\")

set(P4EST_FFLAGS "${CMAKE_Fortran_FLAGS} ${MPI_Fortran_COMPILE_OPTIONS}")
set(P4EST_FFLAGS \"${P4EST_FFLAGS}\")

set(P4EST_FLIBS \"${MPI_Fortran_LIBRARIES}\")

set(P4EST_LDFLAGS \"${MPI_C_LINK_FLAGS}\")
set(P4EST_LIBS \"${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${ZLIB_LIBRARIES} m\")

set(P4EST_ENABLE_BUILD_2D ${Ep4est})
set(P4EST_ENABLE_BUILD_3D ${Ep6est})
set(P4EST_ENABLE_BUILD_P6EST ${Ep8est})

set(P4EST_ENABLE_MEMALIGN 1)

if(MPI_FOUND)
  set(P4EST_ENABLE_MPI ${MPI_FOUND})
  check_symbol_exists(MPI_COMM_TYPE_SHARED mpi.h P4EST_ENABLE_MPICOMMSHARED)
  set(P4EST_ENABLE_MPIIO 1)
  check_symbol_exists(MPI_Init_thread mpi.h P4EST_ENABLE_MPITHREAD)
  check_symbol_exists(MPI_Win_allocate_shared mpi.h P4EST_ENABLE_MPIWINSHARED)
endif(MPI_FOUND)


check_symbol_exists(aligned_alloc stdlib.h P4EST_HAVE_ALIGNED_ALLOC)
check_include_file(arpa/inet.h P4EST_HAVE_ARPA_INET_H)
check_include_file(dlfcn.h P4EST_HAVE_DLFCN_H)
check_symbol_exists(fsync unistd.h P4EST_HAVE_FSYNC)
check_include_file(inttypes.h P4EST_HAVE_INTTYPES_H)

check_symbol_exists(pthread_create pthread.h HAVE_LPTHREAD)
check_symbol_exists(lua_createtable lua.h HAVE_LUA)

check_include_file(memory.h P4EST_HAVE_MEMORY_H)
check_include_file(netinet/in.h P4EST_HAVE_NETINET_IN_H)

check_symbol_exists(posix_memalign stdlib.h P4EST_HAVE_POSIX_MEMALIGN)
check_include_file(stdint.h P4EST_HAVE_STDINT_H)
check_include_file(stdlib.h P4EST_HAVE_STDLIB_H)
check_include_file(strings.h P4EST_HAVE_STRINGS_H)
check_include_file(string.h P4EST_HAVE_STRING_H)
check_include_file(sys/stat.h P4EST_HAVE_SYS_STAT_H)
check_include_file(sys/types.h P4EST_HAVE_SYS_TYPES_H)
check_include_file(unistd.h P4EST_HAVE_UNISTD_H)
check_symbol_exists(adler32_combine zlib.h P4EST_HAVE_ZLIB)

set(P4EST_WITH_BLAS ${BLAS_FOUND})
set(P4EST_WITH_LAPACK ${LAPACK_FOUND})

set(WITH_SC ${sc_ROOT})

configure_file(p4est_config.h.in ${PROJECT_BINARY_DIR}/include/p4est_config.h)
