include(CheckIncludeFile)
include(CheckSymbolExists)
include(ProcessorCount)

# on some platforms e.g. ARM, we have to try a few ways to get CPU count > 1 for multi-core systems
cmake_host_system_information(RESULT Ncpu QUERY NUMBER_OF_PHYSICAL_CORES)
if(Ncpu LESS 2)
  ProcessorCount(n)
  if(n GREATER Ncpu)
    set(Ncpu ${n})
  endif()
endif()
message(STATUS "${PROJECT_NAME} CMake ${CMAKE_VERSION} ${Ncpu} threads")

# --- generate p4est_config.h

set(CMAKE_REQUIRED_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_BINARY_DIR}/include)
set(CMAKE_REQUIRED_LIBRARIES)

if(MPI_FOUND)
  set(CMAKE_REQUIRED_LIBRARIES MPI::MPI_C)
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

set(P4EST_ENABLE_BUILD_2D true CACHE BOOL "p4est is always used")
set(P4EST_ENABLE_BUILD_3D ${Ep6est})
set(P4EST_ENABLE_BUILD_P6EST ${Ep8est})

set(P4EST_ENABLE_MEMALIGN 1)

foreach(t FSYNC INTTYPES_H POSIX_MEMALIGN STDINT_H STDLIB_H STRINGS_H STRING_H SYS_STAT_H SYS_TYPES_H UNISTD_H ZLIB)
  if(DEFINED SC_HAVE_${t})
    set(P4EST_HAVE_${t} ${SC_HAVE_${t}} CACHE BOOL SC_HAVE_${t})
  endif()
endforeach()

foreach(t MPICOMMSHARED MPITHREAD MPIWINSHARED)
  if(DEFINED SC_ENABLE_${t})
    set(P4EST_ENABLE_${t} ${SC_ENABLE_${t}} CACHE BOOL SC_ENABLE_${t})
  endif()
endforeach()

foreach(t NONEED_M NEED_M)
  if(DEFINED SC_${t})
    set(P4EST_${t} ${SC_${t}} CACHE BOOL SC_${t})
  endif()
endforeach()

if(MPI_FOUND)
  set(P4EST_ENABLE_MPI ${MPI_FOUND})
  check_symbol_exists(MPI_COMM_TYPE_SHARED mpi.h P4EST_ENABLE_MPICOMMSHARED)
  set(P4EST_ENABLE_MPIIO 1)
  check_symbol_exists(MPI_Init_thread mpi.h P4EST_ENABLE_MPITHREAD)
  check_symbol_exists(MPI_Win_allocate_shared mpi.h P4EST_ENABLE_MPIWINSHARED)
endif(MPI_FOUND)

check_symbol_exists(sqrt math.h P4EST_NONEED_M)
if(NOT P4EST_NONEED_M)
  set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} m)
  check_symbol_exists(sqrt math.h P4EST_NEED_M)
endif()

check_include_file(arpa/inet.h P4EST_HAVE_ARPA_INET_H)
check_include_file(netinet/in.h P4EST_HAVE_NETINET_IN_H)
if(WIN32 AND NOT P4EST_HAVE_ARPA_INET_H AND NOT P4EST_HAVE_NETINET_IN_H)
  check_include_file(Winsock2.h P4EST_HAVE_WINSOCK2_H)
  set(WINSOCK_LIBRARIES wsock32 ws2_32) # Iphlpapi
endif()

check_include_file(dlfcn.h P4EST_HAVE_DLFCN_H)

check_symbol_exists(fsync unistd.h P4EST_HAVE_FSYNC)
check_include_file(inttypes.h P4EST_HAVE_INTTYPES_H)

check_symbol_exists(pthread_create pthread.h HAVE_LPTHREAD)
check_symbol_exists(lua_createtable lua.h HAVE_LUA)

check_include_file(memory.h P4EST_HAVE_MEMORY_H)

check_symbol_exists(posix_memalign stdlib.h P4EST_HAVE_POSIX_MEMALIGN)
check_include_file(stdint.h P4EST_HAVE_STDINT_H)
check_include_file(stdlib.h P4EST_HAVE_STDLIB_H)
check_include_file(strings.h P4EST_HAVE_STRINGS_H)
check_include_file(string.h P4EST_HAVE_STRING_H)
check_include_file(sys/stat.h P4EST_HAVE_SYS_STAT_H)
check_include_file(sys/types.h P4EST_HAVE_SYS_TYPES_H)
check_include_file(unistd.h P4EST_HAVE_UNISTD_H)

if(ZLIB_FOUND)
  set(CMAKE_REQUIRED_LIBRARIES ZLIB::ZLIB)
  check_symbol_exists(adler32_combine zlib.h P4EST_HAVE_ZLIB)
endif()

configure_file(src/p4est_config.h.in ${PROJECT_BINARY_DIR}/include/p4est_config.h)
