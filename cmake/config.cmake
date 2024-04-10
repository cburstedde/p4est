include(CheckIncludeFile)
include(CheckSymbolExists)
include(ProcessorCount)

# --- retrieve library interface version from configuration file
file(STRINGS config/p4est_soversion.in P4EST_SOVERSION_READ
             REGEX "^[ \t]*P4EST_SOVERSION *= *[0-9:]+")
string(REGEX REPLACE ".*= *([0-9]+):([0-9]+):([0-9]+)" "\\1.\\2.\\3"
             P4EST_SOVERSION ${P4EST_SOVERSION_READ})
message(STATUS "p4est SOVERSION configured as ${P4EST_SOVERSION}")

# on some platforms e.g. ARM, we have to try a few ways to get CPU count > 1 for multi-core systems
cmake_host_system_information(RESULT Ncpu QUERY NUMBER_OF_PHYSICAL_CORES)
if(Ncpu LESS 2)
  ProcessorCount(n)
  if(n GREATER Ncpu)
    set(Ncpu ${n})
  endif()
  set(MPIEXEC_NUMPROC_MAX 1)
else()
  set(MPIEXEC_NUMPROC_MAX 2)
endif()

# --- set global compile environment

# Build all targets with -fPIC so that libsc itself can be linked as a
# shared library, or linked into a shared library.
include(CheckPIESupported)
check_pie_supported()
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# --- generate p4est_config.h

set(CMAKE_REQUIRED_INCLUDES)
set(CMAKE_REQUIRED_LIBRARIES)

set(P4EST_ENABLE_MPI ${SC_ENABLE_MPI} CACHE BOOL "P4EST MPI support" FORCE)
set(P4EST_ENABLE_MPIIO ${SC_ENABLE_MPIIO} CACHE BOOL "P4EST MPI-IO support" FORCE)
set(P4EST_HAVE_ZLIB ${SC_HAVE_ZLIB} CACHE BOOL "P4EST zlib support" FORCE)
set(P4EST_HAVE_JSON ${SC_HAVE_JSON} CACHE BOOL "P4EST json support" FORCE)

if(P4EST_ENABLE_MPI)
  find_package(MPI REQUIRED COMPONENTS C)
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

set(P4EST_CFLAGS "${CMAKE_C_FLAGS}\ ${MPI_C_COMPILE_OPTIONS}")
set(P4EST_CFLAGS \"${P4EST_CFLAGS}\")

set(P4EST_CPPFLAGS \"\")

set(P4EST_FFLAGS "${CMAKE_Fortran_FLAGS}\ ${MPI_Fortran_COMPILE_OPTIONS}")
set(P4EST_FFLAGS \"${P4EST_FFLAGS}\")

set(P4EST_FLIBS \"${MPI_Fortran_LIBRARIES}\")

set(P4EST_LDFLAGS \"${MPI_C_LINK_FLAGS}\")
set(P4EST_LIBS \"${LAPACK_LIBRARIES}\ ${BLAS_LIBRARIES}\ ${ZLIB_LIBRARIES}\ m\")

set(P4EST_ENABLE_BUILD_2D 1 CACHE BOOL "p4est is always used")
set(P4EST_ENABLE_BUILD_3D ${enable_p8est} CACHE BOOL "p8est support" FORCE)
set(P4EST_ENABLE_BUILD_P6EST ${enable_p6est} CACHE BOOL "p6est support" FORCE)

set(P4EST_ENABLE_FILE_DEPRECATED ${enable-file-deprecated})

set(P4EST_ENABLE_MEMALIGN 1)

if(P4EST_ENABLE_MPI)
  set(P4EST_ENABLE_MPICOMMSHARED ${SC_ENABLE_MPICOMMSHARED} CACHE BOOL "P4EST MPI shared memory support" FORCE)
  set(P4EST_ENABLE_MPITHREAD ${SC_ENABLE_MPITHREAD} CACHE BOOL "P4EST MPI thread support" FORCE)
  set(P4EST_ENABLE_MPIWINSHARED ${SC_ENABLE_MPIWINSHARED} CACHE BOOL "P4EST MPI window shared memory support" FORCE)
endif()

set(P4EST_NEED_M ${SC_NEED_M} CACHE BOOL "P4EST needs math -lm")

check_include_file(arpa/inet.h P4EST_HAVE_ARPA_INET_H)
check_include_file(netinet/in.h P4EST_HAVE_NETINET_IN_H)
if(WIN32 AND NOT P4EST_HAVE_ARPA_INET_H AND NOT P4EST_HAVE_NETINET_IN_H)
  check_include_file(Winsock2.h P4EST_HAVE_WINSOCK2_H)
  set(WINSOCK_LIBRARIES wsock32 ws2_32) # Iphlpapi
endif()

set(P4EST_HAVE_INTTYPES_H ${SC_HAVE_INTTYPES_H} CACHE BOOL "platform has inttypes.h")

check_symbol_exists(pthread_create pthread.h HAVE_LPTHREAD)
check_symbol_exists(lua_createtable lua.h HAVE_LUA)

set(P4EST_HAVE_MEMORY_H ${SC_HAVE_MEMORY_H} CACHE BOOL "platform has memory.h")

set(P4EST_HAVE_POSIX_MEMALIGN ${SC_HAVE_POSIX_MEMALIGN} CACHE BOOL "platform has posix_memalign")
set(P4EST_HAVE_STDINT_H ${SC_HAVE_STDINT_H} CACHE BOOL "platform has stdint.h")
set(P4EST_HAVE_STDLIB_H ${SC_HAVE_STDLIB_H} CACHE BOOL "platform has stdlib.h")

check_include_file(strings.h P4EST_HAVE_STRINGS_H)
set(P4EST_HAVE_STRING_H ${SC_HAVE_STRING_H} CACHE BOOL "platform has string.h")
set(P4EST_HAVE_SYS_STAT_H ${SC_HAVE_SYS_STAT_H} CACHE BOOL "platform has sys/stat.h")
set(P4EST_HAVE_SYS_TYPES_H ${SC_HAVE_SYS_TYPES_H} CACHE BOOL "platform has sys/types.h")

set(P4EST_HAVE_UNISTD_H ${SC_HAVE_UNISTD_H} CACHE BOOL "platform has unistd.h")
set(P4EST_HAVE_FSYNC ${SC_HAVE_FSYNC} CACHE BOOL "platform has fsync")
set(P4EST_HAVE_GETOPT_H ${SC_HAVE_GETOPT_H} CACHE BOOL "platform has getopt.h")

if(P4EST_HAVE_ZLIB)
  # ZLIB::ZLIB would be defined in sc/cmake/zlib.cmake IMPORTED INTERFACE GLOBAL via libsc
  if(NOT TARGET ZLIB::ZLIB)
    # libsc didn't build Zlib, so must find an existing Zlib
    find_package(ZLIB REQUIRED)
  endif()
  set(CMAKE_REQUIRED_LIBRARIES ZLIB::ZLIB)
  check_symbol_exists(adler32_combine zlib.h P4EST_HAVE_ZLIB)
  if(vtk_binary)
    set(P4EST_ENABLE_VTK_COMPRESSION 1 CACHE BOOL "p4est VTK compression support")
  endif()
endif()

if(CMAKE_BUILD_TYPE MATCHES "Debug")
  set(P4EST_ENABLE_DEBUG 1)
endif()

configure_file(${CMAKE_CURRENT_LIST_DIR}/p4est_config.h.in ${PROJECT_BINARY_DIR}/include/p4est_config.h)
