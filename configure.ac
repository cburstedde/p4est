dnl
dnl p4est configure.ac
dnl

AC_PREREQ([2.59])
AC_INIT([p4est], [0.1], [ccgo@ices.utexas.edu])
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_LIBOBJ_DIR(src)
AC_CONFIG_SRCDIR([src/p4est.h])
AC_CONFIG_HEADER([src/p4est_config.h])
AM_INIT_AUTOMAKE(1.8.2)

# Configure switches
P4EST_ARG_ENABLE([debug], [enable debug mode], [P4EST_HAVE_DEBUG])
P4EST_ARG_ENABLE([vtk-doubles], [use doubles for vtk file data],
                 [P4EST_VTK_DOUBLES])
P4EST_ARG_DISABLE([vtk-binary], [write vtk ascii file data],
                 [P4EST_VTK_BINARY])
P4EST_ARG_DISABLE([vtk-zlib], [disable zlib compression for vtk binary data],
                 [P4EST_VTK_COMPRESSION])

# General settings.
AC_CANONICAL_HOST
AC_GNU_SOURCE

# Checks for programs.
echo "o----------------------------------"
echo "| Checking for compilers"
echo "o----------------------------------"

ACX_MPI([icc gcc cc])
AC_PROG_RANLIB

# Checks for libraries.
echo "o----------------------------------"
echo "| Checking for libraries"
echo "o----------------------------------"

P4EST_REQUIRE_LIB([m], [fabs])
P4EST_REQUIRE_LIB([z], [compress2])

echo "o----------------------------------"
echo "| Checking headers and functions"
echo "o----------------------------------"

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([execinfo.h getopt.h stdint.h stdlib.h stdio.h ctype.h])

# Checks for functions.
AC_CHECK_FUNCS([backtrace backtrace_symbols])
P4EST_REQUIRE_FUNCS([memset pow sqrt floor])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T

echo "o----------------------------------"
echo "| Results are"
echo "o----------------------------------"
echo "| CC:          $CC"
echo "| CFLAGS:      $CFLAGS"
echo "| CPP:         $CPP"
echo "| CPPFLAGS:    $CPPFLAGS"
echo "| LDFLAGS:     $LDFLAGS"
echo "| LIBS:        $LIBS"
echo "o----------------------------------"

# And the end.
AC_CONFIG_FILES([Makefile
                 src/Makefile
                 src/libb64/Makefile
                 test/Makefile
                 example/Makefile
                 example/read_forest/Makefile
                 example/simple/Makefile
                 example/timings/Makefile
                 doc/Makefile
                ])

AC_OUTPUT
