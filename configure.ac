dnl
dnl This file is part of the p4est library.
dnl

AC_INIT([p4est], [0.2], [ccgo@ices.utexas.edu])
AC_PREREQ([2.59])
AC_CONFIG_HEADER([src/pre_config.h])
AC_CONFIG_SRCDIR([src/p4est.h])
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_MACRO_DIR(config)
AC_CONFIG_MACRO_DIR(sc/config)
AC_PREFIX_DEFAULT([$PWD/install])
AM_INIT_AUTOMAKE([subdir-objects])
AX_PREFIX_CONFIG_H([src/p4est_config.h])

# Configure switches
P4EST_ARG_ENABLE([debug], [enable debug mode], [DEBUG])
P4EST_ARG_ENABLE([log-level], [set the log priority [=PRIORITY]],
                 [LOG_LEVEL])
P4EST_ARG_ENABLE([vtk-doubles], [use doubles for vtk file data],
                 [VTK_DOUBLES])
P4EST_ARG_DISABLE([vtk-binary], [write vtk ascii file data],
                  [VTK_BINARY])
P4EST_ARG_DISABLE([vtk-zlib], [disable zlib compression for vtk binary data],
                  [VTK_COMPRESSION])
P4EST_ARG_DISABLE([resize-realloc], [avoid realloc in array resize],
                  [RESIZE_REALLOC])
P4EST_ARG_DISABLE([2d], [disable the 2D library], [BUILD_2D])
P4EST_ARG_DISABLE([3d], [disable the 3D library], [BUILD_3D])
P4EST_ARG_WITHOUT([sc], [disable builtin libsc, provide your own], [PROVIDE_SC])
P4EST_ARG_WITHOUT([zlib], [disable zlib in builtin libsc, provide your own],
                  [PROVIDE_ZLIB])

# Configure the SC library recursively.
if test "$P4EST_WITH_PROVIDE_SC" = yes; then
  AC_CONFIG_SUBDIRS([sc])
fi

# Checks for programs.
echo "o----------------------------------"
echo "| Checking for compilers"
echo "o----------------------------------"

ACX_MPI([P4EST], [icc gcc cc])
AM_PROG_CC_C_O
AC_PROG_RANLIB
ACX_WITH_LINT

# Checks for libraries.
echo "o----------------------------------"
echo "| Checking for libraries"
echo "o----------------------------------"

SC_REQUIRE_LIB([m], [fabs])
if test "$P4EST_WITH_PROVIDE_ZLIB" = "no" ; then
  SC_REQUIRE_LIB([z], [adler32_combine])
fi

echo "o----------------------------------"
echo "| Checking headers and functions"
echo "o----------------------------------"

# Checks for header files.
AC_CHECK_HEADERS([execinfo.h signal.h])
AC_CHECK_HEADERS([arpa/inet.h netinet/in.h])

# Checks for functions.
AC_CHECK_FUNCS([backtrace backtrace_symbols])
SC_REQUIRE_FUNCS([memset pow sqrt floor])

# Checks for typedefs, structures, and compiler characteristics.
dnl These are done from SC already.
dnl AC_C_CONST
dnl AC_C_INLINE
dnl AC_TYPE_SIZE_T
dnl AC_TYPE_SSIZE_T

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

# Create output files.
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
