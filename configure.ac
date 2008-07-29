dnl
dnl This file is part of the p4est library.
dnl

AC_INIT([p4est],[0.3.pre],[ccgo@ices.utexas.edu])
AC_PREREQ(2.61)
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
P4EST_ARG_ENABLE([stats], [enable statistics], [STATS])
P4EST_ARG_ENABLE([vtk-doubles], [use doubles for vtk file data],
                 [VTK_DOUBLES])
P4EST_ARG_DISABLE([vtk-binary], [write vtk ascii file data],
                  [VTK_BINARY])
P4EST_ARG_DISABLE([vtk-zlib], [disable zlib compression for vtk binary data],
                  [VTK_COMPRESSION])
P4EST_ARG_DISABLE([2d], [disable the 2D library], [BUILD_2D])
P4EST_ARG_DISABLE([3d], [disable the 3D library], [BUILD_3D])

dnl does not AC_DEFINE and AC_CONDITIONAL if --with-sc=<dir> is given
P4EST_ARG_WITHOUT_YES([sc], [\
disable builtin libsc, provide your own. Can use\
 --with-sc=<dir> for <dir>/include, <dir>/lib/libsc.a],
                  [PROVIDE_SC])
P4EST_ARG_WITHOUT([zlib], [disable zlib in builtin libsc, provide your own],
                  [PROVIDE_ZLIB])

# Configure the SC library recursively.
if test "$P4EST_WITH_PROVIDE_SC" = yes ; then
  AC_CONFIG_SUBDIRS([sc])
  P4EST_SC_CPPFLAGS="-I\$(top_builddir)/sc/src -I\$(top_srcdir)/sc/src"
  P4EST_SC_LDADD="\$(top_builddir)/sc/src/libsc.a"
else
  if test "$P4EST_WITH_PROVIDE_SC" != no ; then
    P4EST_SC_DIR="$P4EST_WITH_PROVIDE_SC"
    if test ! -d "$P4EST_SC_DIR" ; then
      AC_MSG_ERROR([Directory "$P4EST_SC_DIR" does not exist])
    fi
    if test ! -f "$P4EST_SC_DIR/include/sc_zlib.h" ; then
      P4EST_WITH_PROVIDE_ZLIB=no
    fi
    P4EST_SC_CPPFLAGS="-I$P4EST_SC_DIR/include"
    P4EST_SC_LDADD="-L$P4EST_SC_DIR/lib -lsc"
  else
    P4EST_WITH_PROVIDE_ZLIB=no
  fi
fi
AC_SUBST([P4EST_SC_CPPFLAGS])
AC_SUBST([P4EST_SC_LDADD])

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

# Checks for header files.
echo "o---------------------------------------"
echo "| Checking headers"
echo "o---------------------------------------"

AC_CHECK_HEADERS([arpa/inet.h netinet/in.h])

# Checks for typedefs, structures, and compiler characteristics.
dnl echo "o---------------------------------------"
dnl echo "| Checking keywords and types"
dnl echo "o---------------------------------------"
dnl
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
