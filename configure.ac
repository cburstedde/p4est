dnl
dnl This file is part of the p4est library.
dnl

AC_INIT([p4est],
        [m4_esyscmd([config/git-version-gen .tarball-version])],
        [ccgo@ices.utexas.edu])
AC_PREREQ(2.61)
AC_CONFIG_HEADERS([src/pre_config.h])
AC_CONFIG_SRCDIR([src/p4est.h])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_MACRO_DIR([config])
AC_CONFIG_MACRO_DIR([sc/config])
AC_PREFIX_DEFAULT([$PWD/local])
AM_INIT_AUTOMAKE([subdir-objects])
SC_PREFIX_CONFIG_H([src/p4est_config.h])

# These variables are only relevant for libsc integration case 2 (see below).
SC_ARG_OVERRIDE_LOAD([P4EST])
SC_ARG_OVERRIDE_VAR([P4EST], [SC_PACKAGE_CPPFLAGS])
SC_ARG_OVERRIDE_VAR([P4EST], [SC_PACKAGE_LDADD])

# Configure switches
P4EST_ARG_ENABLE([debug], [enable debug mode], [DEBUG])
P4EST_ARG_ENABLE([vtk-doubles], [use doubles for vtk file data],
                 [VTK_DOUBLES])
P4EST_ARG_DISABLE([vtk-binary], [write vtk ascii file data],
                  [VTK_BINARY])
P4EST_ARG_DISABLE([vtk-zlib], [disable zlib compression for vtk binary data],
                  [VTK_COMPRESSION])
P4EST_ARG_DISABLE([2d], [disable the 2D library], [BUILD_2D])
P4EST_ARG_DISABLE([3d], [disable the 3D library], [BUILD_3D])

SC_ARG_WITH_BUILTIN_ALL_PREFIX([P4EST])

# There are 4 options for integrating libsc.
# 1. No option or --with-sc=yes to use the source in subdirectory sc.
# 2. Specify --with-sc=<dir> to point at a libsc source directory.
#    Information on the build directory can be supplied through
#    P4EST_SC_PACKAGE_CPPFLAGS and P4EST_SC_PACKAGE_LDADD.
# 3. Specify --with-sc=<dir> to point at a `make install`ed libsc.
# 4. Specify --without-sc or --with-sc=no to expect libsc to be found
#    by setting suitable CPPFLAGS, LDFLAGS and/or LIBS variables.
# Only in case 1 the subdirectory sc is configured and `make dist`ed.
# Using case 2 maintainers can distribute libsc as part of another package.
# In cases 3 and 4 make dist is forbidden.  This is meant for isolated use.
#
# libsc comes with builtin getopt, obstack, zlib and lua alternatives.
# In cases 1 and 2 the configure checks from libsc are repeated by p4est.
# In case 3 the checks are omitted and Makefile.sc.mk is used from <dir>.
# In case 4 the checks are omitted.
# In all cases Makefile.p4est.mk is installed to allow easy use of p4est.
#
dnl Above logic is implemented further down, this is just the option check.
dnl Does not AC_DEFINE and AC_CONDITIONAL if --with-sc=<dir> is given.
P4EST_ARG_WITHOUT_YES([sc],
[disable builtin libsc, provide your own. Can use\
 --with-sc=<dir> with libsc installation path],
                      [PROVIDE_SC])

# Checks for programs.
echo "o----------------------------------"
echo "| Checking for compilers"
echo "o----------------------------------"

SC_MPI([P4EST], [icc gcc cc])
AM_PROG_CC_C_O
AC_PROG_RANLIB
SC_PROG_LINT
SC_C_VERSION

# Checks for libraries.
echo "o----------------------------------"
echo "| Checking for libraries"
echo "o----------------------------------"

SC_REQUIRE_LIB([m], [fabs])

# Checks for header files.
echo "o---------------------------------------"
echo "| Checking headers"
echo "o---------------------------------------"

AC_CHECK_HEADERS([arpa/inet.h netinet/in.h unistd.h])

# Checks for functions.
echo "o---------------------------------------"
echo "| Checking functions"
echo "o---------------------------------------"

AC_CHECK_FUNCS([fsync])

# Check for builtin third-party software.
echo "o---------------------------------------"
echo "| Checking third-party software"
echo "o---------------------------------------"

# Prepare make dist behavior.
P4EST_DIST_ALLOW="yes"

# Configure the SC library recursively.
P4EST_SC_SUBDIR=
P4EST_SC_MK_INCLUDE=
P4EST_SC_USE_MK="no"
P4EST_SC_ENABLE_CHECKS="no"
if test "$P4EST_WITH_PROVIDE_SC" = yes ; then
  # case 1
  AC_CONFIG_SUBDIRS([sc])
  P4EST_SC_SUBDIR=sc
  P4EST_SC_DIR="\$(top_srcdir)/sc"
  P4EST_SC_AMFLAGS="-I $P4EST_SC_DIR/config"
  P4EST_SC_CPPFLAGS="-I\$(top_builddir)/sc/src -I\$(top_srcdir)/sc/src"
  P4EST_SC_LDADD="\$(top_builddir)/sc/src/libsc.a"
  P4EST_SC_MK_INCLUDE="include \${sysconfdir}/Makefile.sc.mk"
  P4EST_SC_ENABLE_CHECKS="yes"
else
  if test "$P4EST_WITH_PROVIDE_SC" != no ; then
    # case 2 and 3
    P4EST_SC_DIR="$P4EST_WITH_PROVIDE_SC"
    SC_CHECK_PACKAGE([P4EST_SC], [true], [true], [true], [true])
    P4EST_SC_AMFLAGS="-I $P4EST_SC_CFG"
    if test "$P4EST_SC_INSTALL" = "no" ; then
      # case 2
      P4EST_SC_CPPFLAGS="$P4EST_SC_PACKAGE_CPPFLAGS -I$P4EST_SC_INC"
      P4EST_SC_LDADD="$P4EST_SC_PACKAGE_LDADD"
      P4EST_SC_ENABLE_CHECKS="yes"
    else
      # case 3
      P4EST_SC_CPPFLAGS="\$(SC_CPPFLAGS)"
      P4EST_SC_LDADD="\$(SC_LDFLAGS) \$(SC_LIBS)"
      P4EST_SC_MK_INCLUDE="include $P4EST_SC_ETC/Makefile.sc.mk"
      P4EST_SC_USE_MK="yes"
      P4EST_DIST_ALLOW="no"
    fi
  else
    # case 4
    P4EST_DIST_ALLOW="no"
  fi
fi

AC_SUBST([P4EST_SC_SUBDIR])
AC_SUBST([P4EST_SC_AMFLAGS])
AC_SUBST([P4EST_SC_CPPFLAGS])
AC_SUBST([P4EST_SC_LDADD])
AC_SUBST([P4EST_SC_MK_INCLUDE])
AM_CONDITIONAL([P4EST_SC_USE_MK], [test "$P4EST_SC_USE_MK" = "yes"])
SC_BUILTIN_ALL_PREFIX([P4EST], [test "$P4EST_SC_ENABLE_CHECKS" = "yes"])

# Implement make dist behavior
AM_CONDITIONAL([P4EST_DIST_ALLOW], [test "$P4EST_DIST_ALLOW" = "yes"])

# Checks for typedefs, structures, and compiler characteristics.
echo "o---------------------------------------"
echo "| Checking keywords and types"
echo "o---------------------------------------"

AC_C_BIGENDIAN
AC_C_CONST
AC_C_INLINE
AC_C_RESTRICT
AC_TYPE_SIZE_T
AC_TYPE_SSIZE_T

AC_DEFINE_UNQUOTED(CC,          ["${CC}"],          [C compiler])
AC_DEFINE_UNQUOTED(C_VERSION,   ["${C_VERSION}"],   [C compiler version])
AC_DEFINE_UNQUOTED(CFLAGS,      ["${CFLAGS}"],      [C compiler flags])
AC_DEFINE_UNQUOTED(CPP,         ["${CPP}"],         [C preprocessor])
AC_DEFINE_UNQUOTED(CPPFLAGS,    ["${CPPFLAGS}"],    [C preprocessor flags])
AC_DEFINE_UNQUOTED(LDFLAGS,     ["${LDFLAGS}"],     [Linker flags])
AC_DEFINE_UNQUOTED(LIBS,        ["${LIBS}"],        [Libraries])

echo "o----------------------------------"
echo "| Results for p4est are"
echo "o----------------------------------"
echo "| CC:          $CC"
echo "| C_VERSION:   $C_VERSION"
echo "| CFLAGS:      $CFLAGS"
echo "| CPP:         $CPP"
echo "| CPPFLAGS:    $CPPFLAGS"
echo "| LDFLAGS:     $LDFLAGS"
echo "| LIBS:        $LIBS"
echo "o----------------------------------"

# Create output files.
SC_ARG_OVERRIDE_SAVE([P4EST])

AC_CONFIG_FILES([Makefile Makefile.p4est.pre])
AC_OUTPUT
