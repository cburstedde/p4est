dnl
dnl This file is part of p4est.
dnl

AC_INIT([p4est],
        [m4_esyscmd([config/git-version-gen .tarball-version])],
        [info@p4est.org])
AC_PREREQ(2.61)
AC_CONFIG_HEADERS([src/pre_config.h])
AC_CONFIG_SRCDIR([src/p4est.h])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([config])
AC_PREFIX_DEFAULT([$PWD/local])
AX_PREFIX_CONFIG_H([src/p4est_config.h])
AM_INIT_AUTOMAKE([parallel-tests subdir-objects])
AM_SILENT_RULES
SC_VERSION([P4EST])

echo "o---------------------------------------"
echo "| Checking options"
echo "o---------------------------------------"

dnl---------------------- HOW SUBPACKAGES WORK --------------------------
dnl
dnl p4est relies on libsc.  In a build from source, libsc usually resides
dnl in the subdirectory sc.  This location can be overridden by the
dnl environment variable P4EST_SC_SOURCE; this situation can arise when
dnl both p4est and libsc are subpackages to yet another software.  The path
dnl in P4EST_SC_SOURCE must be relative to the p4est toplevel directory. All
dnl of this works without specifying any configure command line option.
dnl However, if libsc is already make install'd in the system and should
dnl be used from there, use --with-sc=<path to libsc install directory>.
dnl In this case, p4est expects subdirectories etc, include, lib, and
dnl share/aclocal, which are routinely created by libsc's make install.
dnl
P4EST_ARG_WITH([sc], [path to installed libsc (optional)], [SC])

P4EST_ARG_ENABLE([debug], [enable debug mode (assertions and extra checks)],
                 [DEBUG])
P4EST_ARG_ENABLE([vtk-doubles], [use doubles for vtk file data],
                 [VTK_DOUBLES])
P4EST_ARG_DISABLE([vtk-binary], [write vtk ascii file data],
                  [VTK_BINARY])
P4EST_ARG_DISABLE([vtk-zlib], [disable zlib compression for vtk binary data],
                  [VTK_COMPRESSION])
P4EST_ARG_DISABLE([2d], [disable the 2D library], [BUILD_2D])
P4EST_ARG_DISABLE([3d], [disable the 3D library], [BUILD_3D])

echo "o---------------------------------------"
echo "| Checking MPI and related programs"
echo "o---------------------------------------"

dnl A nonempty second/third argument causes to enable F77/CXX, respectively.
SC_MPI_CONFIG([P4EST], [yes], [yes])
SC_MPI_ENGAGE([P4EST])
# This is needed for compatibility with automake >= 1.12
m4_ifdef([AM_PROG_AR],[AM_PROG_AR])
SC_PROG_LINT
SC_C_VERSION
LT_INIT

echo "o----------------------------------"
echo "| Checking for libraries"
echo "o----------------------------------"

SC_CHECK_LIBRARIES([P4EST])
P4EST_CHECK_LIBRARIES([P4EST])

echo "o---------------------------------------"
echo "| Checking headers"
echo "o---------------------------------------"

AC_CHECK_HEADERS([arpa/inet.h netinet/in.h unistd.h])

echo "o---------------------------------------"
echo "| Checking functions"
echo "o---------------------------------------"

AC_CHECK_FUNCS([fsync])

echo "o----------------------------------"
echo "| Checking subpackages"
echo "o----------------------------------"

P4EST_DIST_ALLOW=
P4EST_SC_SUBDIR=
P4EST_SC_MK_USE=

if test x$P4EST_WITH_SC != xno ; then
  AC_MSG_NOTICE([Using make installed libsc])

  # Verify that we are using a libsc installation
  P4EST_SC_DIR="$P4EST_WITH_SC"
  SC_CHECK_INSTALL([P4EST_SC], [true], [true], [true], [true])

  # Set variables for using the subpackage
  P4EST_SC_AMFLAGS="-I $P4EST_SC_CFG"
  P4EST_SC_MK_USE="yes"
  P4EST_SC_MK_INCLUDE="include $P4EST_SC_ETC/Makefile.sc.mk"
  P4EST_SC_CPPFLAGS="\$(SC_CPPFLAGS)"
  # P4EST_SC_LDADD="\$(SC_LDFLAGS) \$(SC_LIBS)"
  P4EST_SC_LDADD="\$(SC_LDFLAGS) -lsc"
else
  AC_MSG_NOTICE([Building with sc source])

  # Prepare for a build using sc sources
  P4EST_DIST_ALLOW="yes"
  if test -z "$P4EST_SC_SOURCE" ; then
    P4EST_SC_SOURCE="sc"
    P4EST_SC_SUBDIR="sc"
    AC_CONFIG_SUBDIRS([sc])
  fi
  P4EST_SC_AMFLAGS="-I \$(top_srcdir)/$P4EST_SC_SOURCE/config"
  P4EST_SC_MK_INCLUDE="include \${p4est_sysconfdir}/Makefile.sc.mk"
  P4EST_SC_CPPFLAGS="-I\$(top_builddir)/$P4EST_SC_SOURCE/src \
-I\$(top_srcdir)/$P4EST_SC_SOURCE/src"
  P4EST_SC_LDADD="\$(top_builddir)/$P4EST_SC_SOURCE/src/libsc.la"
fi

# Prepare make dist behavior and dependencies.
dnl TODO fix library dependencies
dnl P4EST_DEPENDENCIES=
dnl AC_SUBST([P4EST_DEPENDENCIES])

dnl Make sure we find the m4 macros provided by libsc
AC_SUBST([P4EST_SC_AMFLAGS])

dnl We call make in this subdirectory if not empty
AC_SUBST([P4EST_SC_SUBDIR])

dnl We will need these variables to compile and link with libsc
AM_CONDITIONAL([P4EST_SC_MK_USE], [test -n "$P4EST_SC_MK_USE"])
AC_SUBST([P4EST_SC_MK_INCLUDE])
AC_SUBST([P4EST_SC_CPPFLAGS])
AC_SUBST([P4EST_SC_LDADD])

# Implement make dist behavior and dependencies.
AM_CONDITIONAL([P4EST_DIST_ALLOW], [test "$P4EST_DIST_ALLOW" = "yes"])

# Print summary.

AC_DEFINE_UNQUOTED(CPP,         ["${CPP}"],         [C preprocessor])
AC_DEFINE_UNQUOTED(CPPFLAGS,    ["${CPPFLAGS}"],    [C preprocessor flags])
AC_DEFINE_UNQUOTED(CC,          ["${CC}"],          [C compiler])
dnl AC_DEFINE_UNQUOTED(C_VERSION,   ["${C_VERSION}"],   [C compiler version])
AC_DEFINE_UNQUOTED(CFLAGS,      ["${CFLAGS}"],      [C compiler flags])
AC_DEFINE_UNQUOTED(LDFLAGS,     ["${LDFLAGS}"],     [Linker flags])
AC_DEFINE_UNQUOTED(LIBS,        ["${LIBS}"],        [Libraries])

echo "o----------------------------------"
echo "| Results for p4est are"
echo "o----------------------------------"
echo "| CPP:         $CPP"
echo "| CPPFLAGS:    $CPPFLAGS"
echo "| CC:          $CC"
dnl echo "| C_VERSION:   $C_VERSION"
echo "| CFLAGS:      $CFLAGS"
echo "| LDFLAGS:     $LDFLAGS"
echo "| LIBS:        $LIBS"
echo "o----------------------------------"

# Create output files.
AC_CONFIG_FILES([Makefile Makefile.p4est.pre])
AC_OUTPUT

# Final messages.
SC_FINAL_MESSAGES([P4EST])
