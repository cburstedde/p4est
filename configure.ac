dnl
dnl This file is part of p4est.
dnl

AC_INIT([p4est],
        [m4_esyscmd([build-aux/git-version-gen .tarball-version])],
        [p4est@ins.uni-bonn.de])
AC_PREREQ(2.61)
AC_CONFIG_HEADERS([config/p4est_autotools_config.h])
AC_CONFIG_SRCDIR([src/p4est.h])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([config])
AC_PREFIX_DEFAULT([$PWD/local])
AX_PREFIX_CONFIG_H([config/p4est_config.h])
AM_INIT_AUTOMAKE([parallel-tests subdir-objects])
AM_SILENT_RULES
AM_EXTRA_RECURSIVE_TARGETS([justlibs doxygen])
PKG_INSTALLDIR
SC_VERSION([P4EST])

echo "o---------------------------------------"
echo "| Checking options"
echo "o---------------------------------------"

P4EST_ARG_ENABLE([debug], [enable debug mode (assertions and extra checks)],
                 [DEBUG])
P4EST_ARG_ENABLE([vtk-doubles], [use doubles for vtk file data],
                 [VTK_DOUBLES])
P4EST_ARG_ENABLE([file-deprecated], [use depreacted data file format],
                 [FILE_DEPRECATED])
P4EST_ARG_DISABLE([vtk-binary], [write vtk ascii file data],
                  [VTK_BINARY])
P4EST_ARG_DISABLE([vtk-zlib], [disable zlib compression for vtk binary data],
                  [VTK_COMPRESSION])
P4EST_ARG_DISABLE([2d], [disable the 2D library], [BUILD_2D])
P4EST_ARG_DISABLE([3d], [disable the 3D library], [BUILD_3D])
P4EST_ARG_DISABLE([p6est], [disable hybrid 2D+1D p6est library], [BUILD_P6EST])

echo "o---------------------------------------"
echo "| Checking MPI and related programs"
echo "o---------------------------------------"

dnl A nonempty second/third argument causes to enable F77+F90/CXX, respectively.
SC_MPI_CONFIG([P4EST], [], [])
SC_MPI_ENGAGE([P4EST])
# This is needed for compatibility with automake >= 1.12
m4_ifdef([AM_PROG_AR],[AM_PROG_AR])
LT_INIT

echo "o---------------------------------------"
echo "| Checking libraries"
echo "o---------------------------------------"

SC_CHECK_LIBRARIES([P4EST])
P4EST_CHECK_LIBRARIES([P4EST])

echo "o---------------------------------------"
echo "| Checking headers"
echo "o---------------------------------------"

AC_CHECK_HEADERS([arpa/inet.h netinet/in.h unistd.h])

echo "o---------------------------------------"
echo "| Checking functions"
echo "o---------------------------------------"

echo "o---------------------------------------"
echo "| Checking subpackages"
echo "o---------------------------------------"

P4EST_DIST_DENY=
P4EST_DISTCLEAN=
SC_AS_SUBPACKAGE([P4EST])
AM_CONDITIONAL([P4EST_DIST_DENY], [test "x$P4EST_DIST_DENY" != x])
AC_SUBST([P4EST_DISTCLEAN])

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
echo "| TESTS:       $P4EST_MPIRUN $P4EST_MPI_TEST_FLAGS"
echo "o----------------------------------"

# Create output files.
AC_CONFIG_FILES([Makefile doc/Doxyfile])
AC_CONFIG_FILES([config/Makefile.p4est.pre config/p4est_autotools.pc])
AC_OUTPUT

# Final messages.
SC_FINAL_MESSAGES([P4EST])
P4EST_FINAL_MESSAGES([P4EST])
