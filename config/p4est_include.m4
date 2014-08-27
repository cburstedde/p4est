dnl
dnl p4est_include.m4 - custom macros
dnl

dnl Documentation for macro names: brackets indicate optional arguments

dnl P4EST_ARG_ENABLE(NAME, COMMENT, TOKEN)
dnl Check for --enable/disable-NAME using shell variable P4EST_ENABLE_TOKEN
dnl If shell variable is set beforehand it overrides the option
dnl If enabled, define TOKEN to 1 and set conditional P4EST_TOKEN
dnl Default is disabled
dnl
AC_DEFUN([P4EST_ARG_ENABLE],
         [SC_ARG_ENABLE_PREFIX([$1], [$2], [$3], [P4EST])])

dnl P4EST_ARG_DISABLE(NAME, COMMENT, TOKEN)
dnl Check for --enable/disable-NAME using shell variable P4EST_ENABLE_TOKEN
dnl If shell variable is set beforehand it overrides the option
dnl If enabled, define TOKEN to 1 and set conditional P4EST_TOKEN
dnl Default is enabled
dnl
AC_DEFUN([P4EST_ARG_DISABLE],
         [SC_ARG_DISABLE_PREFIX([$1], [$2], [$3], [P4EST])])

dnl P4EST_ARG_WITH(NAME, COMMENT, TOKEN)
dnl Check for --with/without-NAME using shell variable P4EST_WITH_TOKEN
dnl If shell variable is set beforehand it overrides the option
dnl If with, define TOKEN to 1 and set conditional P4EST_TOKEN
dnl Default is without
dnl
AC_DEFUN([P4EST_ARG_WITH],
         [SC_ARG_WITH_PREFIX([$1], [$2], [$3], [P4EST])])

dnl P4EST_ARG_WITHOUT(NAME, COMMENT, TOKEN)
dnl Check for --with/without-NAME using shell variable P4EST_WITH_TOKEN
dnl If shell variable is set beforehand it overrides the option
dnl If with, define TOKEN to 1 and set conditional P4EST_TOKEN
dnl Default is with
dnl
AC_DEFUN([P4EST_ARG_WITHOUT],
         [SC_ARG_WITHOUT_PREFIX([$1], [$2], [$3], [P4EST])])

dnl P4EST_CHECK_LIBRARIES(PREFIX)
dnl This macro bundles the checks for all libraries and link tests
dnl that are required by p4est.  It can be used by other packages that
dnl link to p4est to add appropriate options to LIBS.
dnl
AC_DEFUN([P4EST_CHECK_LIBRARIES],
[
P4EST_CHECK_METIS([$1])
P4EST_CHECK_PETSC([$1])
])

dnl P4EST_AS_SUBPACKAGE(PREFIX)
dnl Call from a package that is using p4est as a subpackage.
dnl Sets PREFIX_DIST_DENY=yes if p4est is make install'd.
dnl
AC_DEFUN([P4EST_AS_SUBPACKAGE],
         [SC_ME_AS_SUBPACKAGE([$1], [m4_tolower([$1])], [P4EST], [p4est])])

dnl P4EST_FINAL_MESSAGES(PREFIX)
dnl This macro prints messages at the end of the configure run.
dnl
AC_DEFUN([P4EST_FINAL_MESSAGES],
[
if test "x$$1_ENABLE_VTK_COMPRESSION" = xyes && \
   test "x$$1_HAVE_ZLIB" != xyes ; then
AC_MSG_NOTICE([- $1 -------------------------------------------------
VTK compression is enabled, but we did not find a recent zlib.
This is OK if the following does not matter to you:
Calling p4est_vtk_write_file will abort your program.
You can fix this by compiling a working zlib and pointing LIBS to it,
or by using the p4est configure option --disable-vtk-zlib.
])
fi
])
