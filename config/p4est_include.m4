dnl
dnl P4EST acinclude.m4 - custom macros
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

dnl P4EST_ARG_WITH_YES(NAME, COMMENT, TOKEN)
dnl Check for --with/without-NAME using shell variable P4EST_WITH_TOKEN
dnl If shell variable is set beforehand it overrides the option
dnl If with = yes, define TOKEN to 1 and set conditional P4EST_TOKEN
dnl Default is no
dnl
AC_DEFUN([P4EST_ARG_WITH_YES],
         [SC_ARG_WITH_YES_PREFIX([$1], [$2], [$3], [P4EST])])

dnl P4EST_ARG_WITHOUT_YES(NAME, COMMENT, TOKEN)
dnl Check for --with/without-NAME using shell variable P4EST_WITH_TOKEN
dnl If shell variable is set beforehand it overrides the option
dnl If with = yes, define TOKEN to 1 and set conditional P4EST_TOKEN
dnl Default is yes
dnl
AC_DEFUN([P4EST_ARG_WITHOUT_YES],
         [SC_ARG_WITHOUT_YES_PREFIX([$1], [$2], [$3], [P4EST])])

dnl P4EST_CHECK_LIBRARIES(PREFIX)
dnl This macro bundles the checks for all libraries and link tests
dnl that are required by p4est.  It can be used by other packages that
dnl link to p4est to add appropriate options to LIBS.
dnl
AC_DEFUN([P4EST_CHECK_LIBRARIES],
[
P4EST_CHECK_METIS([$1])
])
