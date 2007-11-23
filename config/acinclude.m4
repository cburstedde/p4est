dnl
dnl P4EST acinclude.m4 - custom macros
dnl

dnl Documentation for macro names: brackets indicate optional arguments

dnl P4EST_ARG_ENABLE(NAME, COMMENT, DEFINE)
dnl Check for --enable/disable-NAME using shell variable p4est_enable_NAME
dnl If enabled, define DEFINE to 1
dnl Default is disabled
dnl
AC_DEFUN([P4EST_ARG_ENABLE],
[
AC_ARG_ENABLE([$1], [AC_HELP_STRING([--enable-$1], [$2])],
              [p4est_define_$3="$enableval"], [p4est_define_$3="no"])
if test "$p4est_define_$3" != "no" ; then
  AC_DEFINE([$3], 1, [$2])
fi
])

dnl P4EST_ARG_DISABLE(NAME, COMMENT, DEFINE)
dnl Check for --enable/disable-NAME using shell variable p4est_enable_NAME
dnl If enabled, define DEFINE to 1
dnl Default is enabled
dnl
AC_DEFUN([P4EST_ARG_DISABLE],
[
AC_ARG_ENABLE([$1], [AC_HELP_STRING([--disable-$1], [$2])],
              [p4est_define_$3="$enableval"], [p4est_define_$3="yes"])
if test "$p4est_define_$3" != "no" ; then
  AC_DEFINE([$3], 1, [don't $2])
fi
])

dnl P4EST_ARG_WITH(NAME, COMMENT, CONDITIONAL)
dnl Check for --with/without-NAME using shell variable p4est_with_NAME
dnl If with, define conditional CONDITIONAL
dnl Default is without
dnl
AC_DEFUN([P4EST_ARG_WITH],
    [AC_ARG_WITH([$1], [AC_HELP_STRING([--with-$1], [$2])],
      [p4est_with_$1=$withval], [p4est_with_$1=no])
    AM_CONDITIONAL([$3], [test "$p4est_with_$1" != "no"])])

dnl P4EST_ARG_WITHOUT(NAME, COMMENT, CONDITIONAL)
dnl Check for --with/without-NAME using shell variable p4est_with_NAME
dnl If with, define conditional CONDITIONAL
dnl Default is with
dnl
AC_DEFUN([P4EST_ARG_WITHOUT],
    [AC_ARG_WITH([$1], [AC_HELP_STRING([--without-$1], [$2])],
      [p4est_with_$1=$withval], [p4est_with_$1=yes])
    AM_CONDITIONAL([$3], [test "$p4est_with_$1" = "yes"])])

dnl P4EST_REQUIRE_LIB(LIBRARY, FUNCTION)
dnl Check for FUNCTION in LIBRARY, exit with error if not found
dnl
AC_DEFUN([P4EST_REQUIRE_LIB],
    [AC_CHECK_LIB([$1], [$2], ,
      [AC_MSG_ERROR([Could not find function $2 in library $1])])])

dnl P4EST_REQUIRE_LIB_SEARCH(LIBRARY LIST, FUNCTION)
dnl Check for FUNCTION in any of LIBRARY LIST, exit with error if not found
dnl
AC_DEFUN([P4EST_REQUIRE_LIB_SEARCH],
    [AC_SEARCH_LIBS([$2], [$1], ,
      [AC_MSG_ERROR([Could not find function $2 in any of $1])])])

dnl P4EST_REQUIRE_FUNCS(FUNCTION LIST)
dnl Check for all functions in FUNCTION LIST, exit with error if not found
dnl
AC_DEFUN([P4EST_REQUIRE_FUNCS],
    [
dnl AS_IF([false], [AC_CHECK_FUNCS([$1])]) dnl keep autoscan quiet
     AC_CHECK_FUNCS([$1])
     AC_FOREACH([p4est_thefunc], [$1],
      [AC_CHECK_FUNC([p4est_thefunc], ,
        [AC_MSG_ERROR([Could not find function p4est_thefunc])])])])

dnl EOF acinclude.m4
