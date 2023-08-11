
dnl P4EST_CHECK_GMT(PREFIX)
dnl Check for the GMT library and link a test program
dnl
AC_DEFUN([P4EST_CHECK_GMT], [

AC_MSG_CHECKING([for GMT])

SC_ARG_WITH_PREFIX([gmt], [enable gmt-dependent examples], [GMT], [$1])
if test "x$$1_WITH_GMT" != xno ; then
  $1_GMT_INC="-I/usr/include/gdal -I/usr/include/gmt"
  $1_GMT_LD=
  $1_GMT_LIB="-lgmt"
  if test "x$$1_WITH_GMT" != xyes ; then
    $1_GMT_INC="-I$$1_WITH_GMT/include"
    $1_GMT_LD="-L$$1_WITH_GMT/lib"
  fi
  PRE_GMT_CPPFLAGS="$CPPFLAGS"
  CPPFLAGS="$CPPFLAGS $$1_GMT_INC"
  PRE_GMT_LDFLAGS="$LDFLAGS"
  LDFLAGS="$LDFLAGS $$1_GMT_LD"
  PRE_GMT_LIBS="$LIBS"
  LIBS="$$1_GMT_LIB $LIBS"

  AC_LINK_IFELSE([AC_LANG_PROGRAM(
[[
#undef GMT
#ifdef HAVE_ZLIB
#undef HAVE_ZLIB
#endif
#undef PACKAGE_BUGREPORT
#undef PACKAGE_NAME
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_VERSION
#include <gmt_common_math.h>
#include <gmt_dev.h>
]],
[[
  if (floatAlmostEqualUlpsAndAbs (0., 1., 2., 17)) {
    gmt_grd_shift (NULL, NULL, 18);
  }
]])],,
                 [AC_MSG_ERROR([unable to link GMT:
P4EST_DOUBLE_LINE
The --with-gmt configure option was given, with or without argument.
  Without argument, we expect that CPPFLAGS leads to the header files,
and we expect that -lgmt works with whatever the LIBS variable holds.
If you would like to augment these variables, please add them to the 
configure line with your preferred values.  You may omit these if the
library is installed in the default system search paths.
  With --with-gmt=<path>, we add <path>/lib and <path>/include to search,
which works if this is the actual directory structure of the gmt install.
  Please see the file config.log for the error messages of the link test.
P4EST_DOUBLE_LINE])])
dnl Keep the variables changed as done above
dnl CPPFLAGS="$PRE_GMT_CPPFLAGS"
dnl LDFLAGS="$PRE_GMT_LDFLAGS"
dnl LIBS="$PRE_GMT_LIBS"

  AC_MSG_RESULT([successful])
else
  AC_MSG_RESULT([not used])
fi

dnl No AC_SUBST since we're changing variables directly
dnl AC_SUBST([$1_GMT_LIBS])
dnl AC_SUBST([$1_GMT_INCLUDES])
])
