
dnl Address GMT's configure logic relying on #defines we also define.
AC_DEFUN([P4EST_GMT_UNDEFINE],
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
]])

dnl P4EST_CHECK_GMT(PREFIX)
dnl Check for the GMT library and link a test program
dnl
dnl This macro checks for the GDAL library in passing.
dnl If the GDAL library is configured explicitly, the preprocessor
dnl #define PREFIX_WITH_GDAL is set to 1, otherwise it is undefined.
dnl If the GMT library is configured, it implicitly activates the
dnl default GDAL configuration if not given.  However, in this case,
dnl PREFIX_WITH_GDAL remains undefined.
dnl
dnl To check whether GDAL is configured either way, query
dnl #ifdef PREFIX_WITH_GDAL_ACTIVE.
dnl To check whether GMT is configured, query
dnl #ifdef PREFIX_WITH_GMT.
dnl
AC_DEFUN([P4EST_CHECK_GMT], [

dnl AC_MSG_CHECKING([for GDAL (required by GMT)])
dnl GMT relies on the GDAL library and includes <gdal.h>, which usually
dnl resides not in the standard include path but in a subdirectory.
dnl Thus, add a preprocessor flag if specified.
SC_ARG_WITH_PREFIX([gdal],
                   [Link to GDAL, optionally providing path to installation
                    directory to override defaults],
                   [GDAL], [$1], [[[=DIR]]])
SC_ARG_WITH_PREFIX([gmt],
                   [activate GMT, optionally providing path to installation
                    directory to override defaults],
                   [GMT], [$1], [[[=DIR]]])

# enable convenience functionality: gmt implicitly requests gdal
if test "x$$1_WITH_GMT" != xno ; then
  dnl GMT is generally requested
  if test "x$$1_WITH_GDAL" = xno ; then
    dnl GMT specified and GDAL not specified
    $1_WITH_GDAL=yes
  fi
fi

# configure linking to the GDAL library
if test "x$$1_WITH_GDAL" != xno ; then
  dnl GDAL is generally requested
  AC_DEFINE([WITH_GDAL_ACTIVE], 1, [GDAL is configured and usable])
  if test "x$$1_WITH_GDAL" = xyes ; then
    $1_GDAL_CPP="-I/usr/include/gdal"
    $1_GDAL_LDF=
  else
    $1_GDAL_CPP="-I$$1_WITH_GDAL/include/gdal"
    $1_GDAL_LDF="-L$$1_WITH_GDAL/lib"
  fi
  CPPFLAGS="$CPPFLAGS $$1_GDAL_CPP"
  LDFLAGS="$$1_GDAL_LDF $LDFLAGS"
  SC_SEARCH_LIBS([GDALCreate], [
P4EST_GMT_UNDEFINE
#pragma GCC diagnostic ignored "-Wpedantic"
#include <gdal.h>
  ], [
if (GDALCreate (NULL, "", 0, 1, 2, GDT_Float64, NULL) != NULL) {
  GDALAllRegister ();
}
  ],
  [gdal], [],
  [AC_MSG_ERROR([GDAL library requested but not able to test link.
P4EST_DOUBLE_LINE
The --with-gdal configure option was given, with or without argument or
implicitly by configuring --with-gmt.
  Without argument, we expect that the header files reside under
/usr/include/gdal or are found by the CPPFLAGS include path, and we expect
that -lgdal works with whatever the LIBS variable currently holds.
If you would like to augment these variables, please add them to the
configure line with your preferred values.  You may omit these if the
library is installed in the default system search paths.
  With --with-gdal=<path>, we add <path>/lib and <path>/include/gdal to the
search paths, which works if this is the actual directory structure of the
gdal installation.
  Please see the file config.log for the error messages of the link test.
P4EST_DOUBLE_LINE])])
fi

dnl AC_MSG_CHECKING([for GMT])

if test "x$$1_WITH_GMT" != xno ; then
  dnl GMT is generally requested
  if test "x$$1_WITH_GMT" = xyes ; then
    $1_GMT_CPP="-I/usr/include/gmt"
    $1_GMT_LDF=
  else
    $1_GMT_CPP="-I$$1_WITH_GMT/include/gmt"
    $1_GMT_LDF="-L$$1_WITH_GMT/lib"
  fi
  CPPFLAGS="$CPPFLAGS $$1_GMT_CPP"
  LDFLAGS="$$1_GMT_LDF $LDFLAGS"
  SC_SEARCH_LIBS([floatAlmostEqualUlpsAndAbs], [
P4EST_GMT_UNDEFINE
#pragma GCC diagnostic ignored "-Wpedantic"
#include <gmt_common_math.h>
#include <gmt_dev.h>
  ], [
if (floatAlmostEqualUlpsAndAbs (0., 1., 2., 17)) {
  gmt_grd_shift (NULL, NULL, 18);
}
  ],
  [gmt], [],
  [AC_MSG_ERROR([GMT library requested but not able to test link.
P4EST_DOUBLE_LINE
The --with-gmt configure option was given, with or without argument.
  Without argument, we expect that the header files reside under
/usr/include/gmt or are found by the CPPFLAGS include path, and we expect
that -lgmt works with whatever the LIBS variable currently holds.
If you would like to augment these variables, please add them to the
configure line with your preferred values.  You may omit these if the
library is installed in the default system search paths.
  With --with-gmt=<path>, we add <path>/lib and <path>/include/gmt to the
search paths, which works if this is the actual directory structure of the
gmt installation.
  Please see the file config.log for the error messages of the link test.
P4EST_DOUBLE_LINE])],  [-lm])
fi
])
