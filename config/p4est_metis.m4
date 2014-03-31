
dnl P4EST_CHECK_METIS(PREFIX)
dnl Check for the METIS library and link a test program
dnl
AC_DEFUN([P4EST_CHECK_METIS], [

AC_MSG_CHECKING([for metis])

SC_ARG_WITH_PREFIX([metis], [enable metis-dependent code], [METIS], [$1])
if test "x$$1_WITH_METIS" != xno ; then
  $1_METIS_INC=
  $1_METIS_LD=
  $1_METIS_LIB="-lmetis"
  if test "x$$1_WITH_METIS" != xyes ; then
    $1_METIS_INC="-I$$1_WITH_METIS/include"
    $1_METIS_LD="-L$$1_WITH_METIS/lib"
  fi
  PRE_METIS_CPPFLAGS="$CPPFLAGS"
  CPPFLAGS="$CPPFLAGS $$1_METIS_INC"
  PRE_METIS_LDFLAGS="$LDFLAGS"
  LDFLAGS="$LDFLAGS $$1_METIS_LD"
  PRE_METIS_LIBS="$LIBS"
  LIBS="$$1_METIS_LIB $LIBS"

  AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <metis.h>]],
[[
 int n = 0, xadj, adjncy, adjwgt, vwgt;
 int nparts = 0, options = 0, volume, part, ncon, vsize;
 float tpwgts, ubvec;

 METIS_PartGraphRecursive (&n, &ncon, &xadj, &adjncy, &vwgt, &vsize,
                           &adjwgt, &nparts, &tpwgts, &ubvec, &options,
                           &volume, &part);
]])],,
                 [AC_MSG_ERROR([Unable to link metis])])
dnl Keep the variables changed as done above
dnl CPPFLAGS="$PRE_METIS_CPPFLAGS"
dnl LDFLAGS="$PRE_METIS_LDFLAGS"
dnl LIBS="$PRE_METIS_LIBS"

  AC_MSG_RESULT([successful])
else
  AC_MSG_RESULT([not used])
fi

dnl No AC_SUBST since we're changing variables directly
dnl AC_SUBST([$1_METIS_LIBS])
dnl AC_SUBST([$1_METIS_INCLUDES])
])
