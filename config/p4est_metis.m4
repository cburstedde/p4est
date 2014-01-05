
dnl P4EST_CHECK_METIS(PREXIX)
dnl Check for METIS library and link a test program
dnl
AC_DEFUN([P4EST_CHECK_METIS], [

METIS_LIBS=""
METIS_INCLUDES=""
SC_ARG_WITH_PREFIX([metis], [enable metis-dependent code], [METIS], [$1])
if test "$$1_WITH_METIS" != "no" ; then
  $1_METIS_INC=""
  $1_METIS_LD=""
  $1_METIS_LIB="-lmetis"
  if test "$$1_WITH_METIS" != "yes" ; then
    $1_METIS_INC="-I$$1_WITH_METIS/include"
    $1_METIS_LD="-L$$1_WITH_METIS/lib"
  fi
  METIS_LIBS="$$1_METIS_LD $$1_METIS_LIB"
  METIS_INCLUDES="$$1_METIS_INC"
  PRE_METIS_CFLAGS="$CFLAGS"
  CFLAGS="$CFLAGS $$1_METIS_INC"
  PRE_METIS_LIBS="$LIBS"
  LIBS="$LIBS $$1_METIS_LIB"
  PRE_METIS_LDFLAGS="$LDFLAGS"
  LDFLAGS="$LDFLAGS $$1_METIS_LD"

  AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <metis.h>]],
[int n = 0, xadj, adjncy, adjwgt, vwgt;
 int nparts = 0, options = 0, volume, part, ncon, vsize;
 float tpwgts, ubvec;

 METIS_PartGraphRecursive (&n, &ncon, &xadj, &adjncy, &vwgt, &vsize,
                           &adjwgt, &nparts, &tpwgts, &ubvec, &options,
                           &volume, &part);])],,
                      [AC_MSG_ERROR([Unable to link metis])])
  CFLAGS="$PRE_METIS_CFLAGS"
  LIBS="$PRE_METIS_LIBS"
  LDFLAGS="$PRE_METIS_LDFLAGS"
fi

AC_SUBST([$1_METIS_LIBS])
AC_SUBST([$1_METIS_INCLUDES])
])
