
dnl Check for METIS library and link a test program

AC_DEFUN([P4EST_CHECK_METIS], [

METIS_LIBS=""
METIS_INCLUDES=""
P4EST_ARG_WITH([metis], [enable metis-dependent code], [METIS])
if test "$P4EST_WITH_METIS" != "no" ; then
  P4EST_METIS_INC=""
  P4EST_METIS_LD=""
  P4EST_METIS_LIB="-lmetis"
  if test "$P4EST_WITH_METIS" != "yes" ; then
    P4EST_METIS_INC="-I$P4EST_WITH_METIS/include"
    P4EST_METIS_LD="-L$P4EST_WITH_METIS/lib"
  fi
  METIS_LIBS="$P4EST_METIS_LD $P4EST_METIS_LIB"
  METIS_INCLUDES="$P4EST_METIS_INC"
  PRE_METIS_CFLAGS="$CFLAGS"
  CFLAGS="$CFLAGS $P4EST_METIS_INC"
  PRE_METIS_LIBS="$LIBS"
  LIBS="$LIBS $P4EST_METIS_LIB"
  PRE_METIS_LDFLAGS="$LDFLAGS"
  LDFLAGS="$LDFLAGS $P4EST_METIS_LD"

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

AC_SUBST([METIS_LIBS])
AC_SUBST([METIS_INCLUDES])
])
