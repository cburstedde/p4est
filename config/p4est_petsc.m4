
dnl P4EST_CHECK_PETSC(PREFIX)
dnl Check for the PETSc library and link a test program
dnl
AC_DEFUN([P4EST_CHECK_PETSC], [

AC_MSG_CHECKING([for PETSc])

SC_ARG_WITH_PREFIX([petsc], [enable petsc-dependent code], [PETSC], [$1])
if test "x$$1_WITH_PETSC" != xno ; then
  $1_PETSC_INC=
  $1_PETSC_LD=
  $1_PETSC_LIB="-lpetsc"
  if test "x$$1_WITH_PETSC" != xyes ; then
    $1_PETSC_INC="-I$$1_WITH_PETSC/include"
    $1_PETSC_LD="-L$$1_WITH_PETSC/lib"
  fi
  PRE_PETSC_CPPFLAGS="$CPPFLAGS"
  CPPFLAGS="$CPPFLAGS $$1_PETSC_INC"
  PRE_PETSC_LDFLAGS="$LDFLAGS"
  LDFLAGS="$LDFLAGS $$1_PETSC_LD"
  PRE_PETSC_LIBS="$LIBS"
  LIBS="$$1_PETSC_LIB $LIBS"

  AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <petsc.h>]],
[[
 PetscErrorCode ierr;

 ierr = PetscInitialize (NULL, NULL, NULL, NULL);CHKERRQ(ierr);
 ierr = PetscFinalize();
 return 0;
]])],,
                 [AC_MSG_ERROR([Unable to link petsc])])
dnl Keep the variables changed as done above
dnl CPPFLAGS="$PRE_PETSC_CPPFLAGS"
dnl LDFLAGS="$PRE_PETSC_LDFLAGS"
dnl LIBS="$PRE_PETSC_LIBS"

  AC_MSG_RESULT([successful])
else
  AC_MSG_RESULT([not used])
fi

dnl No AC_SUBST since we're changing variables directly
dnl AC_SUBST([$1_PETSC_LIBS])
dnl AC_SUBST([$1_PETSC_INCLUDES])
])

