
dnl P4EST_CHECK_PETSC(PREFIX)
dnl Check for the PETSc library and link a test program
dnl portions adapted from libmesh petsc.m4
dnl
AC_DEFUN([P4EST_CHECK_PETSC], [

AC_MSG_CHECKING([for PETSc])

SC_ARG_WITH_PREFIX([petsc], [enable PETSc-dependent code], [PETSC], [$1])
$1_PETSC_INCLUDE_DIRS=
$1_PETSC_LINK_LIBS=
if test "x$$1_WITH_PETSC" != xno ; then
  # use the PETSC_DIR environment variable by default
  $1_PETSC_DIR="$PETSC_DIR"
  $1_PETSC_ARCH="$PETSC_ARCH"
  if test "x$$1_WITH_PETSC" != xyes ; then
    $1_PETSC_DIR="$$1_WITH_PETSC"
  fi
  if (test "x$$1_PETSC_DIR" = x); then
    AC_PATH_PROG(PETSCARCH, petscarch)
    if (test "x$PETSCARCH" != x); then
      $1_PETSC_DIR=/usr/lib/petsc
    fi
  fi
  if test ! -r $$1_PETSC_DIR/include/petscversion.h ; then
    AC_MSG_ERROR([Unable to find readable petscversion.h])
  fi
  $1_PETSC_MAJOR=`grep "define PETSC_VERSION_MAJOR" $$1_PETSC_DIR/include/petscversion.h | sed -e "s/#define PETSC_VERSION_MAJOR[ ]*//g"`
  if test "$$1_PETSC_MAJOR" -lt 3 ; then
    AC_MSG_ERROR([PETSc version >= 3.0 required])
  fi
  if test -r $$1_PETSC_DIR/makefile ; then
    $1_PETSC_LINK_LIBS=`make -s -C $$1_PETSC_DIR getlinklibs`
    $1_PETSC_INCLUDE_DIRS=`make -s -C $$1_PETSC_DIR getincludedirs`
  elif test -r $$1_PETSC_DIR/conf/variables ; then
    if ! test -r $$1_PETSC_DIR/conf/rules ; then
      AC_MSG_ERROR([Unable to find $$1_PETSC_DIR/makefile or $$1_PETSC_DIR/conf/rules])
    fi
    cat <<EOF >Makefile_config_petsc
include $$1_PETSC_DIR/conf/variables
include $$1_PETSC_DIR/conf/rules
EOF
    $1_PETSC_LINK_LIBS=`make -s -f Makefile_config_petsc getlinklibs`
    $1_PETSC_INCLUDE_DIRS=`make -s -f Makefile_config_petsc getincludedirs`
    rm -f Makefile_config_petsc
  else 
    AC_MSG_ERROR([Unable to find $$1_PETSC_DIR/makefile or $$1_PETSC_DIR/conf/variables])
  fi
  PRE_PETSC_CPPFLAGS="$CPPFLAGS"
  CPPFLAGS="$CPPFLAGS $$1_PETSC_INCLUDE_DIRS"
  PRE_PETSC_LIBS="$LIBS"
  LIBS="$LIBS $$1_PETSC_LINK_LIBS"

  AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <petsc.h>]],
[[
  PetscErrorCode ierr;

  ierr = PetscInitialize (NULL, NULL, NULL, NULL);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return 0;
]])],,
                 [AC_MSG_ERROR([Unable to link petsc])])
  CPPFLAGS="$PRE_PETSC_CPPFLAGS"
  LIBS="$PRE_PETSC_LIBS"

  AC_MSG_RESULT([successful])
else
  AC_MSG_RESULT([not used])
fi

AC_SUBST([$1_PETSC_INCLUDE_DIRS])
AC_SUBST([$1_PETSC_LINK_LIBS])
])

