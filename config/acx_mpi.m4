dnl This started off as a modified version of the Teuchos config dir
dnl from Trilinos with the following license.
dnl
dnl ***********************************************************************
dnl
dnl                    Teuchos: Common Tools Package
dnl                 Copyright (2004) Sandia Corporation
dnl
dnl Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
dnl license for use of this work by or on behalf of the U.S. Government.
dnl
dnl This library is free software; you can redistribute it and/or modify
dnl it under the terms of the GNU Lesser General Public License as
dnl published by the Free Software Foundation; either version 2.1 of the
dnl License, or (at your option) any later version.
dnl
dnl This library is distributed in the hope that it will be useful, but
dnl WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
dnl Lesser General Public License for more details.
dnl
dnl You should have received a copy of the GNU Lesser General Public
dnl License along with this library; if not, write to the Free Software
dnl Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
dnl USA
dnl Questions? Contact Michael A. Heroux (maherou@sandia.gov)
dnl
dnl ***********************************************************************
dnl
dnl Now almost nothing of the original code is left.
dnl
dnl @synopsis ACX_MPI([list of non-MPI C compilers],
dnl                   [list of non-MPI F77 compilers])
dnl
dnl This macro calls AC_PROG_CC and AC_PROG_F77.
dnl
dnl --enable-mpi           Turn on MPI, use compilers mpicc and mpif77.
dnl --with-mpicc=<...>     Specify MPI C compiler.
dnl --without-mpicc        Do not use a special MPI C compiler.
dnl --with-mpif77=<...>    Specify MPI F77 compiler.
dnl --without-mpif77       Do not use a special MPI F77 compiler.
dnl
dnl All --with* options turn on MPI, so --enable-mpi is then not needed.
dnl
dnl Order of precedence for the selection of MPI compilers
dnl 1. environment variables MPI_CC, MPI_F77
dnl 2. --with-mpicc=<...>, --with-mpif77=<...>
dnl 3. environment variables CC, F77
dnl 4. mpicc, mpif77 unless --without-mpicc, --without-mpif77
dnl
dnl If MPI is turned on HAVE_MPI will be defined for autoconf/automake
dnl and HAVE_MPI will be defined in the config header file.

dnl ACX_MPI_CONFIG
dnl Figure out the MPI configuration
dnl
AC_DEFUN([ACX_MPI_CONFIG],
[
HAVE_PKG_MPI=no
MPI_CC_NONE=
MPI_F77_NONE=

AC_ARG_ENABLE([mpi],
[AC_HELP_STRING([--enable-mpi], [enable MPI])],
[
if test "$enableval" = yes ; then
  HAVE_PKG_MPI=yes
elif test "$enableval" != no ; then
  AC_MSG_ERROR([Please use --enable-mpi without an argument])
fi
])

AC_ARG_WITH([mpicc],
[AC_HELP_STRING([--with-mpicc=MPICC], [specify MPI C compiler])],
[
if test "$withval" = yes ; then
  AC_MSG_ERROR([Please use --with-mpicc=MPICC with a valid mpi C compiler])
else
  HAVE_PKG_MPI=yes
  if test "$withval" = no ; then
    MPI_CC_NONE=yes
  else
    AC_CHECK_PROG([MPI_CC], [$withval], [$withval], [false])
  fi
fi
])

AC_ARG_WITH([mpif77],
[AC_HELP_STRING([--with-mpif77=MPIF77], [specify MPI F77 compiler])],
[
if test "$withval" = yes ; then
  AC_MSG_ERROR([Please use --with-mpif77=MPIF77 with a valid mpi F77 compiler])
else
  HAVE_PKG_MPI=yes
  if test "$withval" = no ; then
    MPI_F77_NONE=yes
  else
    AC_CHECK_PROG([MPI_F77], [$withval], [$withval], [false])
  fi
fi
])

AC_MSG_CHECKING([whether we are using MPI])
AC_MSG_RESULT([$HAVE_PKG_MPI])

if test "$HAVE_PKG_MPI" = yes ; then
  if test -z "$MPI_CC" -a -z "$MPI_CC_NONE" ; then
    MPI_CC=mpicc
  fi
  if test -z "$MPI_F77" -a -z "$MPI_F77_NONE" ; then
    MPI_F77=mpif77
  fi

  if test -n "$MPI_CC" ; then
    CC="$MPI_CC"
  fi
  if test -n "$MPI_F77" ; then
    F77="$MPI_F77"
  fi

  echo "                             MPI_CC set to $MPI_CC"
  echo "                            MPI_F77 set to $MPI_F77"

  AC_DEFINE([HAVE_MPI], 1, [Define to 1 if we are using MPI])
else
  unset MPI_CC
  unset MPI_F77
fi

AM_CONDITIONAL(HAVE_MPI, [test "$HAVE_PKG_MPI" = yes])
])

dnl ACX_MPI_C_COMPILE_AND_LINK([action-if-successful], [action-if-failed])
dnl Compile and link an MPI test program
dnl
AC_DEFUN([ACX_MPI_C_COMPILE_AND_LINK],
[
AC_MSG_CHECKING([compile/link for MPI C program])
AC_LINK_IFELSE([AC_LANG_PROGRAM(
[[#include <mpi.h>]], [[
MPI_Init ((int *) 0, (char ***) 0);
MPI_Finalize ();
]])],
[AC_MSG_RESULT([successful])
 $1],
[AC_MSG_RESULT([failed])
 $2])
])

dnl ACX_MPI_F77_COMPILE_AND_LINK([action-if-successful], [action-if-failed])
dnl Compile and link an MPI test program
dnl
AC_DEFUN([ACX_MPI_F77_COMPILE_AND_LINK],
[
AC_MSG_CHECKING([compile/link for MPI F77 program])
AC_LANG_PUSH([Fortran 77])
AC_LINK_IFELSE([AC_LANG_PROGRAM( ,
[[
      include 'mpif.h'
      integer ierror
      call MPI_INIT(ierror)
      call MPI_FINALIZE(ierror)
]])],
[AC_MSG_RESULT([successful])
 $1],
[AC_MSG_RESULT([failed])
 $2])
AC_LANG_POP
])

dnl ACX_MPI_VERIFY
dnl Verify MPI configuration (creating C and F77 test programs)
dnl
AC_DEFUN([ACX_MPI_VERIFY],
[
if test "$HAVE_PKG_MPI" = yes ; then
  ACX_MPI_C_COMPILE_AND_LINK( , [AC_MSG_ERROR([MPI C test failed])])
  ACX_MPI_F77_COMPILE_AND_LINK( , [AC_MSG_ERROR([MPI F77 test failed])])
fi
])

dnl ACX_MPI
dnl Configure MPI in one line
dnl
AC_DEFUN([ACX_MPI],
[
ACX_MPI_CONFIG
AC_PROG_CC([$1])
AC_PROG_F77([$2])
ACX_MPI_VERIFY
])
