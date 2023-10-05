
dnl P4EST_CPU_SUPPORTS_AVX2(ACTION-IF-FOUND, ACTION-IF-NOT-FOUND)
dnl This macro performs a compile and link test for the AVX2 instruction set.
dnl
AC_DEFUN([P4EST_CPU_SUPPORTS_AVX2], [
   dnl AC_REQUIRE([AC_PROG_CC])
   dnl AC_LANG_PUSH([C])
   dnl AS_VAR_PUSHDEF([simd_cv_check_cpu_init],
   dnl      [AS_TR_SH([simd_cv_cpu_supports_$1])])dnl
   BACKUP_CFLAGS="${CFLAGS}"
   CFLAGS="${BACKUP_CFLAGS} -mavx2"
   AC_RUN_IFELSE(
       [AC_LANG_PROGRAM([[
         #include <immintrin.h>
        ]],
        [__m256i a;
         a = _mm256_set1_epi32 (1);
         _mm256_abs_epi32 (a);
         return 0;
        ])],
        [$1],
        [CFLAGS="$BACKUP_CFLAGS"
         $2],
        [CFLAGS="$BACKUP_CFLAGS"
         $2]
   )
   dnl AC_LANG_POP([C])
   dnl AS_VAR_IF([simd_cv_check_cpu_init], [yes], [$3], [$4])dnl
   dnl AS_VAR_POPDEF([simd_cv_check_cpu_init])dnl
])

dnl P4EST_CHECK_AVX2(PREFIX)
dnl
dnl Check configure line for --enable/disable-avx2, default enabled.
dnl If enabled, run a compile test with an AVX2 program.
dnl Only if that succeeds, shell variable P4EST_ENABLE_AVX2 stays "yes."
dnl On success, we add -mavx2 to the CFLAGS and also
dnl set the #define and AM_CONDITIONAL P4EST_ENABLE_AVX2.
dnl
dnl The argument PREFIX should be set to P4EST for configuring p4est.
dnl If p4est is part of a bigger program, use any application prefix.
dnl
AC_DEFUN([P4EST_CHECK_AVX2], [
  AC_MSG_CHECKING([for AVX2])
  AC_ARG_ENABLE([avx2], [AS_HELP_STRING([--disable-avx2],
                         [disable AVX2 instruction test and feature])],
                , [enableval=yes])
  if test "x$enableval" != xno ; then
    P4EST_CPU_SUPPORTS_AVX2([enableval=yes], [enableval=no])
  fi
  if test "x$enableval" != xno ; then
    AC_DEFINE([ENABLE_AVX2], 1, [Define to 1 if AVX2 instructions are active])
  fi
  AM_CONDITIONAL([$1_ENABLE_AVX2], [test "x$enableval" != xno])
  $1_ENABLE_AVX2="$enableval"
  AC_MSG_RESULT([$enableval])
])
