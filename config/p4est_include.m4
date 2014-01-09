dnl
dnl p4est_include.m4 - custom macros
dnl

dnl Documentation for macro names: brackets indicate optional arguments

dnl P4EST_ARG_ENABLE(NAME, COMMENT, TOKEN)
dnl Check for --enable/disable-NAME using shell variable P4EST_ENABLE_TOKEN
dnl If shell variable is set beforehand it overrides the option
dnl If enabled, define TOKEN to 1 and set conditional P4EST_TOKEN
dnl Default is disabled
dnl
AC_DEFUN([P4EST_ARG_ENABLE],
         [SC_ARG_ENABLE_PREFIX([$1], [$2], [$3], [P4EST])])

dnl P4EST_ARG_DISABLE(NAME, COMMENT, TOKEN)
dnl Check for --enable/disable-NAME using shell variable P4EST_ENABLE_TOKEN
dnl If shell variable is set beforehand it overrides the option
dnl If enabled, define TOKEN to 1 and set conditional P4EST_TOKEN
dnl Default is enabled
dnl
AC_DEFUN([P4EST_ARG_DISABLE],
         [SC_ARG_DISABLE_PREFIX([$1], [$2], [$3], [P4EST])])

dnl P4EST_ARG_WITH(NAME, COMMENT, TOKEN)
dnl Check for --with/without-NAME using shell variable P4EST_WITH_TOKEN
dnl If shell variable is set beforehand it overrides the option
dnl If with, define TOKEN to 1 and set conditional P4EST_TOKEN
dnl Default is without
dnl
AC_DEFUN([P4EST_ARG_WITH],
         [SC_ARG_WITH_PREFIX([$1], [$2], [$3], [P4EST])])

dnl P4EST_ARG_WITHOUT(NAME, COMMENT, TOKEN)
dnl Check for --with/without-NAME using shell variable P4EST_WITH_TOKEN
dnl If shell variable is set beforehand it overrides the option
dnl If with, define TOKEN to 1 and set conditional P4EST_TOKEN
dnl Default is with
dnl
AC_DEFUN([P4EST_ARG_WITHOUT],
         [SC_ARG_WITHOUT_PREFIX([$1], [$2], [$3], [P4EST])])

dnl P4EST_ARG_WITH_YES(NAME, COMMENT, TOKEN)
dnl Check for --with/without-NAME using shell variable P4EST_WITH_TOKEN
dnl If shell variable is set beforehand it overrides the option
dnl If with = yes, define TOKEN to 1 and set conditional P4EST_TOKEN
dnl Default is no
dnl
AC_DEFUN([P4EST_ARG_WITH_YES],
         [SC_ARG_WITH_YES_PREFIX([$1], [$2], [$3], [P4EST])])

dnl P4EST_ARG_WITHOUT_YES(NAME, COMMENT, TOKEN)
dnl Check for --with/without-NAME using shell variable P4EST_WITH_TOKEN
dnl If shell variable is set beforehand it overrides the option
dnl If with = yes, define TOKEN to 1 and set conditional P4EST_TOKEN
dnl Default is yes
dnl
AC_DEFUN([P4EST_ARG_WITHOUT_YES],
         [SC_ARG_WITHOUT_YES_PREFIX([$1], [$2], [$3], [P4EST])])

dnl P4EST_AS_SUBPACKAGE(PREFIX, prefix)
dnl Call from a package that is using p4est as a subpackage.
dnl Sets PREFIX_DIST_DENY="yes" if p4est is make install'd.
dnl
AC_DEFUN([P4EST_AS_SUBPACKAGE],
[
$1_P4EST_SUBDIR=
$1_P4EST_MK_USE=
$1_DISTCLEAN="$$1_DISTCLEAN $1_P4EST_SOURCE.log"

SC_ARG_WITH_PREFIX([p4est], [path to installed p4est (optional)],
                   [P4EST], [$1])
if test x$$1_WITH_P4EST != xno ; then
  AC_MSG_NOTICE([Using make installed p4est])

  # Verify that we are using a p4est installation
  $1_DIST_DENY="yes"
  $1_P4EST_DIR="$$1_WITH_P4EST"
  SC_CHECK_INSTALL([$1_P4EST], [true], [true], [true], [true])

  # Set variables for using the subpackage
  $1_P4EST_AMFLAGS="-I $$1_P4EST_CFG"
  $1_P4EST_MK_USE="yes"
  $1_P4EST_MK_INCLUDE="include $$1_P4EST_ETC/Makefile.p4est.mk"
  $1_P4EST_CPPFLAGS="\$(P4EST_CPPFLAGS)"
  $1_P4EST_LDADD="\$(P4EST_LDFLAGS) -lp4est"
else
  AC_MSG_NOTICE([Building with p4est source])

  # Prepare for a build using p4est sources
  if test -z "$$1_P4EST_SOURCE" ; then
    if test -f "$1_P4EST_SOURCE.log" ; then
      $1_P4EST_SOURCE=`cat $1_P4EST_SOURCE.log`
    else
      $1_P4EST_SOURCE="p4est"
      $1_P4EST_SUBDIR="p4est"
      AC_CONFIG_SUBDIRS([p4est])
    fi
  else
    AC_CONFIG_COMMANDS([$1_P4EST_SOURCE.log],
                       [echo "$$1_P4EST_SOURCE" >$1_P4EST_SOURCE.log])
  fi
  $1_P4EST_AMFLAGS="-I \$(top_srcdir)/$$1_P4EST_SOURCE/config"
  $1_P4EST_MK_INCLUDE="include \${$2_sysconfdir}/Makefile.p4est.mk"
  $1_P4EST_CPPFLAGS="-I\$(top_builddir)/$$1_P4EST_SOURCE/src \
-I\$(top_srcdir)/$$1_P4EST_SOURCE/src"
  $1_P4EST_LDADD="\$(top_builddir)/$$1_P4EST_SOURCE/src/libp4est.la"
fi

dnl Make sure we find the m4 macros provided by p4est
AC_SUBST([$1_P4EST_AMFLAGS])

dnl We call make in this subdirectory if not empty
AC_SUBST([$1_P4EST_SUBDIR])

dnl We will need these variables to compile and link with p4est
AM_CONDITIONAL([$1_P4EST_MK_USE], [test -n "$$1_P4EST_MK_USE"])
AC_SUBST([$1_P4EST_MK_INCLUDE])
AC_SUBST([$1_P4EST_CPPFLAGS])
AC_SUBST([$1_P4EST_LDADD])
])

dnl P4EST_CHECK_LIBRARIES(PREFIX)
dnl This macro bundles the checks for all libraries and link tests
dnl that are required by p4est.  It can be used by other packages that
dnl link to p4est to add appropriate options to LIBS.
dnl
AC_DEFUN([P4EST_CHECK_LIBRARIES],
[
P4EST_CHECK_METIS([$1])
])

dnl P4EST_FINAL_MESSAGES(PREFIX)
dnl This macro prints messages at the end of the configure run.
dnl
AC_DEFUN([P4EST_FINAL_MESSAGES],)
