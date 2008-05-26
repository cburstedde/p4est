
## Makefile.am in example/second

bin_PROGRAMS = p4est_second
p4est_second_SOURCES = second.c

AM_CPPFLAGS = \
        -I@top_builddir@/sc/src -I@top_srcdir@/sc/src \
        -I@top_srcdir@/src

LDADD = @top_builddir@/src/libp4est.a @top_builddir@/src/libb64/libb64.a
