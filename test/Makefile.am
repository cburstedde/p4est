
## Makefile.am in test

bin_PROGRAMS = p4est_test_quadrants p4est_test_file

p4est_test_quadrants_SOURCES = test_quadrants.c
p4est_test_file_SOURCES = test_file.c

TESTS = p4est_test_quadrants p4est_test_file

AM_CPPFLAGS = -I@top_srcdir@/src
LDADD = @top_builddir@/src/libp4est.a
