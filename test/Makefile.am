
## Makefile.am in test

check_PROGRAMS = \
	p4est_test_pqueue p4est_test_hash \
	p4est_test_file p4est_test_quadrants p4est_test_balance

p4est_test_pqueue_SOURCES = test_pqueue.c
p4est_test_hash_SOURCES = test_hash.c

p4est_test_file_SOURCES = test_file.c
p4est_test_quadrants_SOURCES = test_quadrants.c
p4est_test_balance_SOURCES = test_balance.c

TESTS = \
	p4est_test_pqueue p4est_test_hash \
	p4est_test_file p4est_test_quadrants p4est_test_balance

AM_CPPFLAGS = -I@top_srcdir@/src
LDADD = @top_builddir@/src/libp4est.a
