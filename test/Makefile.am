
## Makefile.am in test

noinst_PROGRAMS = \
	p4est_bug_wp

check_PROGRAMS = \
	p4est_test_comm p4est_test_pqueue p4est_test_hash \
	p4est_test_file p4est_test_quadrants p4est_test_balance \
	p4est_test_partition p4est_test_order

p4est_test_comm_SOURCES = test_comm.c
p4est_test_pqueue_SOURCES = test_pqueue.c
p4est_test_hash_SOURCES = test_hash.c
p4est_test_file_SOURCES = test_file.c
p4est_test_quadrants_SOURCES = test_quadrants.c
p4est_test_balance_SOURCES = test_balance.c
p4est_test_partition_SOURCES = test_partition.c
p4est_test_order_SOURCES = test_order.c

p4est_bug_wp_SOURCES = bug_weighted_partition.c

TESTS = \
	p4est_test_comm p4est_test_pqueue p4est_test_hash \
	p4est_test_file p4est_test_quadrants p4est_test_balance \
	p4est_test_partition p4est_test_order

AM_CPPFLAGS = \
        -I@top_builddir@/sc/src -I@top_srcdir@/sc/src \
        -I@top_srcdir@/src

LDADD = @top_builddir@/src/libp4est.a
