
# This file is part of the p4est library
# Makefile.am test
# included non-recursively from toplevel directory

p4est_test_programs = \
        test/p4est_test_comm test/p4est_test_pqueue test/p4est_test_hash \
        test/p4est_test_file test/p4est_test_quadrants \
        test/p4est_test_balance test/p4est_test_partition \
        test/p4est_test_order

noinst_PROGRAMS += \
        test/p4est_bug_wp

check_PROGRAMS += $(p4est_test_programs)

test_p4est_test_comm_SOURCES = test/test_comm.c
test_p4est_test_pqueue_SOURCES = test/test_pqueue.c
test_p4est_test_hash_SOURCES = test/test_hash.c
test_p4est_test_file_SOURCES = test/test_file.c
test_p4est_test_quadrants_SOURCES = test/test_quadrants.c
test_p4est_test_balance_SOURCES = test/test_balance.c
test_p4est_test_partition_SOURCES = test/test_partition.c
test_p4est_test_order_SOURCES = test/test_order.c

test_p4est_bug_wp_SOURCES = test/bug_weighted_partition.c

TESTS += $(p4est_test_programs)

LINT_CSOURCES += \
        $(test_p4est_test_comm_SOURCES) \
        $(test_p4est_test_pqueue_SOURCES) \
        $(test_p4est_test_hash_SOURCES) \
        $(test_p4est_test_file_SOURCES) \
        $(test_p4est_test_quadrants_SOURCES) \
        $(test_p4est_test_balance_SOURCES) \
        $(test_p4est_test_partition_SOURCES) \
        $(test_p4est_test_order_SOURCES) \
        $(test_p4est_bug_wp_SOURCES)
