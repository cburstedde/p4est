
# This file is part of the p4est library
# Makefile.am test
# included non-recursively from toplevel directory

p4est_test_programs =
if P4EST_BUILD_2D
p4est_test_programs += \
        test/p4est_test_comm test/p4est_test_hash \
        test/p4est_test_file test/p4est_test_quadrants \
        test/p4est_test_balance test/p4est_test_partition \
        test/p4est_test_order test/p4est_test_coarsen \
        test/p4est_test_valid test/p4est_test_balance_type \
        test/p4est_test_loadsave test/p4est_test_iterate \
        test/p4est_test_lnodes
noinst_PROGRAMS += \
        test/p4est_bug_wp
endif
if P4EST_BUILD_3D
p4est_test_programs += \
        test/p8est_test_quadrants test/p8est_test_balance \
        test/p8est_test_partition test/p8est_test_coarsen \
        test/p8est_test_valid test/p8est_test_balance_type \
        test/p8est_test_face_transform test/p8est_test_edge_face_corners \
        test/p8est_test_periodic test/p8est_test_geometry \
        test/p8est_test_loadsave test/p8est_test_brick \
        test/p8est_test_iterate test/p8est_test_lnodes
endif

check_PROGRAMS += $(p4est_test_programs)

test_p4est_test_comm_SOURCES = test/test_comm.c
test_p4est_test_hash_SOURCES = test/test_hash.c
test_p4est_test_file_SOURCES = test/test_file.c
test_p4est_test_quadrants_SOURCES = test/test_quadrants2.c
test_p4est_test_balance_SOURCES = test/test_balance2.c
test_p4est_test_partition_SOURCES = test/test_partition2.c
test_p4est_test_order_SOURCES = test/test_order.c
test_p4est_test_coarsen_SOURCES = test/test_coarsen2.c
test_p4est_test_valid_SOURCES = test/test_valid2.c
test_p4est_test_balance_type_SOURCES = test/test_balance_type2.c
test_p4est_test_loadsave_SOURCES = test/test_loadsave2.c
test_p4est_test_iterate_SOURCES = test/test_iterate2.c
test_p4est_test_lnodes_SOURCES = test/test_lnodes2.c

test_p4est_bug_wp_SOURCES = test/bug_weighted_partition.c

test_p8est_test_quadrants_SOURCES = test/test_quadrants3.c
test_p8est_test_balance_SOURCES = test/test_balance3.c
test_p8est_test_partition_SOURCES = test/test_partition3.c
test_p8est_test_coarsen_SOURCES = test/test_coarsen3.c
test_p8est_test_valid_SOURCES = test/test_valid3.c
test_p8est_test_balance_type_SOURCES = test/test_balance_type3.c
test_p8est_test_face_transform_SOURCES = test/test_face_transform3.c
test_p8est_test_edge_face_corners_SOURCES = test/test_edge_face_corners3.c
test_p8est_test_periodic_SOURCES = test/test_periodic3.c
test_p8est_test_geometry_SOURCES = test/test_geometry3.c
test_p8est_test_loadsave_SOURCES = test/test_loadsave3.c
test_p8est_test_brick_SOURCES = test/test_brick3.c
test_p8est_test_iterate_SOURCES = test/test_iterate3.c
test_p8est_test_lnodes_SOURCES = test/test_lnodes3.c

TESTS += $(p4est_test_programs)

LINT_CSOURCES += \
        $(test_p4est_test_comm_SOURCES) \
        $(test_p4est_test_hash_SOURCES) \
        $(test_p4est_test_file_SOURCES) \
        $(test_p4est_test_quadrants_SOURCES) \
        $(test_p4est_test_balance_SOURCES) \
        $(test_p4est_test_partition_SOURCES) \
        $(test_p4est_test_order_SOURCES) \
        $(test_p4est_test_coarsen_SOURCES) \
        $(test_p4est_test_valid_SOURCES) \
        $(test_p4est_test_balance_type_SOURCES) \
        $(test_p4est_test_loadsave_SOURCES) \
        $(test_p4est_test_iterate_SOURCES) \
        $(test_p4est_test_lnodes_SOURCES) \
        $(test_p4est_bug_wp_SOURCES) \
        $(test_p8est_test_quadrants_SOURCES) \
        $(test_p8est_test_balance_SOURCES) \
        $(test_p8est_test_partition_SOURCES) \
        $(test_p8est_test_coarsen_SOURCES) \
        $(test_p8est_test_valid_SOURCES) \
        $(test_p8est_test_balance_type_SOURCES) \
        $(test_p8est_test_face_transform_SOURCES) \
        $(test_p8est_test_edge_face_corners_SOURCES) \
        $(test_p8est_test_periodic_SOURCES) \
        $(test_p8est_test_geometry_SOURCES) \
        $(test_p8est_test_loadsave_SOURCES) \
        $(test_p8est_test_brick_SOURCES) \
        $(test_p8est_test_iterate_SOURCES) \
        $(test_p8est_test_lnodes_SOURCES)
