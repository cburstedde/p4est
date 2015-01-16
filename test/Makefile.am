
# This file is part of p4est.
# Makefile.am test
# included non-recursively from toplevel directory

p4est_test_programs =
if P4EST_ENABLE_BUILD_2D
p4est_test_programs += \
        test/p4est_test_comm test/p4est_test_hash \
        test/p4est_test_quadrants test/p4est_test_balance \
        test/p4est_test_partition test/p4est_test_coarsen \
        test/p4est_test_valid test/p4est_test_balance_type \
        test/p4est_test_loadsave test/p4est_test_order \
        test/p4est_test_ghost \
        test/p4est_test_iterate test/p4est_test_lnodes \
        test/p4est_test_search test/p4est_test_brick \
        test/p4est_test_partition_corr \
        test/p4est_test_conn_complete test/p4est_test_balance_seeds \
        test/p4est_test_wrap test/p4est_test_replace test/p4est_test_join \
        test/p4est_test_conn_reduce test/p4est_test_plex \
        test/p4est_test_connrefine
if P4EST_WITH_METIS
p4est_test_programs += \
        test/p4est_test_reorder
endif
endif
if P4EST_ENABLE_BUILD_3D
p4est_test_programs += \
        test/p8est_test_quadrants test/p8est_test_balance \
        test/p8est_test_partition test/p8est_test_coarsen \
        test/p8est_test_valid test/p8est_test_balance_type \
        test/p8est_test_face_transform test/p8est_test_edge_face_corners \
        test/p8est_test_periodic test/p8est_test_loadsave \
        test/p8est_test_ghost \
        test/p8est_test_iterate test/p8est_test_lnodes \
        test/p8est_test_search test/p8est_test_brick \
        test/p8est_test_partition_corr \
        test/p8est_test_conn_complete test/p8est_test_balance_seeds \
        test/p8est_test_wrap test/p8est_test_replace test/p8est_test_join \
        test/p8est_test_conn_reduce test/p8est_test_plex \
        test/p8est_test_connrefine
if P4EST_WITH_METIS
p4est_test_programs += \
        test/p8est_test_reorder
endif
endif

check_PROGRAMS += $(p4est_test_programs)

test_p4est_test_comm_SOURCES = test/test_comm.c
test_p4est_test_hash_SOURCES = test/test_hash.c
test_p4est_test_quadrants_SOURCES = test/test_quadrants2.c
test_p4est_test_balance_SOURCES = test/test_balance2.c
test_p4est_test_partition_SOURCES = test/test_partition2.c
test_p4est_test_order_SOURCES = test/test_order.c
test_p4est_test_coarsen_SOURCES = test/test_coarsen2.c
test_p4est_test_valid_SOURCES = test/test_valid2.c
test_p4est_test_balance_type_SOURCES = test/test_balance_type2.c
test_p4est_test_loadsave_SOURCES = test/test_loadsave2.c
test_p4est_test_ghost_SOURCES = test/test_ghost2.c
test_p4est_test_iterate_SOURCES = test/test_iterate2.c
test_p4est_test_lnodes_SOURCES = test/test_lnodes2.c
test_p4est_test_search_SOURCES = test/test_search2.c
test_p4est_test_brick_SOURCES = test/test_brick2.c
test_p4est_test_partition_corr_SOURCES = test/test_partition_corr2.c
test_p4est_test_conn_complete_SOURCES = test/test_conn_complete2.c
test_p4est_test_balance_seeds_SOURCES = test/test_balance_seeds2.c
test_p4est_test_wrap_SOURCES = test/test_wrap2.c
test_p4est_test_replace_SOURCES = test/test_replace2.c
test_p4est_test_join_SOURCES = test/test_join2.c
test_p4est_test_conn_reduce_SOURCES = test/test_conn_reduce2.c
test_p4est_test_plex_SOURCES = test/test_plex2.c
test_p4est_test_connrefine_SOURCES = test/test_connrefine2.c

test_p8est_test_quadrants_SOURCES = test/test_quadrants3.c
test_p8est_test_balance_SOURCES = test/test_balance3.c
test_p8est_test_partition_SOURCES = test/test_partition3.c
test_p8est_test_coarsen_SOURCES = test/test_coarsen3.c
test_p8est_test_valid_SOURCES = test/test_valid3.c
test_p8est_test_balance_type_SOURCES = test/test_balance_type3.c
test_p8est_test_face_transform_SOURCES = test/test_face_transform3.c
test_p8est_test_edge_face_corners_SOURCES = test/test_edge_face_corners3.c
test_p8est_test_periodic_SOURCES = test/test_periodic3.c
test_p8est_test_loadsave_SOURCES = test/test_loadsave3.c
test_p8est_test_ghost_SOURCES = test/test_ghost3.c
test_p8est_test_brick_SOURCES = test/test_brick3.c
test_p8est_test_iterate_SOURCES = test/test_iterate3.c
test_p8est_test_lnodes_SOURCES = test/test_lnodes3.c
test_p8est_test_search_SOURCES = test/test_search3.c
test_p8est_test_partition_corr_SOURCES = test/test_partition_corr3.c
test_p8est_test_conn_complete_SOURCES = test/test_conn_complete3.c
test_p8est_test_balance_seeds_SOURCES = test/test_balance_seeds3.c
test_p8est_test_wrap_SOURCES = test/test_wrap3.c
test_p8est_test_replace_SOURCES = test/test_replace3.c
test_p8est_test_join_SOURCES = test/test_join3.c
test_p8est_test_conn_reduce_SOURCES = test/test_conn_reduce3.c
test_p8est_test_plex_SOURCES = test/test_plex3.c
test_p8est_test_connrefine_SOURCES = test/test_connrefine3.c

test_p4est_test_plex_CPPFLAGS = @P4EST_PETSC_INCLUDE_DIRS@ $(AM_CPPFLAGS) 
test_p8est_test_plex_CPPFLAGS = @P4EST_PETSC_INCLUDE_DIRS@ $(AM_CPPFLAGS) 
test_p4est_test_plex_LDADD = @P4EST_PETSC_LINK_LIBS@ $(LDADD) 
test_p8est_test_plex_LDADD = @P4EST_PETSC_LINK_LIBS@ $(LDADD) 

if P4EST_WITH_METIS
test_p4est_test_reorder_SOURCES = test/test_reorder2.c
test_p8est_test_reorder_SOURCES = test/test_reorder3.c
endif

TESTS += $(p4est_test_programs)

LINT_CSOURCES += \
        $(test_p4est_test_comm_SOURCES) \
        $(test_p4est_test_hash_SOURCES) \
        $(test_p4est_test_quadrants_SOURCES) \
        $(test_p4est_test_balance_SOURCES) \
        $(test_p4est_test_partition_SOURCES) \
        $(test_p4est_test_order_SOURCES) \
        $(test_p4est_test_coarsen_SOURCES) \
        $(test_p4est_test_valid_SOURCES) \
        $(test_p4est_test_balance_type_SOURCES) \
        $(test_p4est_test_loadsave_SOURCES) \
        $(test_p4est_test_ghost_SOURCES) \
        $(test_p4est_test_iterate_SOURCES) \
        $(test_p4est_test_lnodes_SOURCES) \
        $(test_p4est_test_search_SOURCES) \
        $(test_p4est_test_brick_SOURCES) \
        $(test_p4est_test_partition_corr_SOURCES) \
        $(test_p4est_test_reorder_SOURCES) \
        $(test_p4est_test_balance_seeds_SOURCES) \
        $(test_p4est_test_wrap_SOURCES) \
        $(test_p4est_test_replace_SOURCES) \
        $(test_p4est_test_join_SOURCES) \
        $(test_p4est_test_conn_reduce_SOURCES) \
        $(test_p4est_test_plex_SOURCES) \
        $(test_p4est_test_connrefine_SOURCES) \
        $(test_p8est_test_quadrants_SOURCES) \
        $(test_p8est_test_balance_SOURCES) \
        $(test_p8est_test_partition_SOURCES) \
        $(test_p8est_test_coarsen_SOURCES) \
        $(test_p8est_test_valid_SOURCES) \
        $(test_p8est_test_balance_type_SOURCES) \
        $(test_p8est_test_face_transform_SOURCES) \
        $(test_p8est_test_edge_face_corners_SOURCES) \
        $(test_p8est_test_periodic_SOURCES) \
        $(test_p8est_test_loadsave_SOURCES) \
        $(test_p8est_test_ghost_SOURCES) \
        $(test_p8est_test_brick_SOURCES) \
        $(test_p8est_test_iterate_SOURCES) \
        $(test_p8est_test_lnodes_SOURCES) \
        $(test_p8est_test_search_SOURCES) \
        $(test_p8est_test_reorder_SOURCES) \
        $(test_p8est_test_partition_corr_SOURCES) \
        $(test_p8est_test_balance_seeds_SOURCES) \
        $(test_p8est_test_wrap_SOURCES) \
        $(test_p8est_test_replace_SOURCES) \
        $(test_p8est_test_join_SOURCES) \
        $(test_p8est_test_conn_reduce_SOURCES) \
        $(test_p8est_test_plex_SOURCES) \
        $(test_p8est_test_connrefine_SOURCES)

if P4EST_WITH_METIS
LINT_CSOURCES += \
        $(test_p4est_test_reorder_SOURCES) \
        $(test_p8est_test_reorder_SOURCES)
endif
