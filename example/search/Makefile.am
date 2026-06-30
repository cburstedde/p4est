
# This file is part of p4est.
# Makefile.am in example/search
# included non-recursively from toplevel directory

if P4EST_ENABLE_BUILD_2D
bin_PROGRAMS += \
        example/search/p4est_count_quadrants \
        example/search/p4est_search_replicated

example_search_p4est_count_quadrants_SOURCES = \
        example/search/count_quadrants2.c
example_search_p4est_search_replicated_SOURCES = \
        example/search/search_replicated2.c
endif

if P4EST_ENABLE_BUILD_3D
bin_PROGRAMS += \
        example/search/p8est_count_quadrants \
        example/search/p8est_search_replicated

example_search_p8est_count_quadrants_SOURCES = \
        example/search/count_quadrants3.c
example_search_p8est_search_replicated_SOURCES = \
        example/search/search_replicated3.c
endif
