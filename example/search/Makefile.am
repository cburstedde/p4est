
# This file is part of p4est.
# Makefile.am in example/timings
# included non-recursively from toplevel directory

if P4EST_ENABLE_BUILD_2D
bin_PROGRAMS += example/search/p4est_count_quadrants

example_search_p4est_count_quadrants_SOURCES = \
        example/search/p4est_count_quadrants.c

LINT_CSOURCES += \
        $(example_search_p4est_count_quadrants_SOURCES)
endif

if P4EST_ENABLE_BUILD_3D
bin_PROGRAMS += example/search/p8est_count_quadrants

example_search_p8est_count_quadrants_SOURCES = \
        example/search/p8est_count_quadrants.c

LINT_CSOURCES += \
        $(example_search_p8est_count_quadrants_SOURCES)
endif