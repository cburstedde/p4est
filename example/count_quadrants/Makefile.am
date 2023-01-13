
# This file is part of p4est.
# Makefile.am in example/timings
# included non-recursively from toplevel directory

if P4EST_ENABLE_BUILD_2D
bin_PROGRAMS += example/count_quadrants/p4est_count_quadrants

example_count_quadrants_p4est_count_quadrants_SOURCES = \
        example/count_quadrants/p4est_count_quadrants.c

LINT_CSOURCES += \
        $(example_count_quadrants_p4est_count_quadrants_SOURCES)
endif

if P4EST_ENABLE_BUILD_3D
bin_PROGRAMS += example/count_quadrants/p8est_count_quadrants

example_count_quadrants_p8est_count_quadrants_SOURCES = \
        example/count_quadrants/p8est_count_quadrants.c

LINT_CSOURCES += \
        $(example_count_quadrants_p8est_count_quadrants_SOURCES)
endif