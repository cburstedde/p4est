
# This file is part of p4est.
# Makefile.am in example/quadrant_counter
# included non-recursively from toplevel directory

if P4EST_ENABLE_BUILD_2D
bin_PROGRAMS += example/quadrant_counter/p4est_quad_count
example_quadrant_counter_p4est_quad_count_SOURCES = example/quadrant_counter/quad_count2.c

LINT_CSOURCES += $(example_quadrant_counter_p4est_quad_count_SOURCES)
endif

if P4EST_ENABLE_BUILD_3D
bin_PROGRAMS += example/quadrant_counter/p8est_quad_count
example_quadrant_counter_p8est_quad_count_SOURCES = example/quadrant_counter/quad_count3.c

LINT_CSOURCES += $(example_quadrant_counter_p8est_quad_count_SOURCES)
endif