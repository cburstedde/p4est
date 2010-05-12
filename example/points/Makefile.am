
# This file is part of p4est.
# Makefile.am in example/points
# included non-recursively from toplevel directory

if P4EST_ENABLE_BUILD_2D
bin_PROGRAMS += example/points/p4est_points
example_points_p4est_points_SOURCES = example/points/points2.c

LINT_CSOURCES += $(example_points_p4est_points_SOURCES)
endif

if P4EST_ENABLE_BUILD_3D
bin_PROGRAMS += example/points/p8est_points
example_points_p8est_points_SOURCES = example/points/points3.c

LINT_CSOURCES += $(example_points_p8est_points_SOURCES)
endif
