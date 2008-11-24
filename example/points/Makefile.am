
# This file is part of the p4est library
# Makefile.am in example/points
# included non-recursively from toplevel directory

bin_PROGRAMS += \
        example/points/p4est_points \
        example/points/p8est_points

example_points_p4est_points_SOURCES = example/points/points2.c
example_points_p8est_points_SOURCES = example/points/points3.c

LINT_CSOURCES += \
        $(example_points_p4est_points_SOURCES) \
        $(example_points_p8est_points_SOURCES)
