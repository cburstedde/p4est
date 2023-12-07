
# This file is part of p4est.
# Makefile.am in example/gmt
# included non-recursively from toplevel directory

if P4EST_ENABLE_BUILD_2D
bin_PROGRAMS += example/gmt/p4est_gmt

example_gmt_p4est_gmt_SOURCES = \
        example/gmt/gmt2.c \
        example/gmt/gmt_models.c example/gmt/gmt_models.h

LINT_CSOURCES += \
        $(example_gmt_p4est_gmt_SOURCES)

bin_PROGRAMS += example/gmt/p4est_sphere_preprocessing

example_gmt_p4est_sphere_preprocessing_SOURCES = example/gmt/sphere_preprocessing.c

LINT_CSOURCES += \
        $(example_gmt_p4est_sphere_preprocessing_SOURCES)

endif

if P4EST_ENABLE_BUILD_3D
endif
