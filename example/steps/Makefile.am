
# This file is part of p4est.
# Makefile.am in example/steps
# included non-recursively from toplevel directory

if P4EST_ENABLE_BUILD_2D
bin_PROGRAMS += example/steps/p4est_step1
example_steps_p4est_step1_SOURCES = example/steps/p4est_step1.c

LINT_CSOURCES += $(example_steps_p4est_step1_SOURCES)
endif

if P4EST_ENABLE_BUILD_3D
bin_PROGRAMS += example/steps/p8est_step1
example_steps_p8est_step1_SOURCES = example/steps/p8est_step1.c

LINT_CSOURCES += $(example_steps_p8est_step1_SOURCES)
endif

EXTRA_DIST += example/steps/hw32.h
