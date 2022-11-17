
# This file is part of p4est.
# Makefile.am in example/multi
# included non-recursively from toplevel directory

if P4EST_ENABLE_BUILD_2D
bin_PROGRAMS += example/multi/p4est_overlap
example_multi_p4est_overlap_SOURCES = example/multi/overlap2.c

LINT_CSOURCES += $(example_multi_p4est_overlap_SOURCES)
endif

if P4EST_ENABLE_BUILD_3D
bin_PROGRAMS += example/multi/p8est_overlap
example_multi_p8est_overlap_SOURCES = example/multi/overlap3.c

LINT_CSOURCES += $(example_multi_p8est_overlap_SOURCES)
endif
