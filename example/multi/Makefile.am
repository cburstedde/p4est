
# This file is part of p4est.
# Makefile.am in example/multi
# included non-recursively from toplevel directory

if P4EST_ENABLE_BUILD_2D
bin_PROGRAMS += \
  example/multi/p4est_overlap \
  example/multi/p4est_overset
example_multi_p4est_overlap_SOURCES = example/multi/overlap2.c
example_multi_p4est_overset_SOURCES = example/multi/overset2.c

LINT_CSOURCES += \
  $(example_multi_p4est_overlap_SOURCES) \
  $(example_multi_p4est_overset_SOURCES)
endif

if P4EST_ENABLE_BUILD_3D
bin_PROGRAMS += \
  example/multi/p8est_overlap \
  example/multi/p8est_overset
example_multi_p8est_overlap_SOURCES = example/multi/overlap3.c
example_multi_p8est_overset_SOURCES = example/multi/overset3.c

LINT_CSOURCES += \
  $(example_multi_p8est_overlap_SOURCES) \
  $(example_multi_p8est_overset_SOURCES)
endif
