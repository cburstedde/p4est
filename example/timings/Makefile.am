
# This file is part of the p4est library
# Makefile.am in example/timings
# included non-recursively from toplevel directory

if P4EST_BUILD_2D
bin_PROGRAMS += example/timings/p4est_timings
example_timings_p4est_timings_SOURCES = example/timings/timings2.c

LINT_CSOURCES += $(example_timings_p4est_timings_SOURCES)
endif

if P4EST_BUILD_3D
bin_PROGRAMS += example/timings/p8est_timings
example_timings_p8est_timings_SOURCES = example/timings/timings3.c

LINT_CSOURCES += $(example_timings_p8est_timings_SOURCES)
endif
