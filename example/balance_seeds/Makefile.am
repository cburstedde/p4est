
# This file is part of p4est.
# Makefile.am in example/timings
# included non-recursively from toplevel directory

if P4EST_ENABLE_BUILD_2D
bin_PROGRAMS += example/balance_seeds/p4est_balance_seeds

example_balance_seeds_p4est_balance_seeds_SOURCES = \
        example/balance_seeds/balance_seeds2.c

LINT_CSOURCES += \
        $(example_balance_seeds_p4est_balance_seeds_SOURCES)
endif

if P4EST_ENABLE_BUILD_3D
bin_PROGRAMS += example/balance_seeds/p8est_balance_seeds

example_balance_seeds_p8est_balance_seeds_SOURCES = \
        example/balance_seeds/balance_seeds3.c

LINT_CSOURCES += \
        $(example_balance_seeds_p8est_balance_seeds_SOURCES)
endif
