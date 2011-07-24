
# This file is part of p4est.
# Makefile.am in example/timings
# included non-recursively from toplevel directory

if P4EST_ENABLE_BUILD_2D
bin_PROGRAMS += example/balance_test/p4est_balance_test

example_balance_test_p4est_balance_test_SOURCES = \
        example/balance_test/balance_test2.c

LINT_CSOURCES += \
        $(example_balance_test_p4est_balance_test_SOURCES)
endif

if P4EST_ENABLE_BUILD_3D
bin_PROGRAMS += example/balance_test/p8est_balance_test

example_balance_test_p8est_balance_test_SOURCES = \
        example/balance_test/balance_test3.c

LINT_CSOURCES += \
        $(example_balance_test_p8est_balance_test_SOURCES)
endif
