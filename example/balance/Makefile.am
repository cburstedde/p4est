
# This file is part of p4est.
# Makefile.am in example/balance
# included non-recursively from toplevel directory

if P4EST_ENABLE_BUILD_2D
bin_PROGRAMS += example/balance/p4est_balance_seeds
bin_PROGRAMS += example/balance/p4est_balance_corner

example_balance_p4est_balance_seeds_SOURCES = \
        example/balance/balance_seeds2.c
example_balance_p4est_balance_corner_SOURCES = \
        example/balance/balance_corner2.c
endif

if P4EST_ENABLE_BUILD_3D
bin_PROGRAMS += example/balance/p8est_balance_seeds
bin_PROGRAMS += example/balance/p8est_balance_corner

example_balance_p8est_balance_seeds_SOURCES = \
        example/balance/balance_seeds3.c
example_balance_p8est_balance_corner_SOURCES = \
        example/balance/balance_corner3.c
endif
