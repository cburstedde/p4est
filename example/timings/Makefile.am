
# This file is part of p4est.
# Makefile.am in example/timings
# included non-recursively from toplevel directory

if P4EST_ENABLE_BUILD_2D
bin_PROGRAMS += \
        example/timings/p4est_timings \
        example/timings/p4est_bricks \
        example/timings/p4est_loadconn

example_timings_p4est_timings_SOURCES = example/timings/timings2.c
example_timings_p4est_bricks_SOURCES = example/timings/bricks2.c
example_timings_p4est_loadconn_SOURCES = example/timings/loadconn2.c
endif

if P4EST_ENABLE_BUILD_3D
bin_PROGRAMS += \
        example/timings/p8est_timings \
        example/timings/p8est_bricks \
        example/timings/p8est_loadconn \
        example/timings/p8est_tsearch

example_timings_p8est_timings_SOURCES = example/timings/timings3.c
example_timings_p8est_bricks_SOURCES = example/timings/bricks3.c
example_timings_p8est_loadconn_SOURCES = example/timings/loadconn3.c
example_timings_p8est_tsearch_SOURCES = example/timings/tsearch3.c
endif

EXTRA_DIST += example/timings/timana.awk example/timings/timana.sh
