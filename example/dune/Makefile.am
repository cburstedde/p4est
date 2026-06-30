
# This file is part of p4est.
# Makefile.am in example/dune
# included non-recursively from toplevel directory

if P4EST_ENABLE_BUILD_2D
bin_PROGRAMS += \
  example/dune/p4est_dune_interface
example_dune_p4est_dune_interface_SOURCES = \
  example/dune/dune_interface2.c
endif

if P4EST_ENABLE_BUILD_3D
bin_PROGRAMS += \
  example/dune/p8est_dune_interface
example_dune_p8est_dune_interface_SOURCES = \
  example/dune/dune_interface3.c
endif
