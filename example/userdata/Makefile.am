
# This file is part of p4est.
# Makefile.am in example/userdata
# included non-recursively from toplevel directory

if P4EST_ENABLE_BUILD_2D
bin_PROGRAMS += example/userdata/p4est_userdata
example_userdata_p4est_userdata_SOURCES = \
  example/userdata/userdata2.c example/userdata/userdata_global.h
endif

if P4EST_ENABLE_BUILD_3D
bin_PROGRAMS += example/userdata/p8est_userdata
example_userdata_p8est_userdata_SOURCES = \
  example/userdata/userdata3.c example/userdata/userdata_global.h
endif
