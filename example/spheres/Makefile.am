
# This file is part of p4est.
# Makefile.am in example/spheres
# included non-recursively from toplevel directory

if P4EST_ENABLE_BUILD_2D

# Some of the sphere logic is generic enough for the library
libp4est_installed_headers += \
  example/spheres/p4est_spheres.h
libp4est_compiled_sources += \
  example/spheres/p4est_spheres.c

# Specific code for this example program
bin_PROGRAMS += \
  example/spheres/p4est_spheres
example_spheres_p4est_spheres_SOURCES = \
  example/spheres/spheres_global.h \
  example/spheres/spheres2.c
endif

if P4EST_ENABLE_BUILD_3D

# Some of the sphere logic is generic enough for the library
libp4est_installed_headers += \
  example/spheres/p8est_spheres.h example/spheres/p4est_to_p8est_spheres.h
libp4est_compiled_sources += \
  example/spheres/p8est_spheres.c

# Specific code for this example program
bin_PROGRAMS += \
  example/spheres/p8est_spheres
example_spheres_p8est_spheres_SOURCES = \
  example/spheres/spheres_global.h \
  example/spheres/spheres3.c
endif

## from example/steps for future reference
##
## EXTRA_DIST +=
##
## dist_p4estdata_DATA +=  \
##   example/spheres/hole_2d_cubit.inp \
##   example/spheres/hole_2d_cubit.jou \
##   example/spheres/hole_3d_cubit.inp \
##   example/spheres/hole_3d_cubit.jou \
##   example/spheres/hole_2d_gmsh.inp \
##   example/spheres/hole_2d_gmsh.geo \
##   example/spheres/hole_3d_gmsh.inp \
##   example/spheres/hole_3d_gmsh.geo
