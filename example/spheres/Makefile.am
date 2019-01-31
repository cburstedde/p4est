
# This file is part of p4est.
# Makefile.am in example/spheres
# included non-recursively from toplevel directory

if P4EST_ENABLE_BUILD_2D
bin_PROGRAMS += \
  example/spheres/p4est_spheres
example_spheres_p4est_spheres_SOURCES = \
  example/spheres/spheres_global.h \
  example/spheres/spheres2.c

LINT_CSOURCES += $(example_spheres_p4est_spheres_SOURCES)
endif

if P4EST_ENABLE_BUILD_3D
bin_PROGRAMS += \
  example/spheres/p8est_spheres
example_spheres_p8est_spheres_SOURCES = \
  example/spheres/spheres_global.h \
  example/spheres/spheres3.c

LINT_CSOURCES += $(example_spheres_p8est_spheres_SOURCES)
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
