
# This file is part of p4est.
# Makefile.am in example/particles
# included non-recursively from toplevel directory

if P4EST_ENABLE_BUILD_2D
bin_PROGRAMS += \
  example/particles/p4est_particles
example_particles_p4est_particles_SOURCES = \
  example/particles/particles_global.h \
  example/particles/particles2.c
endif

if P4EST_ENABLE_BUILD_3D
bin_PROGRAMS += \
  example/particles/p8est_particles
example_particles_p8est_particles_SOURCES = \
  example/particles/particles_global.h \
  example/particles/particles3.c
endif

## from example/steps for future reference
##
## EXTRA_DIST +=
##
## dist_p4estdata_DATA +=  \
##   example/particles/hole_2d_cubit.inp \
##   example/particles/hole_2d_cubit.jou \
##   example/particles/hole_3d_cubit.inp \
##   example/particles/hole_3d_cubit.jou \
##   example/particles/hole_2d_gmsh.inp \
##   example/particles/hole_2d_gmsh.geo \
##   example/particles/hole_3d_gmsh.inp \
##   example/particles/hole_3d_gmsh.geo
