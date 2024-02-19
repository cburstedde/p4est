
# This file is part of p4est.
# Makefile.am in example/steps
# included non-recursively from toplevel directory

if P4EST_ENABLE_BUILD_2D
bin_PROGRAMS += \
  example/steps/p4est_step1 \
  example/steps/p4est_step2 \
  example/steps/p4est_step3 \
  example/steps/p4est_step4 \
  example/steps/p4est_step5
example_steps_p4est_step1_SOURCES = example/steps/p4est_step1.c
example_steps_p4est_step2_SOURCES = example/steps/p4est_step2.c
example_steps_p4est_step3_SOURCES = example/steps/p4est_step3.c
example_steps_p4est_step4_SOURCES = example/steps/p4est_step4.c
example_steps_p4est_step5_SOURCES = example/steps/p4est_step5.c
endif

if P4EST_ENABLE_BUILD_3D
bin_PROGRAMS += \
  example/steps/p8est_step1 \
  example/steps/p8est_step2 \
  example/steps/p8est_step3 \
  example/steps/p8est_step4 \
  example/steps/p8est_step5
example_steps_p8est_step1_SOURCES = example/steps/p8est_step1.c
example_steps_p8est_step2_SOURCES = example/steps/p8est_step2.c
example_steps_p8est_step3_SOURCES = example/steps/p8est_step3.c
example_steps_p8est_step4_SOURCES = example/steps/p8est_step4.c
example_steps_p8est_step5_SOURCES = example/steps/p8est_step5.c
endif

EXTRA_DIST += \
  example/steps/hw32.h

dist_p4estdata_DATA +=  \
  example/steps/hole_2d_cubit.inp \
  example/steps/hole_2d_cubit.jou \
  example/steps/hole_3d_cubit.inp \
  example/steps/hole_3d_cubit.jou \
  example/steps/hole_2d_gmsh.inp \
  example/steps/hole_2d_gmsh.geo \
  example/steps/hole_3d_gmsh.inp \
  example/steps/hole_3d_gmsh.geo
