
# This file is part of p4est.
# Makefile.am in example/mesh
# included non-recursively from toplevel directory

if P4EST_ENABLE_BUILD_2D
bin_PROGRAMS += example/mesh/p4est_mesh
example_mesh_p4est_mesh_SOURCES = example/mesh/mesh2.c
endif

if P4EST_ENABLE_BUILD_3D
bin_PROGRAMS += example/mesh/p8est_mesh \
	example/mesh/p8est_periodicity

example_mesh_p8est_mesh_SOURCES = example/mesh/mesh3.c
example_mesh_p8est_periodicity_SOURCES = example/mesh/periodicity3.c
endif
