
# This file is part of the p4est Library
# Makefile.am in src

SUBDIRS = . libb64 iniparser

lib_LIBRARIES = libp4est.a
libp4est_a_SOURCES = \
        p4est_log.h p4est_log.c \
        p4est_types.h p4est_base.h p4est_base.c p4est_mpi_dummy.h \
        p4est_obstack.h p4est_obstack.c p4est_avl.h p4est_avl.c \
        p4est_memory.h p4est_memory.c p4est_options.h p4est_options.c \
        p4est_connectivity.h p4est_connectivity.c p4est.h p4est.c \
        p4est_algorithms.h p4est_algorithms.c \
        p4est_communication.h p4est_communication.c \
        p4est_file.h p4est_file.c p4est_vtk.h p4est_vtk.c \
        p4est_mesh.h p4est_mesh.c

include_HEADERS = \
        p4est_log.h p4est_types.h p4est_base.h p4est_mpi_dummy.h \
        p4est_obstack.h p4est_avl.h p4est_memory.h \
        p4est_connectivity.h p4est.h \
        p4est_algorithms.h p4est_communication.h \
        p4est_file.h p4est_vtk.h p4est_mesh.h

AM_CPPFLAGS = \
        -I@top_builddir@/sc/src -I@top_srcdir@/sc/src \
        -I@top_srcdir@/src/libb64 -I@top_srcdir@/src/iniparser
