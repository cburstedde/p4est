
# This file is part of the p4est library
# Makefile.am in src
# included non-recursively from toplevel directory

libp4est_internal_headers =
libp4est_generated_headers = src/p4est_config.h
libp4est_installed_headers = \
        src/p4est_types.h \
        src/p4est_base.h \
        src/p4est_connectivity.h src/p4est.h \
        src/p4est_algorithms.h src/p4est_communication.h \
        src/p4est_file.h src/p4est_vtk.h src/p4est_mesh.h
libp4est_compiled_sources = \
        src/p4est_base.c \
        src/p4est_connectivity.c src/p4est.c \
        src/p4est_algorithms.c src/p4est_communication.c \
        src/p4est_file.c src/p4est_vtk.c src/p4est_mesh.c

lib_LIBRARIES += src/libp4est.a
src_libp4est_a_SOURCES = \
        $(libp4est_internal_headers) \
        $(libp4est_installed_headers) \
        $(libp4est_compiled_sources)

src_libp4est_a_CPPFLAGS = \
        $(AM_CPPFLAGS) -I@top_srcdir@/sc/libb64

nodist_include_HEADERS += $(libp4est_generated_headers)
include_HEADERS += $(libp4est_installed_headers)

AM_CPPFLAGS += \
        -I@top_srcdir@/src -I@top_builddir@/sc/src -I@top_srcdir@/sc/src

LDADD += @top_builddir@/src/libp4est.a @top_builddir@/sc/src/libsc.a

LINT_CSOURCES += $(libp4est_compiled_sources)
