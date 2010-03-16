
# This file is part of the p4est library
# Makefile.am in src
# included non-recursively from toplevel directory

libp4est_generated_headers = src/p4est_config.h
libp4est_installed_headers = src/p4est_base.h
libp4est_internal_headers =
libp4est_compiled_sources = src/p4est_base.c
if P4EST_ENABLE_BUILD_2D
libp4est_installed_headers += \
        src/p4est_connectivity.h src/p4est.h \
        src/p4est_bits.h src/p4est_search.h \
        src/p4est_algorithms.h src/p4est_communication.h \
        src/p4est_ghost.h src/p4est_nodes.h src/p4est_vtk.h \
        src/p4est_points.h src/p4est_geometry.h \
        src/p4est_iterate.h src/p4est_lnodes.h
libp4est_compiled_sources += \
        src/p4est_connectivity.c src/p4est.c \
        src/p4est_bits.c src/p4est_search.c \
        src/p4est_algorithms.c src/p4est_communication.c \
        src/p4est_ghost.c src/p4est_nodes.c src/p4est_vtk.c \
        src/p4est_points.c \
        src/p4est_iterate.c src/p4est_lnodes.c
endif
if P4EST_ENABLE_BUILD_3D
libp4est_installed_headers += \
        src/p4est_to_p8est.h \
        src/p8est_connectivity.h src/p8est.h \
        src/p8est_bits.h src/p8est_search.h \
        src/p8est_algorithms.h src/p8est_communication.h \
        src/p8est_ghost.h src/p8est_nodes.h src/p8est_vtk.h \
        src/p8est_trilinear.h \
        src/p8est_points.h src/p8est_geometry.h \
        src/p8est_iterate.h src/p8est_lnodes.h
libp4est_compiled_sources += \
        src/p8est_connectivity.c src/p8est.c \
        src/p8est_bits.c src/p8est_search.c \
        src/p8est_algorithms.c src/p8est_communication.c \
        src/p8est_ghost.c src/p8est_nodes.c src/p8est_vtk.c \
        src/p8est_trilinear.c src/p8est_trilinear_nodes.c \
        src/p8est_points.c src/p8est_geometry.c \
        src/p8est_iterate.c src/p8est_lnodes.c src/p8est_trilinear_lnodes.c
endif
P4EST_CPPFLAGS =

# build static library
lib_LIBRARIES += src/libp4est.a
src_libp4est_a_SOURCES = \
        $(libp4est_installed_headers) \
        $(libp4est_internal_headers) \
        $(libp4est_compiled_sources)
src_libp4est_a_CPPFLAGS = $(AM_CPPFLAGS) $(P4EST_CPPFLAGS)
src_libp4est_a_DEPENDENCIES = @P4EST_DEPENDENCIES@

# build shared library
if P4EST_ENABLE_SHARED
lib_LTLIBRARIES += src/libp4est.la
src_libp4est_la_SOURCES = \
        $(libp4est_internal_headers) \
        $(libp4est_compiled_sources)
src_libp4est_la_CPPFLAGS = $(AM_CPPFLAGS) $(P4EST_CPPFLAGS)
src_libp4est_la_DEPENDENCIES = @P4EST_DEPENDENCIES@
LDADD += @top_builddir@/src/libp4est.la @P4EST_SC_LDADD@
else
LDADD += @top_builddir@/src/libp4est.a @P4EST_SC_LDADD@
endif

nodist_include_HEADERS += $(libp4est_generated_headers)
include_HEADERS += $(libp4est_installed_headers)

if P4EST_SC_USE_MK
@P4EST_SC_MK_INCLUDE@
endif

AM_CPPFLAGS += -I@top_srcdir@/src @P4EST_SC_CPPFLAGS@

LINT_CSOURCES += $(libp4est_compiled_sources)
