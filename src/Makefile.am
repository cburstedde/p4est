
# This file is part of p4est.
# Makefile.am in src
# included non-recursively from toplevel directory

libp4est_generated_headers = config/p4est_config.h
libp4est_installed_headers = src/p4est_base.h
libp4est_internal_headers =
libp4est_compiled_sources = src/p4est_base.c
if P4EST_ENABLE_BUILD_2D
libp4est_installed_headers += \
        src/p4est_connectivity.h src/p4est.h src/p4est_extended.h \
        src/p4est_bits.h src/p4est_search.h src/p4est_build.h \
        src/p4est_algorithms.h src/p4est_communication.h \
        src/p4est_ghost.h src/p4est_nodes.h src/p4est_vtk.h \
        src/p4est_points.h src/p4est_geometry.h \
        src/p4est_iterate.h src/p4est_lnodes.h src/p4est_mesh.h \
        src/p4est_balance.h src/p4est_io.h \
        src/p4est_wrap.h src/p4est_plex.h \
        src/p4est_empty.h
libp4est_compiled_sources += \
        src/p4est_connectivity.c src/p4est.c \
        src/p4est_bits.c src/p4est_search.c src/p4est_build.c \
        src/p4est_algorithms.c src/p4est_communication.c \
        src/p4est_ghost.c src/p4est_nodes.c src/p4est_vtk.c \
        src/p4est_points.c src/p4est_geometry.c \
        src/p4est_iterate.c src/p4est_lnodes.c src/p4est_mesh.c \
        src/p4est_balance.c src/p4est_io.c \
        src/p4est_connrefine.c \
        src/p4est_wrap.c src/p4est_plex.c \
        src/p4est_empty.c
endif
if P4EST_ENABLE_BUILD_3D
libp4est_installed_headers += \
        src/p4est_to_p8est.h \
        src/p8est_connectivity.h src/p8est.h src/p8est_extended.h \
        src/p8est_bits.h src/p8est_search.h src/p8est_build.h \
        src/p8est_algorithms.h src/p8est_communication.h \
        src/p8est_ghost.h src/p8est_nodes.h src/p8est_vtk.h \
        src/p8est_points.h src/p8est_geometry.h \
        src/p8est_iterate.h src/p8est_lnodes.h src/p8est_mesh.h \
        src/p8est_tets_hexes.h src/p8est_balance.h src/p8est_io.h \
        src/p8est_wrap.h src/p8est_plex.h \
        src/p8est_empty.h src/p4est_to_p8est_empty.h
libp4est_compiled_sources += \
        src/p8est_connectivity.c src/p8est.c \
        src/p8est_bits.c src/p8est_search.c src/p8est_build.c  \
        src/p8est_algorithms.c src/p8est_communication.c \
        src/p8est_ghost.c src/p8est_nodes.c src/p8est_vtk.c \
        src/p8est_points.c src/p8est_geometry.c \
        src/p8est_iterate.c src/p8est_lnodes.c src/p8est_mesh.c \
        src/p8est_tets_hexes.c src/p8est_balance.c src/p8est_io.c \
        src/p8est_connrefine.c \
        src/p8est_wrap.c src/p8est_plex.c \
        src/p8est_empty.c
endif
if P4EST_ENABLE_BUILD_2D
if P4EST_ENABLE_BUILD_3D
if P4EST_ENABLE_BUILD_P6EST
libp4est_installed_headers += \
        src/p6est.h src/p6est_ghost.h src/p6est_lnodes.h \
        src/p6est_profile.h src/p6est_vtk.h \
        src/p6est_extended.h src/p6est_communication.h \
        src/p6est_empty.h
libp4est_compiled_sources += \
        src/p6est.c src/p6est_ghost.c src/p6est_lnodes.c \
        src/p6est_profile.c src/p6est_vtk.c \
        src/p6est_communication.c \
        src/p6est_empty.c
endif
endif
endif

# this variable is used for headers that are not publicly installed
P4EST_CPPFLAGS =

# read the .so version from configuration file
include config/p4est_soversion.in

justlibs-local: src/libp4est.la
lib_LTLIBRARIES += src/libp4est.la
src_libp4est_la_SOURCES = \
        $(libp4est_internal_headers) \
        $(libp4est_compiled_sources)
src_libp4est_la_CPPFLAGS = $(AM_CPPFLAGS) $(P4EST_CPPFLAGS)
## This is the official API versioning scheme of libtool.  Please see:
## Read https://www.gnu.org/software/libtool/manual/libtool.html#Versioning
src_libp4est_la_LDFLAGS = -version-info $(P4EST_SOVERSION)
src_libp4est_la_LIBADD = @P4EST_SC_LIBADD@
EXTRA_src_libp4est_la_DEPENDENCIES = @P4EST_SC_EDEPS@

## Code libsc library installation path into p4est executables
AM_LDFLAGS += @P4EST_SC_RPATH@
LDADD += src/libp4est.la @P4EST_SC_LDADD@

nodist_include_HEADERS += $(libp4est_generated_headers)
include_HEADERS += $(libp4est_installed_headers)

AM_CPPFLAGS += -I@top_srcdir@/src @P4EST_SC_CPPFLAGS@
