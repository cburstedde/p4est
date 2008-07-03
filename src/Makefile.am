
# This file is part of the p4est library
# Makefile.am in src
# included non-recursively from toplevel directory

libp4est_internal_headers =
libp4est_generated_headers = src/p4est_config.h
libp4est_installed_headers = src/p4est_base.h
libp4est_compiled_sources = src/p4est_base.c
if P4EST_BUILD_2D
libp4est_installed_headers += \
        src/p4est_connectivity.h src/p4est.h src/p4est_bits.h \
        src/p4est_algorithms.h src/p4est_communication.h \
        src/p4est_file.h src/p4est_vtk.h src/p4est_mesh.h
libp4est_compiled_sources += \
        src/p4est_connectivity.c src/p4est.c src/p4est_bits.c \
        src/p4est_algorithms.c src/p4est_communication.c \
        src/p4est_file.c src/p4est_vtk.c src/p4est_mesh.c
endif
if P4EST_BUILD_3D
libp4est_internal_headers += src/p4est_to_p8est.h
libp4est_installed_headers += \
        src/p8est_connectivity.h src/p8est.h src/p8est_bits.h \
        src/p8est_algorithms.h src/p8est_communication.h \
        src/p8est_vtk.h src/p8est_mesh.h src/p8est_trilinear.h
#        src/p8est_file.h
libp4est_compiled_sources += \
        src/p8est_connectivity.c src/p8est.c src/p8est_bits.c \
        src/p8est_algorithms.c src/p8est_communication.c \
        src/p8est_vtk.c src/p8est_mesh.c src/p8est_trilinear.c
#        src/p8est_file.c
endif

lib_LIBRARIES += src/libp4est.a
src_libp4est_a_SOURCES = \
        $(libp4est_internal_headers) \
        $(libp4est_installed_headers) \
        $(libp4est_compiled_sources)

nodist_include_HEADERS += $(libp4est_generated_headers)
include_HEADERS += $(libp4est_installed_headers)

AM_CPPFLAGS += -I@top_srcdir@/src @P4EST_SC_CPPFLAGS@
LDADD += @top_builddir@/src/libp4est.a @P4EST_SC_LDADD@

LINT_CSOURCES += $(libp4est_compiled_sources)
