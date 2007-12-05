
## Makefile.am in src

SUBDIRS = . libb64

lib_LIBRARIES = libp4est.a
libp4est_a_SOURCES = \
	p4est_base.h p4est_base.c p4est_mpi_dummy.h \
	p4est_obstack.h p4est_obstack.c p4est_memory.h p4est_memory.c \
	p4est_connectivity.h p4est_connectivity.c \
	p4est_file.h p4est_file.c p4est.h p4est.c \
	p4est_algorithms.h p4est_algorithms.c \
	p4est_vtk.h p4est_vtk.c

include_HEADERS = \
	p4est_mpi_dummy.h p4est_obstack.h p4est_memory.h \
	p4est_connectivity.h p4est_file.h p4est_algorithms.h \
	p4est_vtk.h p4est.h

AM_CPPFLAGS = -I@top_srcdir@/src/libb64

install-data-hook:
	perl -pi -e 's/(#define\s+P4EST_CONFIG_INSTALLED)\s+\d/\1 1/' \
		$(DESTDIR)$(includedir)/p4est.h
if HAVE_MPI
	perl -pi -e 's/(#define\s+P4EST_CONFIG_MPI)\s+\d/\1 1/' \
		$(DESTDIR)$(includedir)/p4est.h
endif
