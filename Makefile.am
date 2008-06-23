
# This file is part of the p4est library
# Makefile.am in toplevel directory

ACLOCAL_AMFLAGS = -I config -I sc/config

# initialize empty variables
AM_CPPFLAGS =
DISTCLEANFILES =
EXTRA_DIST =
LDADD =
LINT_CSOURCES =
MAINTAINERCLEANFILES =
TESTS =
bin_PROGRAMS =
check_PROGRAMS =
include_HEADERS =
lib_LIBRARIES =
nodist_include_HEADERS =
noinst_HEADERS =
noinst_PROGRAMS =

# recurse only into subpackages
SUBDIRS =
if P4EST_PROVIDE_SC
SUBDIRS += sc
endif
DIST_SUBDIRS = $(SUBDIRS)

# handle toplevel directory
EXTRA_DIST += bootstrap p4estindent
DISTCLEANFILES += _configs.sed src/p4est_config.h
MAINTAINERCLEANFILES += ChangeLog
.PHONY: lint ChangeLog

# non-recursive build
include src/Makefile.am
include test/Makefile.am
include example/read_forest/Makefile.am
include example/second/Makefile.am
include example/simple/Makefile.am
include example/timings/Makefile.am
include doc/Makefile.am

ALL_LINT_FLAGS = $(LINT_FLAGS) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) \
                 $(MPI_INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_CFLAGS) \
                 $(src_libp4est_a_CPPFLAGS)

lint:
if LINT
	for f in $(LINT_CSOURCES) ; do \
		$(LINT) $(ALL_LINT_FLAGS) $(top_srcdir)/$$f || \
		echo "Lint check failed for $$f" ; \
	done
endif LINT

ChangeLog:
	@(GIT_DIR=@top_srcdir@/.git git-log | @top_srcdir@/sc/git2cl \
		> .ChangeLog.tmp && \
		mv .ChangeLog.tmp ChangeLog; rm -f .ChangeLog.tmp) || \
	(touch ChangeLog; \
		echo 'Warning not in git dir; ChangeLog could be wrong.' >&2)

dist-hook: ChangeLog
if !P4EST_PROVIDE_SC
	@echo "-----------------------------------------------------"
	@echo "make dist requires the builtin libsc to be configured"
	@echo "-----------------------------------------------------"
	@exit 1
endif

clean-local:
	rm -f *vtu
