
# This file is part of the p4est library
# Makefile.am in toplevel directory

ACLOCAL_AMFLAGS = -I config -I @P4EST_SC_CONFIG@

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

# set test environment
TESTS_ENVIRONMENT = @P4EST_MPI_TESTS_ENVIRONMENT@

# recurse only into subpackages
SUBDIRS =
LINT_SUBDIRS =
if P4EST_PROVIDE_SC
SUBDIRS += sc
LINT_SUBDIRS += sc
endif
DIST_SUBDIRS = $(SUBDIRS)

# handle toplevel directory
EXTRA_DIST += bootstrap p4estindent config/git-version-gen bugs
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
	@for subdir in $(LINT_SUBDIRS) ; do \
		echo "Making $@ in $$subdir"; \
		(cd $$subdir && $(MAKE) $(AM_MAKEFLAGS) lint) ; \
	done
	for f in $(LINT_CSOURCES) ; do \
		$(LINT) $(ALL_LINT_FLAGS) $(top_srcdir)/$$f || \
		echo "Lint check failed for $$f" ; \
	done
endif LINT

ChangeLog:
	(GIT_DIR=@top_srcdir@/.git git log | @P4EST_SC_DIR@/config/git2cl \
		> .ChangeLog.tmp && \
		mv .ChangeLog.tmp ChangeLog; rm -f .ChangeLog.tmp) || \
	(touch ChangeLog; \
		echo 'Warning not in git dir; ChangeLog could be wrong.' >&2)

dist-hook: ChangeLog
	echo $(VERSION) > $(distdir)/.tarball-version
	(test "x$(VERSION)" = "x`@top_srcdir@/config/git-version-gen @top_srcdir@/.tarball-version`" || (echo "Stale version;" && echo "Please run:" && echo "     cd @top_srcdir@ && ./bootstrap" && echo "before make dist" && exit 1 ))
if !P4EST_WITH_SC
	@echo "-----------------------------------------------------"
	@echo "make dist does not work with --without-sc"
	@echo "-----------------------------------------------------"
	@exit 1
endif

clean-local:
	rm -f *vtu *.p4c
