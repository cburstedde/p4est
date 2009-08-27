
# This file is part of the p4est library
# Makefile.am in toplevel directory

ACLOCAL_AMFLAGS = -I config @P4EST_SC_AMFLAGS@

# initialize empty variables
AM_CPPFLAGS =
CLEANFILES =
DISTCLEANFILES =
EXTRA_DIST =
LDADD =
LINT_CSOURCES =
TESTS =
bin_PROGRAMS =
check_PROGRAMS =
include_HEADERS =
lib_LIBRARIES =
nodist_include_HEADERS =
noinst_HEADERS =
noinst_PROGRAMS =
sysconf_DATA =

# use this if you want to link in p4est without autotools
sysconf_DATA += Makefile.p4est.mk
CLEANFILES += Makefile.p4est.mk
Makefile.p4est.mk : Makefile.p4est.pre
	cat $< | \
        sed -e 's,{\(\(.*prefix\|sysconfdir\)\)},{p4est_\1},g' \
            -e 's,^\(\(.*prefix\|sysconfdir\) *=\),p4est_\1,g' > $@

# install p4est m4 macros in the correct directory
p4estaclocaldir = $(datadir)/aclocal
dist_p4estaclocal_DATA = \
        config/p4est_include.m4

# set test environment
TESTS_ENVIRONMENT = @P4EST_MPI_TESTS_ENVIRONMENT@

# recurse only into subpackages
SUBDIRS = @P4EST_SC_SUBDIR@
DIST_SUBDIRS = $(SUBDIRS)

# handle toplevel directory
EXTRA_DIST += \
        bootstrap p4estindent config/git-version-gen config/git2cl bugs doc
DISTCLEANFILES += \
        _configs.sed src/p4est_config.h P4EST.override.pre P4EST.override
.PHONY: lint ChangeLog

# non-recursive build
include src/Makefile.am
include test/Makefile.am
include example/geometry/Makefile.am
include example/points/Makefile.am
## include example/read_forest/Makefile.am
include example/second/Makefile.am
include example/simple/Makefile.am
include example/timings/Makefile.am

ALL_LINT_FLAGS = $(LINT_FLAGS) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) \
                 $(MPI_INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_CFLAGS) \
                 $(src_libp4est_a_CPPFLAGS)

lint:
if LINT
	@for subdir in $(SUBDIRS) ; do \
		echo "Making $@ in $$subdir"; \
		(cd $$subdir && $(MAKE) $(AM_MAKEFLAGS) lint) ; \
	done
	for f in $(LINT_CSOURCES) ; do \
		$(LINT) $(ALL_LINT_FLAGS) $(top_srcdir)/$$f || \
		echo "Lint check failed for $$f" ; \
	done
else
	@echo "Static source code check disabled by configure"
endif

ChangeLog:
	(GIT_DIR=@top_srcdir@/.git git log > .ChangeLog.tmp && \
         cat .ChangeLog.tmp | @top_srcdir@/config/git2cl > ChangeLog) ; \
        rm -f .ChangeLog.tmp

dist-hook:
	echo $(VERSION) > $(distdir)/.tarball-version
	test "x$(VERSION)" = "x`@top_srcdir@/config/git-version-gen\
              @top_srcdir@/.tarball-version`" || \
        ((echo "Stale version;" ; echo "Please run:" ; \
          echo "     cd @top_srcdir@ && ./bootstrap" ; \
          echo "before make dist") 1>&2 ; exit 1)
if !P4EST_DIST_ALLOW
	@echo "-----------------------------------------------------"
	@echo "make dist does not work with external libsc"
	@echo "-----------------------------------------------------"
	@exit 1
endif

clean-local:
	rm -f *vtu *.p4c *.p4p *.p8c *.p8p
