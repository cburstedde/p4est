
# This file is part of p4est.
# Makefile.am in toplevel directory

ACLOCAL_AMFLAGS = -I config @P4EST_SC_AMFLAGS@
if P4EST_SC_MK_USE
@P4EST_SC_MK_INCLUDE@
endif

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
lib_LTLIBRARIES =
nodist_include_HEADERS =
noinst_HEADERS =
noinst_PROGRAMS =
sysconf_DATA =
dist_p4estdata_DATA =

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
        config/p4est_include.m4 config/p4est_metis.m4

# install p4est data files in the correct directory
p4estdatadir = $(datadir)/data

# setup test environment
if P4EST_MPIRUN
LOG_COMPILER = @P4EST_MPIRUN@
AM_LOG_FLAGS = @P4EST_MPI_TEST_FLAGS@
endif

# recurse only into subpackages
SUBDIRS = @P4EST_SC_SUBDIR@
DIST_SUBDIRS = $(SUBDIRS)

# handle toplevel directory
EXTRA_DIST += \
        bootstrap p4estindent build-aux/git-version-gen build-aux/git2cl bugs doc
DISTCLEANFILES += \
        _configs.sed src/p4est_config.h @P4EST_DISTCLEAN@
.PHONY: lint ChangeLog

# non-recursive build
include src/Makefile.am
include test/Makefile.am
include example/mesh/Makefile.am
include example/points/Makefile.am
include example/simple/Makefile.am
include example/steps/Makefile.am
include example/tetgen/Makefile.am
include example/timings/Makefile.am
include example/balance_seeds/Makefile.am

# lint static syntax checker
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

# revision control and ChangeLog
ChangeLog:
	(GIT_DIR=@top_srcdir@/.git git log > .ChangeLog.tmp && \
         cat .ChangeLog.tmp | @top_srcdir@/build-aux/git2cl > ChangeLog) ; \
        rm -f .ChangeLog.tmp

dist-hook:
	echo $(VERSION) > $(distdir)/.tarball-version
	test "x$(VERSION)" = "x`@top_srcdir@/build-aux/git-version-gen\
              @top_srcdir@/.tarball-version`" || \
        ((echo "Stale version;" ; echo "Please run:" ; \
          echo "     cd @top_srcdir@ && ./bootstrap" ; \
          echo "before make dist") 1>&2 ; rm -r $(distdir) ; exit 1)
if P4EST_DIST_DENY
	@echo "-----------------------------------------------------"
	@echo "make dist does not work with external libsc"
	@echo "-----------------------------------------------------"
	rm -rf $(distdir)
	@exit 1
endif

clean-local: clean-local-p6est
	rm -f ChangeLog *vtu *.visit *.p4c *.p4p *.p8c *.p8p

maintainer-clean-local:
	rm -rf doxygen

doxygen: Doxyfile
	doxygen

.PHONY: doxygen
