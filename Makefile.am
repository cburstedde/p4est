
# This file is part of p4est.
# Makefile.am in toplevel directory

ACLOCAL_AMFLAGS = -I config @P4EST_SC_AMFLAGS@
if P4EST_SC_MK_USE
@P4EST_SC_MK_INCLUDE@
endif

# initialize empty variables
AM_CPPFLAGS =
BUILT_SOURCES =
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

# pkg-config configuration file
pkgconfig_DATA = config/p4est_autotools.pc

# use this if you want to link in p4est without autotools
sysconf_DATA += config/Makefile.p4est.mk
CLEANFILES += config/Makefile.p4est.mk
config/Makefile.p4est.mk : config/Makefile.p4est.pre
	cat $< | \
        sed -e 's,{\(\(.*prefix\|sysconfdir\)\)},{p4est_\1},g' \
            -e 's,^\(\(.*prefix\|sysconfdir\) *=\),p4est_\1,g' > $@

# install p4est m4 macros in the correct directory
p4estaclocaldir = $(datadir)/aclocal
dist_p4estaclocal_DATA = \
        config/p4est_include.m4 config/p4est_metis.m4 config/p4est_petsc.m4

# install p4est data files in the correct directory
p4estdatadir = $(datadir)/data

# setup test environment
LOG_COMPILER = @P4EST_MPIRUN@ @P4EST_MPI_TEST_FLAGS@ \
               @P4EST_VALGRIND@ @P4EST_VALGRIND_FLAGS@

# recurse only into subpackages
SUBDIRS = @P4EST_SC_SUBDIR@
DIST_SUBDIRS = $(SUBDIRS)

# handle toplevel directory
EXTRA_DIST += \
        bootstrap p4estindent build-aux/git-version-gen build-aux/git2cl doc \
        cmake CMakeLists.txt src/CMakeLists.txt \
        example/CMakeLists.txt test/CMakeLists.txt CMakePresets.json

DISTCLEANFILES += @P4EST_DISTCLEAN@ \
        _configs.sed config/p4est_config.h config/p4est_autotools.pc
.PHONY: ChangeLog

# non-recursive build
include src/Makefile.am
include test/Makefile.am
include example/mesh/Makefile.am
include example/particles/Makefile.am
include example/points/Makefile.am
include example/simple/Makefile.am
include example/spheres/Makefile.am
include example/steps/Makefile.am
include example/tetgen/Makefile.am
include example/timings/Makefile.am
include example/balance_seeds/Makefile.am

## This was only used for our lint code, which needs to be replaced.
##
## # lint static syntax checker
## ALL_LINT_FLAGS = $(LINT_FLAGS) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) \
##                  $(MPI_INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_CFLAGS) \
##                  $(src_libp4est_a_CPPFLAGS)
## lint:
## if LINT
## 	@for subdir in $(SUBDIRS) ; do \
## 		echo "Making $@ in $$subdir"; \
## 		(cd $$subdir && $(MAKE) $(AM_MAKEFLAGS) lint) ; \
## 	done
## 	for f in $(LINT_CSOURCES) ; do \
## 		$(LINT) $(ALL_LINT_FLAGS) $(top_srcdir)/$$f || \
## 		echo "Lint check failed for $$f" ; \
## 	done
## else
## 	@echo "Static source code check disabled by configure"
## endif
## .PHONY: lint

# revision control and ChangeLog
ChangeLog:
	(GIT_DIR=@top_srcdir@/.git git log > .ChangeLog.tmp && \
         cat .ChangeLog.tmp | @top_srcdir@/build-aux/git2cl > ChangeLog) ; \
        rm -f .ChangeLog.tmp

dist-hook:
	echo $(VERSION) > $(distdir)/.tarball-version
	(GITGEN_VERSION=`(cd @top_srcdir@ && ./build-aux/git-version-gen\
                .tarball-version)` ; \
	 test "x$(VERSION)" = "x$${GITGEN_VERSION}" || \
         ((echo "Stale version;"; echo $(VERSION); echo "$${GITGEN_VERSION}"; \
           echo "Please run:" ; \
           echo "     (cd @top_srcdir@ && ./bootstrap)" ; \
           echo "before make dist") 1>&2 ; rm -r $(distdir) ; exit 1))
if P4EST_DIST_DENY
	@echo "-----------------------------------------------------"
	@echo "make dist does not work with external libsc"
	@echo "-----------------------------------------------------"
	rm -rf $(distdir)
	@exit 1
endif

install-data-hook:
	cd $(DESTDIR)$(pkgconfigdir) && \
              mv -f p4est_autotools.pc "p4est-$(VERSION).pc"

uninstall-hook:
	cd $(DESTDIR)$(pkgconfigdir) && rm -f "p4est-$(VERSION).pc"

## From the GNU automake documentation "What gets Cleaned"
## -------------------------------------------------------

## If make built it, and it is commonly something that one would want to rebuild
## (for instance, a .o file), then mostlyclean should delete it.
## Otherwise, if make built it, then clean should delete it.
clean-local:
	rm -f ChangeLog *vtu *.visit *.p4c *.p4p *.p8c *.p8p *.p4d *.p8d

## If configure built it, then distclean should delete it.
distclean-local:

## If the maintainer built it (for instance, a .info file), then
## maintainer-clean should delete it. However maintainer-clean should not delete
## anything that needs to exist in order to run ‘./configure && make’.
maintainer-clean-local:
	rm -rf doc/{html,latex,man,*.doxygen.tags}

doxygen-local: doc/Doxyfile
	cd doc && doxygen Doxyfile

.PHONY: doxygen
