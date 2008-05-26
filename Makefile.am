
# This file is part of the p4est Library
# Makefile.am in toplevel directory

ACLOCAL_AMFLAGS = -I config -I sc/config

SUBDIRS =
if P4EST_PROVIDE_SC
SUBDIRS += sc
endif
SUBDIRS += src test example doc

EXTRA_DIST = bootstrap p4estindent

MAINTAINERCLEANFILES = ChangeLog

.PHONY: lint ChangeLog

ALL_LINT_FLAGS = $(LINT_FLAGS) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) \
                 $(MPI_INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_CFLAGS) \
                 $(src_libsc_a_CPPFLAGS) -DP4EST_SPLINT

lint:
if LINT
	for f in $(LINT_CSOURCES) ; do \
		$(LINT) $(ALL_LINT_FLAGS) $(top_srcdir)/$$f || \
		echo "Lint check failed for $$f" ; \
	done
endif LINT

ChangeLog:
	(GIT_DIR=@top_srcdir@/.git git-log | @top_srcdir@/git2cl \
		> .ChangeLog.tmp && \
		mv .ChangeLog.tmp ChangeLog; rm -f .ChangeLog.tmp) || \
	(touch ChangeLog; \
		echo 'Warning not in git dir; ChangeLog could be wrong.' >&2)

dist-hook: ChangeLog
