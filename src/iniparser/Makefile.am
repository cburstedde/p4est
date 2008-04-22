
# Makefile.am in src/iniparser

bin_PROGRAMS = p4est_iniexample
p4est_iniexample_SOURCES = iniexample.c

lib_LIBRARIES = libiniparser.a
libiniparser_a_SOURCES = dictionary.h dictionary.c iniparser.h iniparser.c

EXTRA_DIST = AUTHORS LICENSE twisted.ini

LDADD = @top_builddir@/src/iniparser/libiniparser.a
