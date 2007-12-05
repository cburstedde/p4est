
# Makefile.am in src/libb64

bin_PROGRAMS = p4est_b64enc p4est_b64dec
p4est_b64enc_SOURCES = b64enc.c
p4est_b64dec_SOURCES = b64dec.c

lib_LIBRARIES = libb64.a
libb64_a_SOURCES = cencode.c cdecode.c libb64.h

EXTRA_DIST = LICENSE CHANGELOG

LDADD = @top_builddir@/src/libb64/libb64.a
