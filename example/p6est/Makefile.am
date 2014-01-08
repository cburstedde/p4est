
if P4EST_ENABLE_BUILD_P6EST
libp4est_installed_headers += \
        example/p6est/p6est.h example/p6est/p6est_ghost.h example/p6est/p6est_lnodes.h \
        example/p6est/p6est_vtk.h example/p6est/p6est_io.h
libp4est_compiled_sources += \
        example/p6est/p6est.c example/p6est/p6est_ghost.c example/p6est/p6est_lnodes.c \
        example/p6est/p6est_vtk.c example/p6est/p6est_io.c
AM_CPPFLAGS += -I@top_srcdir@/example/p6est

p6est_test_programs = \
                      example/p6est/test/p6est_test_new_destroy

example_p6est_test_p6est_test_new_destroy_SOURCES = example/p6est/test/test_new_destroy.c
TESTS += $(p6est_test_programs)
check_PROGRAMS += $(p6est_test_programs)
endif
