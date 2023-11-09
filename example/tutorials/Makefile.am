
# This file is part of p4est.
# Makefile.am in example/tutorials/
# included non-recursively from toplevel directory

bin_PROGRAMS += example/tutorials/p4est_tutorials_t0_new
bin_PROGRAMS += example/tutorials/p4est_tutorials_t1_conn_tree
bin_PROGRAMS += example/tutorials/p4est_tutorials_t1_conn_fileio
bin_PROGRAMS += example/tutorials/p4est_tutorials_t2_forest

example_tutorials_p4est_tutorials_t0_new_SOURCES = example/tutorials/t0_new.c
example_tutorials_p4est_tutorials_t1_conn_tree_SOURCES = example/tutorials/t1_conn_tree.c
example_tutorials_p4est_tutorials_t1_conn_fileio_SOURCES = example/tutorials/t1_conn_fileio.c
example_tutorials_p4est_tutorials_t2_forest_SOURCES = example/tutorials/t2_forest.c

LINT_CSOURCES += $(example_tutorials_p4est_tutorials_t0_new_SOURCES)
LINT_CSOURCES += $(example_tutorials_p4est_tutorials_t1_conn_tree_SOURCES)
LINT_CSOURCES += $(example_tutorials_p4est_tutorials_t2_forest_SOURCES)
