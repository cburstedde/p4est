/*
  This file is part of p4est.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2010 The University of Texas System
  Additional copyright (C) 2011 individual authors
  Written by Carsten Burstedde, Lucas C. Wilcox, and Tobin Isaac

  p4est is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  p4est is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with p4est; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/* Purpose of this file is to test a connectivity that throws an assertion. */

#ifndef P4_TO_P8
#include <p4est_connectivity.h>
#include <p4est_extended.h>
#include <p4est_vtk.h>
#else
#include <p8est_connectivity.h>
#include <p8est_extended.h>
#include <p8est_vtk.h>
#endif

#ifdef P4_TO_P8

static p4est_connectivity_t *
p8est_connectivity_new_3flat (void)
{
  const p4est_topidx_t num_vertices = 14;
  const p4est_topidx_t num_trees = 3;
  const p4est_topidx_t num_ett = 0;
  const p4est_topidx_t num_ctt = 0;
/* *INDENT-OFF* */
  const double        vertices[14 * 3] = {
    0., 0., 0.,
    1., 0., 0.,
    2., 0., 0.,
    0., 0., 1.,
    1., 0., 1.,
    2., 0., 1.,
    0., 1., 0.,
    0., 2., 0.,
    0., 1., 1.,
    0., 2., 1.,
    1., 1., 0.,
    2., 2., 0., 
    1., 1., 1.,
    2., 2., 1.
  };
  const p4est_topidx_t tree_to_vertex[3 * 8] = {
     0,  1,  6, 10,  3,  4,  8, 12,
     1,  2, 10, 11,  4,  5, 12, 13,
    10, 11,  6,  7, 12, 13,  8,  9
  };
  const p4est_topidx_t tree_to_tree[3 * 6] = {
    0, 1, 0, 2, 0, 0,
    0, 1, 1, 2, 1, 1,
    0, 2, 1, 2, 2, 2
  };
  const int8_t        tree_to_face[3 * 6] = {
    0, 0, 2, 6, 4, 5,
    1, 1, 2, 2, 4, 5,
    9, 1, 3, 3, 4, 5
  };
/* *INDENT-ON* */

  return p4est_connectivity_new_copy (num_vertices, num_trees, 0, 0,
                                      vertices, tree_to_vertex,
                                      tree_to_tree, tree_to_face,
                                      NULL, &num_ett, NULL, NULL,
                                      NULL, &num_ctt, NULL, NULL);
}

static void
test_conn_which (sc_MPI_Comm mpicomm, p4est_connectivity_t * conn,
                 const char *cname)
{
  int                 retval;
  int                 mpiret;
  int                 mpirank;
  char                fname[BUFSIZ];
  p4est_connectivity_t *conn2;
  p4est_t            *p4est;

  /* save and load connectivity to file */
  mpiret = sc_MPI_Comm_rank (mpicomm, &mpirank);
  SC_CHECK_MPI (mpiret);
  if (mpirank == 0) {
    snprintf (fname, BUFSIZ, "%s_test_%s.p8c", P4EST_STRING, cname);

    retval = p4est_connectivity_save (fname, conn);
    SC_CHECK_ABORTF (!retval, "Failed saving %s", cname);

    conn2 = p4est_connectivity_load (fname, NULL);
    SC_CHECK_ABORTF (conn2 != NULL, "Failed loading %s", cname);
    SC_CHECK_ABORTF (p4est_connectivity_is_equal (conn, conn2),
                     "Unequal %s", cname);

    p4est_connectivity_destroy (conn2);
  }

  /* create a forest for visualization */
  snprintf (fname, BUFSIZ, "%s_test_%s", P4EST_STRING, cname);
  p4est = p4est_new_ext (mpicomm, conn, 0, 0, 1, 0, NULL, NULL);
  p4est_vtk_write_file (p4est, NULL, fname);
  p4est_destroy (p4est);

  /* destroy connectivity */
  p4est_connectivity_destroy (conn);
}

static void
test_conn_3flat (sc_MPI_Comm mpicomm)
{
  test_conn_which (mpicomm, p8est_connectivity_new_3flat (), "3flat");
}

static void
test_conn_3edge (sc_MPI_Comm mpicomm)
{
  test_conn_which (mpicomm, p8est_connectivity_new_edge (), "3edge");
}

#endif /* P4_TO_P8 */

int
main (int argc, char **argv)
{
  int                 mpiret;
  sc_MPI_Comm         mpicomm;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  mpicomm = sc_MPI_COMM_WORLD;
  sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

#ifdef P4_TO_P8
  test_conn_3flat (mpicomm);
  test_conn_3edge (mpicomm);
#endif

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
