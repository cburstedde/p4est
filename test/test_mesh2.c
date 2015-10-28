/*
  This file is part of p4est.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2010 The University of Texas System
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

#ifndef P4_TO_P8
#include <p4est_extended.h>
#include <p4est_mesh.h>
#else /* !P4_TO_P8 */
#include <p8est_extended.h>
#include <p8est_mesh.h>
#endif /* !P4_TO_P8 */

/* Function for testing p4est-mesh for a single tree scenario
 *
 * \param [in] p4est     The forest.
 * \param [in] conn      The connectivity structure
 * \param [in] periodic  Flag for checking if we have periodic boundaries
 * \returns 0 for success, -1 for failure
 */
int
test_mesh_one_tree (p4est_t * p4est, p4est_connectivity_t * conn,
                    int8_t periodic, sc_MPI_Comm mpicomm)
{
  /* ensure that we have null pointers at beginning and end of function */
  P4EST_ASSERT (p4est == NULL);
  P4EST_ASSERT (conn == NULL);

  /* create connectivity */
#ifndef P4_TO_P8
  conn =
    periodic == 1 ? p4est_connectivity_new_periodic () :
    p4est_connectivity_new_unitsquare ();
#else /* !P4_TO_P8 */
  conn =
    periodic == 1 ? p8est_connectivity_new_periodic () :
    p8est_connectivity_new_unitcube ();
#endif /* !P4_TO_P8 */

  /* setup p4est */
  int                 minLevel = 1;
  p4est = p4est_new_ext (mpicomm, conn, 0, minLevel, 0, 0, 0, 0);

  /* cleanup */
  p4est_destroy (p4est);
  p4est_connectivity_destroy (conn);

  conn = 0;
  p4est = 0;

  P4EST_ASSERT (p4est == NULL);
  P4EST_ASSERT (conn == NULL);

  return 0;
}

/* Function for testing p4est-mesh for multiple trees in a brick scenario
 *
 * \param [in] p4est     The forest.
 * \param [in] conn      The connectivity structure
 * \param [in] periodic  Flag for checking if we have periodic boundaries
 * \returns 0 for success, -1 for failure
 */
int
test_mesh_multiple_trees_brick (p4est_t * p4est, p4est_connectivity_t * conn,
                                int8_t periodic, sc_MPI_Comm mpicomm)
{
  /* ensure that we have null pointers at beginning and end of function */
  P4EST_ASSERT (p4est == NULL);
  P4EST_ASSERT (conn == NULL);

  /* create connectivity */
#ifndef P4_TO_P8
  conn = p4est_connectivity_new_brick (1, 1, periodic, periodic);
#else /* !P4_TO_P8 */
  conn = p8est_connectivity_new_brick (1, 1, 1, periodic, periodic, periodic);
#endif /* !P4_TO_P8 */

  /* setup p4est */
  int                 minLevel = 1;
  p4est = p4est_new_ext (mpicomm, conn, 0, minLevel, 0, 0, 0, 0);
  /* cleanup */
  p4est_destroy (p4est);
  p4est_connectivity_destroy (conn);

  conn = 0;
  p4est = 0;

  P4EST_ASSERT (p4est == NULL);
  P4EST_ASSERT (conn == NULL);

  return 0;
}

/* Function for testing p4est-mesh for multiple trees in a non-brick scenario
 *
 * \param [in] p4est     The forest.
 * \param [in] conn      The connectivity structure
 * \param [in] periodic  Flag for checking if we have periodic boundaries
 * \returns 0 for success, -1 for failure
 */
int
test_mesh_multiple_trees_nonbrick (p4est_t * p4est,
                                   p4est_connectivity_t * conn,
                                   int8_t periodic, sc_MPI_Comm mpicomm)
{
  return 0;
}

int
main (int argc, char **argv)
{
  sc_MPI_Comm         mpicomm;
  int                 mpiret;
  int                 mpisize, mpirank;
  p4est_t            *p4est;
  p4est_connectivity_t *conn;
  int8_t              periodic_boundaries;

  /* initialize MPI */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpicomm = sc_MPI_COMM_WORLD;
  mpiret = sc_MPI_Comm_size (mpicomm, &mpisize);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (mpicomm, &mpirank);
  SC_CHECK_MPI (mpiret);

  sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  p4est = 0;
  conn = 0;
  /* test both periodic and non-periodic boundaries */
  /* test one tree */
  periodic_boundaries = 0;
  test_mesh_one_tree (p4est, conn, periodic_boundaries, mpicomm);

  periodic_boundaries = 1;
  test_mesh_one_tree (p4est, conn, periodic_boundaries, mpicomm);

  /* test multiple trees; brick */
  periodic_boundaries = 0;
  test_mesh_multiple_trees_brick (p4est, conn, periodic_boundaries, mpicomm);

  periodic_boundaries = 1;
  test_mesh_multiple_trees_brick (p4est, conn, periodic_boundaries, mpicomm);

  /* test multiple trees; non-brick */
  periodic_boundaries = 0;
  test_mesh_multiple_trees_nonbrick (p4est, conn, periodic_boundaries,
                                     mpicomm);

  periodic_boundaries = 1;
  test_mesh_multiple_trees_nonbrick (p4est, conn, periodic_boundaries,
                                     mpicomm);

  /* exit */
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
