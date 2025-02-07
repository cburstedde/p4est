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

#include <p4est.h>
#include <p4est_vtk.h>
#include <p4est_extended.h>

/* refinement and coarsen level initialization */
static const int    refine_level = 6;
static const int    coarsen_level = 4;
static const int    partition_ext = 1;

/* refinement function */
static int
refine_fn (p4est_t * p4est, p4est_topidx_t which_tree,
           p4est_quadrant_t * quadrant)
{
  double              x_center, y_center, radius;

  /* Calculate the center coordinates of the quadrant */
  x_center =
    (double) (quadrant->x +
              (P4EST_QUADRANT_LEN (quadrant->level) / 2)) / P4EST_ROOT_LEN;
  y_center =
    (double) (quadrant->y +
              (P4EST_QUADRANT_LEN (quadrant->level) / 2)) / P4EST_ROOT_LEN;

  /* Define the radius of the circle */
  radius = 0.3;
  /* If the refinement level  */
  if ((int) quadrant->level <= refine_level) {
    /* Check if the center of the quadrant lies within the circle of radius 0.1 */
    if ((x_center - 0.5) * (x_center - 0.5) +
        (y_center - 0.5) * (y_center - 0.5) < radius * radius) {
      /* The center is within the circle, so refine this quadrant */
      return 1;
    }
  }
  return 0;
}

/* coarsen function */
static int
coarsen_fn (p4est_t * p4est, p4est_topidx_t which_tree,
            p4est_quadrant_t * quadrant[])
{
  /* TO DO: check whole family of quadrants */
  /* If the current level is greater than coarsen level, do the coarserning */
  if ((int) quadrant[0]->level >= coarsen_level) {
    P4EST_GLOBAL_LDEBUGF ("coarsening level [%d][%d]\n",
                          coarsen_level, (int) quadrant[0]->level);
    return 1;
  }
  P4EST_GLOBAL_LDEBUGF ("coarsening level [%d][%d]\n",
                        coarsen_level, (int) quadrant[0]->level);
  return 0;
}

int
main (int argc, char **argv)
{
  p4est_connectivity_t *conn;
  p4est_t            *p4est;
  sc_MPI_Comm         mpicomm = sc_MPI_COMM_WORLD;
  int                 mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  /* Exercise 1 */
  conn = p4est_connectivity_new_unitsquare ();
  p4est = p4est_new (mpicomm, conn, 0, NULL, NULL);
  p4est_vtk_write_file (p4est, NULL, P4EST_STRING "_unitsquare_new");

  /* Exercise 2 */
  p4est_refine (p4est, 1, refine_fn, NULL);
  p4est_vtk_write_file (p4est, NULL, P4EST_STRING "_unitsquare_refine");

  p4est_partition (p4est, 0, NULL);
  p4est_vtk_write_file (p4est, NULL, P4EST_STRING "_unitsquare_partition");

  /* Exercise 3
     Here are two options */
  p4est_coarsen (p4est, 1, coarsen_fn, NULL);

  if (partition_ext) {
    p4est_partition_ext (p4est, 1, NULL);
    p4est_vtk_write_file (p4est, NULL, P4EST_STRING "_unitsquare_balance");
  }
  else {
    p4est_partition_ext (p4est, 1, NULL);
  }

#if 0
  p4est_balance (p4est, P4EST_CONNECT_FULL, NULL);
#endif
  p4est_vtk_write_file (p4est, NULL,
                        P4EST_STRING "_unitsquare_partition_ext");

  p4est_destroy (p4est);
  p4est_connectivity_destroy (conn);

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
