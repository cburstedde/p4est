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

#ifdef P4_TO_P8
#include <p8est_algorithms.h>
#include <p8est_bits.h>
#include <p8est_extended.h>
#include <p8est_vtk.h>
#else
#include <p4est_algorithms.h>
#include <p4est_bits.h>
#include <p4est_extended.h>
#include <p4est_vtk.h>
#endif

/* typedefs */
typedef struct
{
  p4est_topidx_t      a;
}
user_data_t;

/* global variables */
static int          coarsen_all = 1;

/* functions */
static void
init_fn (p4est_t * p4est, p4est_topidx_t which_tree,
         p4est_quadrant_t * quadrant)
{
  user_data_t        *data = (user_data_t *) quadrant->p.user_data;

  data->a = which_tree;
}

static int
refine_fn (p4est_t * p4est, p4est_topidx_t which_tree,
           p4est_quadrant_t * quadrant)
{
  if (quadrant->level >= 6) {
    return 0;
  }
#ifdef P4_TO_P8
  if (quadrant->level >= 5 && quadrant->z <= P4EST_QUADRANT_LEN (3)) {
    return 0;
  }
#endif

  if (quadrant->x == P4EST_LAST_OFFSET (2) &&
      quadrant->y == P4EST_LAST_OFFSET (2)) {
    return 1;
  }
  if (quadrant->x >= P4EST_QUADRANT_LEN (2)) {
    return 0;
  }

  return 1;
}

static int
coarsen_fn (p4est_t * p4est, p4est_topidx_t which_tree,
            p4est_quadrant_t * q[])
{
  SC_CHECK_ABORT (p4est_quadrant_is_familypv (q), "Coarsen invocation");

  return coarsen_all || q[0]->y >= P4EST_ROOT_LEN / 2;
}

/* main */
int
main (int argc, char **argv)
{
  int                 rank, num_procs, mpiret, i;
  sc_MPI_Comm         mpicomm = sc_MPI_COMM_WORLD;
  p4est_t            *p4est_1tree, *p4est_ntrees;
  p4est_connectivity_t *connectivity_1tree, *connectivity_ntrees;

  /* initialize MPI and p4est internals */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_size (mpicomm, &num_procs);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (mpicomm, &rank);
  SC_CHECK_MPI (mpiret);

  sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  /* create connectivity */
#ifdef P4_TO_P8
  connectivity_1tree = p8est_connectivity_new_unitcube ();
  connectivity_ntrees = p8est_connectivity_new_twocubes ();
#else
  connectivity_1tree = p4est_connectivity_new_unitsquare ();
  connectivity_ntrees = p4est_connectivity_new_corner ();
#endif

  /* create p4est structure */
  p4est_1tree = p4est_new_ext (mpicomm, connectivity_1tree, 15, 0, 0,
                               sizeof (user_data_t), init_fn, NULL);

  p4est_ntrees = p4est_new_ext (mpicomm, connectivity_ntrees, 15, 0, 0,
                                sizeof (user_data_t), init_fn, NULL);

  /* write output: new */
  p4est_vtk_write_file (p4est_1tree, NULL,
                        P4EST_STRING "_partition_corr_1tree_new");
  p4est_vtk_write_file (p4est_ntrees, NULL,
                        P4EST_STRING "_partition_corr_ntrees_new");

  /* refine */
  p4est_refine (p4est_1tree, 1, refine_fn, init_fn);
  p4est_refine (p4est_ntrees, 1, refine_fn, init_fn);

  /* write output: refined */
  p4est_vtk_write_file (p4est_1tree, NULL,
                        P4EST_STRING "_partition_corr_1tree_refined");
  p4est_vtk_write_file (p4est_ntrees, NULL,
                        P4EST_STRING "_partition_corr_ntrees_refined");

  /* run partition and coarsen till one quadrant per tree remains */
  i = 0;
  while (p4est_1tree->global_num_quadrants > 1 && i <= P4EST_MAXLEVEL) {
    (void) p4est_partition_ext (p4est_1tree, 1, NULL);
    p4est_coarsen (p4est_1tree, 0, coarsen_fn, init_fn);
    i++;
  }
  SC_CHECK_ABORT (p4est_1tree->global_num_quadrants == 1,
                  "coarsest forest with one tree was not achieved");

  i = 0;
  while (p4est_ntrees->global_num_quadrants > connectivity_ntrees->num_trees
         && i <= P4EST_MAXLEVEL) {
    (void) p4est_partition_ext (p4est_ntrees, 1, NULL);
    p4est_coarsen (p4est_ntrees, 0, coarsen_fn, init_fn);
    i++;
  }
  SC_CHECK_ABORT (p4est_ntrees->global_num_quadrants
                  == connectivity_ntrees->num_trees,
                  "coarsest forest with multiple trees was not achieved");

  /* run partition on coarse forest (one quadrant per tree) once again */
  (void) p4est_partition_ext (p4est_1tree, 1, NULL);
  (void) p4est_partition_ext (p4est_ntrees, 1, NULL);

  /* write output: coarsened */
  p4est_vtk_write_file (p4est_1tree, NULL,
                        P4EST_STRING "_partition_corr_1tree_coarsened");
  p4est_vtk_write_file (p4est_ntrees, NULL,
                        P4EST_STRING "_partition_corr_ntrees_coarsened");

  /* destroy the p4est and its connectivity structure */
  p4est_destroy (p4est_1tree);
  p4est_destroy (p4est_ntrees);
  p4est_connectivity_destroy (connectivity_1tree);
  p4est_connectivity_destroy (connectivity_ntrees);

  /* clean up and exit */
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
