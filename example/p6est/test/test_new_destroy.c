/*
  This file is part of p4est.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2014 The University of Texas System
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

#include <p4est_bits.h>
#include <p6est.h>
#include <p6est_ghost.h>
#include <p6est_vtk.h>

char               *TEST_USER_POINTER;

static int          refine_level = -1;

/* To define a p6est_refine_column_t, all we have to do is take a p4est refine
 * function ... */
static int
p4est_refine_fn (p4est_t * p4est, p4est_topidx_t which_tree,
                 p4est_quadrant_t * quadrant)
{
  int                 cid;

  if (quadrant->level >= refine_level) {
    return 0;
  }
  if (which_tree == 2 || which_tree == 3) {
    return 0;
  }

  cid = p4est_quadrant_child_id (quadrant);

  if (cid == P4EST_CHILDREN - 1 ||
      (quadrant->x >= P4EST_LAST_OFFSET (P4EST_MAXLEVEL - 2) &&
       quadrant->y >= P4EST_LAST_OFFSET (P4EST_MAXLEVEL - 2)
#ifdef P4_TO_P8
       && quadrant->z >= P4EST_LAST_OFFSET (P4EST_MAXLEVEL - 2)
#endif
      )) {
    return 1;
  }
  if (quadrant->level == 1 && cid == 2) {
    return 1;
  }
  if (quadrant->x == P4EST_QUADRANT_LEN (2) &&
      quadrant->y == P4EST_LAST_OFFSET (2)) {
    return 1;
  }
  if (quadrant->y >= P4EST_QUADRANT_LEN (2)) {
    return 0;
  }

  return 1;
}

/* and wrap it.*/
static int
refine_column_fn (p6est_t * p6est, p4est_topidx_t which_tree,
                  p4est_quadrant_t * column)
{
  return p4est_refine_fn (p6est->columns, which_tree, column);
}

static int
refine_layer_fn (p6est_t * p6est, p4est_topidx_t which_tree,
                 p4est_quadrant_t * column, p2est_quadrant_t * layer)
{
  p4est_topidx_t      tohash[4];
  unsigned            hash;

  tohash[0] = (p4est_topidx_t) column->x;
  tohash[1] = (p4est_topidx_t) column->y;
  tohash[2] = (p4est_topidx_t) layer->z;
  tohash[3] = (((p4est_topidx_t) column->level) << 16) |
    ((p4est_topidx_t) layer->level);

  hash = p4est_topidx_hash4 (tohash);

  return (layer->level < refine_level && !((int) hash % 3));
}

void
init_fn (p6est_t * p6est, p4est_topidx_t which_tree,
         p4est_quadrant_t * col, p2est_quadrant_t * layer)
{
  SC_CHECK_ABORT (p6est->user_pointer == TEST_USER_POINTER,
                  "user_pointer corruption\n");
}

static int
coarsen_column_fn (p6est_t * p6est, p4est_topidx_t which_tree,
                   p4est_quadrant_t * column[])
{
  return 1;
}

static int
weight_fn (p6est_t * p6est, p4est_topidx_t which_tree,
           p4est_quadrant_t * col, p2est_quadrant_t * layer)
{
  return 1;
}

static int
coarsen_layer_fn (p6est_t * p6est, p4est_topidx_t which_tree,
                  p4est_quadrant_t * column, p2est_quadrant_t * layers[])
{
  return 1;
}

int
main (int argc, char **argv)
{
  MPI_Comm            mpicomm = MPI_COMM_WORLD;
  p4est_connectivity_t *conn4;
  p6est_connectivity_t *conn;
  p6est_t            *p6est, *copy_p6est;
  p6est_ghost_t      *ghost;
  double              height[3] = { 0., 0., 0.1 };
  int                 mpiret;

  mpiret = MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  SC_CHECK_ABORTF (argc == 2,
                   "Usage:\n%s NAME\n"
                   "  NAME=<corner|cubed|disk|periodic|rotwrap|star|unit>\n",
                   argv[0]);

  conn4 = p4est_connectivity_new_byname (argv[1]);
  conn = p6est_connectivity_new (conn4, NULL, height);

  p4est_connectivity_destroy (conn4);

  p6est = p6est_new (mpicomm, conn, 4, init_fn, TEST_USER_POINTER);
  p6est_destroy (p6est);

  p6est = p6est_new_ext (mpicomm, conn, 0, 1, 2, 1, 3, init_fn,
                         TEST_USER_POINTER);

  refine_level = 3;
  p6est_refine_columns (p6est, 1, refine_column_fn, init_fn);
  refine_level = 4;
  p6est_refine_layers (p6est, 1, refine_layer_fn, init_fn);
  refine_level = 5;
  p6est_refine_columns (p6est, 1, refine_column_fn, init_fn);

  copy_p6est = p6est_copy (p6est, 1);
  p6est_coarsen_columns (copy_p6est, 1, coarsen_column_fn, init_fn);
  p6est_vtk_write_file (copy_p6est, "p6est_test_coarsen_columns");
  p6est_destroy (copy_p6est);

  copy_p6est = p6est_copy (p6est, 1);
  p6est_coarsen_layers (copy_p6est, 0, coarsen_layer_fn, init_fn);
  p6est_vtk_write_file (copy_p6est, "p6est_test_coarsen_layers");
  p6est_destroy (copy_p6est);

  p6est_partition (p6est, weight_fn);
  p6est_partition (p6est, NULL);

  p6est_vtk_write_file (p6est, "p6est_test_new_destroy");

  ghost = p6est_ghost_new (p6est, P4EST_CONNECT_FACE);
  p6est_ghost_destroy (ghost);
  ghost = p6est_ghost_new (p6est, P4EST_CONNECT_FULL);
  p6est_ghost_destroy (ghost);

  p6est_destroy (p6est);

  p6est_connectivity_destroy (conn);

  /* exit */
  sc_finalize ();

  mpiret = MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
