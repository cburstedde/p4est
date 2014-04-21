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

#include <p4est_algorithms.h>
#include <p4est_communication.h>
#include <p4est_extended.h>

typedef struct
{
  p4est_topidx_t      a;
}
user_data_t;

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
  if (quadrant->x == P4EST_LAST_OFFSET (2) &&
      quadrant->y == P4EST_LAST_OFFSET (2)) {
    return 1;
  }
  if (quadrant->x >= P4EST_QUADRANT_LEN (2)) {
    return 0;
  }

  return 1;
}

static void
test_pertree (p4est_t * p4est)
{
  p4est_topidx_t      num_trees;
  p4est_gloidx_t     *pertree;

  num_trees = p4est->connectivity->num_trees;
  P4EST_ASSERT ((size_t) num_trees == p4est->trees->elem_count);
  pertree = P4EST_ALLOC (p4est_gloidx_t, num_trees + 1);
  p4est_comm_count_pertree (p4est, pertree);
  SC_CHECK_ABORT (pertree[num_trees] == p4est->global_num_quadrants,
                  "pertree check failed");
  P4EST_FREE (pertree);
}

int
main (int argc, char **argv)
{
  int                 num_procs, rank;
  int                 mpiret;
  sc_MPI_Comm         mpicomm;
  p4est_t            *p4est;
  p4est_connectivity_t *connectivity;
  p4est_gloidx_t      qglobal, qlocal, qbegin, qend, qsum;
  int                 i;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpicomm = sc_MPI_COMM_WORLD;
  mpiret = sc_MPI_Comm_rank (mpicomm, &rank);
  SC_CHECK_MPI (mpiret);

  sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT);

  /* create connectivity and forest structures */
  connectivity = p4est_connectivity_new_corner ();
  p4est = p4est_new_ext (mpicomm, connectivity, 15, 0, 0,
                         sizeof (user_data_t), init_fn, NULL);

  num_procs = p4est->mpisize;

  /* test tree counting */
  test_pertree (p4est);

  /* refine and balance to make the number of elements interesting */
  p4est_refine (p4est, 1, refine_fn, init_fn);

  /* test tree counting */
  test_pertree (p4est);

  /* Check the global number of elements */
  qlocal = p4est->local_num_quadrants;
  qglobal = -13473829;
  mpiret =
    sc_MPI_Allreduce (&qlocal, &qglobal, 1, P4EST_MPI_GLOIDX, sc_MPI_SUM,
                      p4est->mpicomm);
  SC_CHECK_MPI (mpiret);
  SC_CHECK_ABORT (qglobal == p4est->global_num_quadrants,
                  "wrong number of p4est->global_num_quadrants");

  /* Check the number of elements per proc */
  qsum = 0;
  for (i = 0; i < num_procs; ++i) {
    if (i == rank) {
      qlocal = p4est->local_num_quadrants;
    }
    else {
      qlocal = 0;
    }

    qglobal = qlocal;
    mpiret = sc_MPI_Bcast (&qglobal, 1, P4EST_MPI_GLOIDX, i, p4est->mpicomm);
    SC_CHECK_MPI (mpiret);

    qsum += qglobal;
    qbegin = p4est->global_first_quadrant[i];
    qend = p4est->global_first_quadrant[i + 1];
    SC_CHECK_ABORT (qglobal == qend - qbegin,
                    "wrong number in p4est->global_first_quadrant");
  }
  SC_CHECK_ABORT (qsum == p4est->global_num_quadrants,
                  "Wrong number after quadrant counting");

  /* clean up and exit */
  p4est_destroy (p4est);
  p4est_connectivity_destroy (connectivity);
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
