/*
  This file is part of p4est.
  p4est is a C library to manage a parallel collection of quadtrees and/or
  octrees.

  Copyright (C) 2007 Carsten Burstedde, Lucas Wilcox.

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <p4est_base.h>
#include <p4est_algorithms.h>
#include <p4est_vtk.h>

typedef struct
{
  int32_t             a;
}
user_data_t;

static void
init_fn (p4est_t * p4est, int32_t which_tree, p4est_quadrant_t * quadrant)
{
  user_data_t        *data = quadrant->p.user_data;

  data->a = which_tree;
}

static int
refine_fn (p4est_t * p4est, int32_t which_tree, p4est_quadrant_t * quadrant)
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

int
main (int argc, char **argv)
{
  int                 num_procs, rank;
  int                 mpiret;
  MPI_Comm            mpicomm;
  p4est_t            *p4est;
  p4est_connectivity_t *connectivity;
  int64_t             qglobal, qlocal, qbegin, qend;
  int                 i;

  mpiret = MPI_Init (&argc, &argv);
  P4EST_CHECK_MPI (mpiret);
  mpicomm = MPI_COMM_WORLD;
  mpiret = MPI_Comm_rank (mpicomm, &rank);
  P4EST_CHECK_MPI (mpiret);

  sc_init (rank, NULL, NULL, NULL, SC_LP_DEFAULT);

  /* create connectivity and forest structures */
  connectivity = p4est_connectivity_new_corner ();
  p4est = p4est_new (mpicomm, connectivity, sizeof (user_data_t), init_fn);

  num_procs = p4est->mpisize;

  /* refine and balance to make the number of elements interesting */
  p4est_refine (p4est, refine_fn, init_fn);

  /* Check the global number of elements */
  qlocal = p4est->local_num_quadrants;
  qglobal = qlocal;
#ifdef P4EST_MPI
  mpiret = MPI_Allreduce (&qlocal, &qglobal, 1, MPI_LONG_LONG, MPI_SUM,
                          p4est->mpicomm);
  P4EST_CHECK_MPI (mpiret);
#endif
  P4EST_CHECK_ABORT (qglobal == p4est->global_num_quadrants,
                     "wrong number of p4est->global_num_quadrants");

  /* Check the number of elements per proc */
  for (i = 0; i < num_procs; ++i) {
    if (i == rank) {
      qlocal = p4est->local_num_quadrants;
    }
    else {
      qlocal = 0;
    }

    qglobal = qlocal;
#ifdef P4EST_MPI
    mpiret = MPI_Bcast (&qglobal, 1, MPI_LONG_LONG, i, p4est->mpicomm);
    P4EST_CHECK_MPI (mpiret);
#endif
    qbegin = (i == 0) ? 0 : p4est->global_last_quad_index[i - 1];
    qend = p4est->global_last_quad_index[i];
    P4EST_CHECK_ABORT (qglobal == qend - qbegin + 1,
                       "wrong number in p4est->global_last_quad_index");
  }

  /* clean up and exit */
  p4est_destroy (p4est);
  p4est_connectivity_destroy (connectivity);
  sc_finalize ();

  mpiret = MPI_Finalize ();
  P4EST_CHECK_MPI (mpiret);

  return 0;
}

/* EOF test_comm.c */
