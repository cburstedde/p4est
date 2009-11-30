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

#include <p4est_algorithms.h>

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

int
main (int argc, char **argv)
{
  int                 num_procs, rank;
  int                 mpiret;
  MPI_Comm            mpicomm;
  p4est_t            *p4est;
  p4est_connectivity_t *connectivity;
  p4est_gloidx_t      qglobal, qlocal, qbegin, qend, qsum;
  int                 i;

  mpiret = MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpicomm = MPI_COMM_WORLD;
  mpiret = MPI_Comm_rank (mpicomm, &rank);
  SC_CHECK_MPI (mpiret);

  sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT);

  /* create connectivity and forest structures */
  connectivity = p4est_connectivity_new_corner ();
  p4est = p4est_new (mpicomm, connectivity, 15,
                     sizeof (user_data_t), init_fn, NULL);

  num_procs = p4est->mpisize;

  /* refine and balance to make the number of elements interesting */
  p4est_refine (p4est, 1, refine_fn, init_fn);

  /* Check the global number of elements */
  qlocal = p4est->local_num_quadrants;
  qglobal = -13473829;
  mpiret = MPI_Allreduce (&qlocal, &qglobal, 1, P4EST_MPI_GLOIDX, MPI_SUM,
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
    mpiret = MPI_Bcast (&qglobal, 1, P4EST_MPI_GLOIDX, i, p4est->mpicomm);
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

  mpiret = MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
