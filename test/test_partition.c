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
#include <p4est_base.h>
#include <p4est_vtk.h>

typedef struct
{
  int32_t             a;
  int64_t             sum;
}
user_data_t;

static void
init_fn (p4est_t * p4est, int32_t which_tree, p4est_quadrant_t * quadrant)
{
  user_data_t        *data = quadrant->user_data;

  data->a = which_tree;
  data->sum = quadrant->x + quadrant->y + quadrant->level;
}

static int
refine_fn (p4est_t * p4est, int32_t which_tree, p4est_quadrant_t * quadrant)
{
  if (quadrant->level >= 6) {
    return 0;
  }
  if (quadrant->x == (1 << (P4EST_MAXLEVEL)) - (1 << (P4EST_MAXLEVEL - 2)) &&
      quadrant->y == (1 << (P4EST_MAXLEVEL)) - (1 << (P4EST_MAXLEVEL - 2))) {
    return 1;
  }
  if (quadrant->x >= (1 << (P4EST_MAXLEVEL - 2))) {
    return 0;
  }

  return 1;
}

int
main (int argc, char **argv)
{
#ifdef HAVE_MPI
  int                 mpiret;
#endif
  MPI_Comm            mpicomm;
  p4est_t            *p4est;
  p4est_connectivity_t *connectivity;
  int32_t             num_procs, rank;
  int32_t             i;
  int32_t             num_quadrants_on_last;
  int32_t            *num_quadrants_in_proc;
  int32_t            *num_quadrants_in_proc_check;
  int                 t, q;
  p4est_quadrant_t   *quad;
  p4est_tree_t       *tree;
  user_data_t        *user_data;
  int64_t             sum;
  unsigned            crc;

  mpicomm = MPI_COMM_NULL;
#ifdef HAVE_MPI
  mpiret = MPI_Init (&argc, &argv);
  P4EST_CHECK_MPI (mpiret);
  mpicomm = MPI_COMM_WORLD;
#endif

  /* create connectivity and forest structures */
  connectivity = p4est_connectivity_new_corner ();
  p4est = p4est_new (mpicomm, NULL, connectivity,
                     sizeof (user_data_t), init_fn);

  num_procs = p4est->mpisize;
  rank = p4est->mpirank;

  num_quadrants_in_proc = P4EST_ALLOC (int32_t, num_procs);
  P4EST_CHECK_ALLOC (num_quadrants_in_proc);

  /* refine and balance to make the number of elements interesting */
  p4est_refine (p4est, refine_fn, init_fn);

  /* Set an arbitrary partition.
   *
   * Since this is just a test we assume the global number of
   * quadrants will fit in an int32_t
   */
  num_quadrants_on_last = (int32_t) p4est->global_num_quadrants;
  for (i = 0; i < num_procs - 1; ++i) {
    num_quadrants_in_proc[i] = i + 1;
    num_quadrants_on_last -= i + 1;
  }
  num_quadrants_in_proc[num_procs - 1] = num_quadrants_on_last;
  P4EST_CHECK_ABORT (num_quadrants_on_last > 0,
                     "Negative number of quadrants on the last processor");

  /* Save a checksum of the original forest */
  crc = p4est_checksum (p4est);

  /* partition the forest */
  p4est_partition_given (p4est, num_quadrants_in_proc);

  /* Double check that we didn't loose any quads */
  P4EST_CHECK_ABORT (crc == p4est_checksum (p4est),
                     "bad checksum, missing a quad");

  /* count the actual number of quadrants per proc */
  P4EST_CHECK_ABORT (num_quadrants_in_proc[rank]
                     == p4est->local_num_quadrants,
                     "partition failed, wrong number of quadrants");

  for (t = p4est->first_local_tree; t <= p4est->last_local_tree; ++t) {
    tree = p4est_array_index (p4est->trees, t);
    for (q = 0; q < tree->quadrants->elem_count; ++q) {
      quad = p4est_array_index (tree->quadrants, q);
      user_data = (user_data_t *) quad->user_data;
      sum = quad->x + quad->y + quad->level;

      P4EST_CHECK_ABORT (user_data->a == t, "bad user_data, a");
      P4EST_CHECK_ABORT (user_data->sum == sum, "bad user_data, sum");
    }
  }

  /* clean up and exit */
  P4EST_FREE (num_quadrants_in_proc);
  p4est_destroy (p4est);
  p4est_connectivity_destroy (connectivity);
  p4est_memory_check ();

#ifdef HAVE_MPI
  mpiret = MPI_Finalize ();
  P4EST_CHECK_MPI (mpiret);
#endif

  return 0;
}

/* EOF simple.c */
