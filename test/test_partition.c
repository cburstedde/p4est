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
  int64_t             sum;
}
user_data_t;

static int          weight_counter;
static int          weight_index;

static void
init_fn (p4est_t * p4est, int32_t which_tree, p4est_quadrant_t * quadrant)
{
  user_data_t        *data = quadrant->p.user_data;

  data->a = which_tree;
  data->sum = quadrant->x + quadrant->y + quadrant->level;
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

static int
weight_one (p4est_t * p4est, int32_t which_tree, p4est_quadrant_t * quadrant)
{
  return 1;
}

static int
weight_once (p4est_t * p4est, int32_t which_tree, p4est_quadrant_t * quadrant)
{
  if (weight_counter++ == weight_index) {
    return 1;
  }

  return 0;
}

int
main (int argc, char **argv)
{
  int                 rank = 0;
  int                 num_procs;
#ifdef P4EST_MPI
  int                 mpiret;
#endif
  MPI_Comm            mpicomm;
  p4est_t            *p4est, *copy;
  p4est_connectivity_t *connectivity;
  int                 i;
  int                 t, q;
  int32_t             num_quadrants_on_last;
  int32_t            *num_quadrants_in_proc;
  p4est_quadrant_t   *quad;
  p4est_tree_t       *tree;
  user_data_t        *user_data;
  int64_t             sum;
  unsigned            crc;

  mpicomm = MPI_COMM_NULL;
#ifdef P4EST_MPI
  mpiret = MPI_Init (&argc, &argv);
  P4EST_CHECK_MPI (mpiret);
  mpicomm = MPI_COMM_WORLD;
  mpiret = MPI_Comm_rank (mpicomm, &rank);
  P4EST_CHECK_MPI (mpiret);
#endif
  p4est_init (stdout, rank, NULL, NULL);

  /* create connectivity and forest structures */
  connectivity = p4est_connectivity_new_corner ();
  p4est = p4est_new (mpicomm, connectivity, sizeof (user_data_t), init_fn);

  num_procs = p4est->mpisize;

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

  /* check user data content */
  for (t = p4est->first_local_tree; t <= p4est->last_local_tree; ++t) {
    tree = p4est_array_index (p4est->trees, t);
    for (q = 0; q < tree->quadrants.elem_count; ++q) {
      quad = p4est_array_index (&tree->quadrants, q);
      user_data = (user_data_t *) quad->p.user_data;
      sum = quad->x + quad->y + quad->level;

      P4EST_CHECK_ABORT (user_data->a == t, "bad user_data, a");
      P4EST_CHECK_ABORT (user_data->sum == sum, "bad user_data, sum");
    }
  }

  /* do a weighted partition with uniform weights */
  p4est_partition (p4est, weight_one);
  P4EST_CHECK_ABORT (crc == p4est_checksum (p4est),
                     "bad checksum after uniformly weighted partition");

  /* copy the p4est */
  copy = p4est_copy (p4est, 1);
  P4EST_CHECK_ABORT (crc == p4est_checksum (copy), "bad checksum after copy");

  /* do a weighted partition with many zero weights */
  weight_counter = 0;
  weight_index = (rank == 1) ? 1342 : 0;
  p4est_partition (copy, weight_once);
  P4EST_CHECK_ABORT (crc == p4est_checksum (copy),
                     "bad checksum after unevenly weighted partition 1");

  /* do a weighted partition with many zero weights */
  weight_counter = 0;
  weight_index = 0;
  p4est_partition (copy, weight_once);
  P4EST_CHECK_ABORT (crc == p4est_checksum (copy),
                     "bad checksum after unevenly weighted partition 2");

  /* do a weighted partition with many zero weights */
  weight_counter = 0;
  weight_index =
    (rank == num_procs - 1) ? (copy->local_num_quadrants - 1) : 0;
  p4est_partition (copy, weight_once);
  P4EST_CHECK_ABORT (crc == p4est_checksum (copy),
                     "bad checksum after unevenly weighted partition 3");

  /* check user data content */
  for (t = copy->first_local_tree; t <= copy->last_local_tree; ++t) {
    tree = p4est_array_index (copy->trees, t);
    for (q = 0; q < tree->quadrants.elem_count; ++q) {
      quad = p4est_array_index (&tree->quadrants, q);
      user_data = (user_data_t *) quad->p.user_data;
      sum = quad->x + quad->y + quad->level;

      P4EST_CHECK_ABORT (user_data->a == t, "bad user_data, a");
      P4EST_CHECK_ABORT (user_data->sum == sum, "bad user_data, sum");
    }
  }

  /* clean up and exit */
  P4EST_FREE (num_quadrants_in_proc);
  p4est_destroy (p4est);
  p4est_destroy (copy);
  p4est_connectivity_destroy (connectivity);
  p4est_memory_check ();

#ifdef P4EST_MPI
  mpiret = MPI_Finalize ();
  P4EST_CHECK_MPI (mpiret);
#endif

  return 0;
}

/* EOF simple.c */
