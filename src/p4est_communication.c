/*
  This file is part of p4est.
  p4est is a C library to manage a parallel collection of quadtrees and/or
  octrees.

  Copyright (C) 2007,2008 Carsten Burstedde, Lucas Wilcox.

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

#include <p4est_communication.h>
#include <p4est_algorithms.h>
#include <p4est_base.h>

void
p4est_comm_count_quadrants (p4est_t * p4est)
{
#ifdef HAVE_MPI
  int                 mpiret;
#endif
  int64_t             qlocal = p4est->local_num_quadrants;
  int64_t            *global_last_quad_index = p4est->global_last_quad_index;
  int32_t             rank = p4est->mpirank;
  int32_t             num_procs = p4est->mpisize;
  int                 i;

  global_last_quad_index[rank] = qlocal;
#ifdef HAVE_MPI
  if (p4est->mpicomm != MPI_COMM_NULL) {
    mpiret = MPI_Allgather (&qlocal, 1, MPI_LONG_LONG,
                            global_last_quad_index, 1, MPI_LONG_LONG,
                            p4est->mpicomm);
    P4EST_CHECK_MPI (mpiret);
  }
#endif

  /* Subtract 1 from the first index since we are zero based */
  --global_last_quad_index[0];
  for (i = 1; i < num_procs; ++i) {
    global_last_quad_index[i] += global_last_quad_index[i - 1];
  }
  /* Add 1 to the last index since we are zero based */
  p4est->global_num_quadrants = global_last_quad_index[num_procs - 1] + 1;
}

void
p4est_comm_global_partition (p4est_t * p4est)
{
  const int           num_procs = p4est->mpisize;
  const int32_t       num_trees = p4est->connectivity->num_trees;
#ifdef HAVE_MPI
  int                 i;
  int                 mpiret;
  int32_t             input[3];
  int32_t            *pi;
  const int32_t       first_tree = p4est->first_local_tree;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *quadrant;
#endif

  p4est->global_first_indices[0] = 0;
  p4est->global_first_indices[1] = 0;
  p4est->global_first_indices[2] = 0;

#ifdef HAVE_MPI
  if (p4est->mpicomm != MPI_COMM_NULL) {
    if (first_tree < 0) {
      /* i don't have any quadrants, send negative values */
      P4EST_ASSERT (first_tree == -1 && p4est->last_local_tree == -2);
      input[0] = input[1] = input[2] = -1;
    }
    else {
      /* send values corresponding to my first quadrant */
      tree = p4est_array_index (p4est->trees, first_tree);
      quadrant = p4est_array_index (&tree->quadrants, 0);
      input[0] = first_tree;
      input[1] = quadrant->x;
      input[2] = quadrant->y;
    }
    mpiret = MPI_Allgather (input, 3, MPI_INT,
                            p4est->global_first_indices, 3, MPI_INT,
                            p4est->mpicomm);
    P4EST_CHECK_MPI (mpiret);

    /* correct for processors that don't have any quadrants */
    for (i = num_procs - 1; i >= 0; --i) {
      pi = &p4est->global_first_indices[3 * i];
      if (pi[0] < 0) {
        P4EST_ASSERT (pi[1] < 0 && pi[2] < 0);
        memcpy (pi, &pi[3], 3 * sizeof (*pi));
      }
      P4EST_ASSERT (pi[0] >= 0 && pi[1] >= 0 && pi[2] >= 0);
    }
  }
#endif

  p4est->global_first_indices[3 * num_procs + 0] = num_trees;
  p4est->global_first_indices[3 * num_procs + 1] = 0;
  p4est->global_first_indices[3 * num_procs + 2] = 0;
}

int32_t
p4est_comm_find_owner (p4est_t * p4est, int32_t which_tree,
                       const p4est_quadrant_t * q, int32_t guess)
{
  const int           num_procs = p4est->mpisize;
  const int32_t      *global_first_indices = p4est->global_first_indices;
  int32_t             proc_low, proc_high;
  int32_t             ctree;
  p4est_quadrant_t    cur;

  proc_low = 0;
  proc_high = num_procs - 1;
  cur.level = P4EST_MAXLEVEL;

  for (;;) {
    P4EST_ASSERT (proc_low <= proc_high);
    P4EST_ASSERT (0 <= proc_low && proc_low < num_procs);
    P4EST_ASSERT (0 <= proc_high && proc_high < num_procs);
    P4EST_ASSERT (proc_low <= guess && guess <= proc_high);

    /* check if q is on a lower processor than guess */
    ctree = global_first_indices[3 * guess + 0];
    cur.x = global_first_indices[3 * guess + 1];
    cur.y = global_first_indices[3 * guess + 2];
    if (which_tree < ctree ||
        (which_tree == ctree &&
         (p4est_quadrant_compare (q, &cur) < 0 &&
          (q->x != cur.x || q->y != cur.y)))) {
      proc_high = guess - 1;
      guess = (proc_low + proc_high + 1) / 2;
      continue;
    }

    /* check if q is on a higher processor than guess */
    ctree = global_first_indices[3 * (guess + 1) + 0];
    cur.x = global_first_indices[3 * (guess + 1) + 1];
    cur.y = global_first_indices[3 * (guess + 1) + 2];
    if (which_tree > ctree ||
        (which_tree == ctree &&
         (p4est_quadrant_compare (&cur, q) <= 0 ||
          (q->x == cur.x && q->y == cur.y)))) {
      proc_low = guess + 1;
      guess = (proc_low + proc_high) / 2;
      continue;
    }

    /* otherwise guess is the correct processor */
    break;
  }

  /* make sure we found a valid processor with nonzero quadrant count */
  P4EST_ASSERT (0 <= guess && guess < num_procs);
  P4EST_ASSERT (memcmp (&global_first_indices[3 * guess],
                        &global_first_indices[3 * (guess + 1)],
                        3 * sizeof (int32_t)) != 0);
  return guess;
}

/* EOF p4est_communication.h */
