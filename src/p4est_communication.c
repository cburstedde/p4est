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

#include <p4est_communication.h>
#include <p4est_algorithms.h>
#include <p4est_base.h>

void
p4est_comm_count_quadrants (p4est_t * p4est)
{
#ifdef HAVE_MPI
  int                 mpiret;
#endif
  int64_t             qlocal;

  qlocal = p4est->local_num_quadrants;
  p4est->global_num_quadrants = qlocal;

#ifdef HAVE_MPI
  if (p4est->mpicomm != MPI_COMM_NULL) {
    mpiret = MPI_Allreduce (&qlocal, &p4est->global_num_quadrants,
                            1, MPI_LONG_LONG, MPI_SUM, p4est->mpicomm);
    P4EST_CHECK_MPI (mpiret);
  }
#endif
}

void
p4est_comm_global_partition (p4est_t * p4est)
{
#ifdef HAVE_MPI
  int                 mpiret;
  int32_t             input[3];
  p4est_tree_t       *tree;
  p4est_quadrant_t   *quadrant;
#endif

  p4est->global_first_indices[0] = 0;
  p4est->global_first_indices[1] = 0;
  p4est->global_first_indices[2] = 0;

#ifdef HAVE_MPI
  if (p4est->mpicomm != MPI_COMM_NULL) {
    tree = p4est_array_index (p4est->trees, p4est->first_local_tree);
    quadrant = p4est_array_index (tree->quadrants, 0);

    input[0] = p4est->first_local_tree;
    input[1] = quadrant->x;
    input[2] = quadrant->y;
    mpiret = MPI_Allgather (input, 3, MPI_INT,
                            p4est->global_first_indices, 3, MPI_INT,
                            p4est->mpicomm);
    P4EST_CHECK_MPI (mpiret);
  }
#endif

  p4est->global_first_indices[3 * p4est->mpisize + 0] =
    p4est->connectivity->num_trees;
  p4est->global_first_indices[3 * p4est->mpisize + 1] = 0;
  p4est->global_first_indices[3 * p4est->mpisize + 2] = 0;
}

int32_t
p4est_comm_find_owner (p4est_t * p4est, int32_t which_tree,
                       const p4est_quadrant_t * q, int32_t guess)
{
  int32_t             proc_low, proc_high;
  int32_t             ctree;
  p4est_quadrant_t    cur;

  proc_low = 0;
  proc_high = p4est->mpisize - 1;
  cur.level = P4EST_MAXLEVEL;
  
  for (;;) {
    P4EST_ASSERT (proc_low <= proc_high);
    P4EST_ASSERT (0 <= proc_low && proc_low < p4est->mpisize);
    P4EST_ASSERT (0 <= proc_high && proc_high < p4est->mpisize);
    P4EST_ASSERT (proc_low <= guess && guess <= proc_high);

    /* check if q is on a lower processor than guess */
    ctree = p4est->global_first_indices[3 * guess + 0];
    cur.x = p4est->global_first_indices[3 * guess + 1];
    cur.y = p4est->global_first_indices[3 * guess + 2];
    if (which_tree < ctree ||
        (which_tree == ctree &&
         (p4est_quadrant_compare (q, &cur) < 0 &&
          (q->x != cur.x || q->y != cur.y)))) {
      proc_high = guess - 1;
      guess = (proc_low + proc_high + 1) / 2;
      continue;
    }

    /* check if q is on a higher processor than guess */
    ctree = p4est->global_first_indices[3 * (guess + 1) + 0];
    cur.x = p4est->global_first_indices[3 * (guess + 1) + 1];
    cur.y = p4est->global_first_indices[3 * (guess + 1) + 2];
    if (which_tree > ctree ||
        (which_tree == ctree &&
         (p4est_quadrant_compare (&cur, q) <= 0))) {
      proc_low = guess + 1;
      guess = (proc_low + proc_high) / 2;
      continue;
    }

    /* otherwise guess is the correct processor */
    break;
  }

  return guess;
}

/* EOF p4est_communication.h */
