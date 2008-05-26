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

#include <p4est_base.h>
#include <p4est_communication.h>
#include <p4est_algorithms.h>

void
p4est_comm_count_quadrants (p4est_t * p4est)
{
#ifdef P4EST_MPI
  int                 mpiret;
#endif
  p4est_gloidx_t      qlocal = p4est->local_num_quadrants;
  p4est_gloidx_t     *global_last_quad_index = p4est->global_last_quad_index;
  int                 i;
  const int           rank = p4est->mpirank;
  const int           num_procs = p4est->mpisize;

  global_last_quad_index[rank] = qlocal;
#ifdef P4EST_MPI
  if (p4est->mpicomm != MPI_COMM_NULL) {
    mpiret = MPI_Allgather (&qlocal, 1, P4EST_MPI_GLOIDX,
                            global_last_quad_index, 1, P4EST_MPI_GLOIDX,
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
#ifdef P4EST_MPI
  int                 i;
  int                 mpiret;
  const int32_t       first_tree = p4est->first_local_tree;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *quadrant;
  p4est_position_t   *pi, input;
#endif

  SC_BZERO (&p4est->global_first_position[0], 1);
  SC_BZERO (&p4est->global_first_position[num_procs], 1);

  p4est->global_first_position[num_procs].which_tree = num_trees;
  p4est->global_first_position[num_procs].x = 0;
  p4est->global_first_position[num_procs].y = 0;

#ifdef P4EST_MPI
  if (p4est->mpicomm != MPI_COMM_NULL) {
    SC_BZERO (&input, 1);
    if (first_tree < 0) {
      /* i don't have any quadrants, send negative values */
      P4EST_ASSERT (first_tree == -1 && p4est->last_local_tree == -2);
      input.which_tree = -1;
      input.x = input.y = -1;
    }
    else {
      /* send values corresponding to my first quadrant */
      tree = sc_array_index (p4est->trees, first_tree);
      quadrant = sc_array_index (&tree->quadrants, 0);
      input.which_tree = first_tree;
      input.x = quadrant->x;
      input.y = quadrant->y;
    }
    mpiret = MPI_Allgather (&input, sizeof (p4est_position_t), MPI_BYTE,
                            p4est->global_first_position,
                            sizeof (p4est_position_t), MPI_BYTE,
                            p4est->mpicomm);
    P4EST_CHECK_MPI (mpiret);

    /* correct for processors that don't have any quadrants */
    for (i = num_procs - 1; i >= 0; --i) {
      pi = &p4est->global_first_position[i];
      if (pi->which_tree < 0) {
        P4EST_ASSERT (pi->x == -1 && pi->y == -1);
        memcpy (pi, pi + 1, sizeof (p4est_position_t));
      }
      P4EST_ASSERT (pi->which_tree >= 0 && pi->x >= 0 && pi->y >= 0);
    }
  }
#endif
}

int
p4est_comm_find_owner (p4est_t * p4est, p4est_locidx_t which_tree,
                       const p4est_quadrant_t * q, int guess)
{
  const int           num_procs = p4est->mpisize;
  const p4est_position_t *global_first_position =
    p4est->global_first_position;
  int                 proc_low, proc_high;
  p4est_locidx_t      ctree;
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
    ctree = global_first_position[guess].which_tree;
    cur.x = global_first_position[guess].x;
    cur.y = global_first_position[guess].y;
    if (which_tree < ctree ||
        (which_tree == ctree &&
         (p4est_quadrant_compare (q, &cur) < 0 &&
          (q->x != cur.x || q->y != cur.y)))) {
      proc_high = guess - 1;
      guess = (proc_low + proc_high + 1) / 2;
      continue;
    }

    /* check if q is on a higher processor than guess */
    ctree = global_first_position[guess + 1].which_tree;
    cur.x = global_first_position[guess + 1].x;
    cur.y = global_first_position[guess + 1].y;
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
  P4EST_ASSERT (memcmp (&global_first_position[guess],
                        &global_first_position[guess + 1],
                        sizeof (p4est_position_t)) != 0);
  return guess;
}

/* EOF p4est_communication.h */
