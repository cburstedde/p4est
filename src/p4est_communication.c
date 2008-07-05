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

#ifdef P4_TO_P8
#include <p8est_communication.h>
#include <p8est_bits.h>
#else
#include <p4est_communication.h>
#include <p4est_bits.h>
#endif /* !P4_TO_P8 */

void
p4est_comm_count_quadrants (p4est_t * p4est)
{
  int                 mpiret;
  p4est_gloidx_t      qlocal = p4est->local_num_quadrants;
  p4est_gloidx_t     *global_last_quad_index = p4est->global_last_quad_index;
  int                 i;
  const int           rank = p4est->mpirank;
  const int           num_procs = p4est->mpisize;

  global_last_quad_index[rank] = qlocal;

  if (p4est->mpicomm != MPI_COMM_NULL) {
    mpiret = MPI_Allgather (&qlocal, 1, P4EST_MPI_GLOIDX,
                            global_last_quad_index, 1, P4EST_MPI_GLOIDX,
                            p4est->mpicomm);
    SC_CHECK_MPI (mpiret);
  }

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
  const p4est_topidx_t num_trees = p4est->connectivity->num_trees;
  int                 i;
  int                 mpiret;
  const p4est_topidx_t first_tree = p4est->first_local_tree;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *quadrant;
  p4est_quadrant_t   *pi, input;

  SC_BZERO (&p4est->global_first_position[0], 1);
  SC_BZERO (&p4est->global_first_position[num_procs], 1);

  p4est->global_first_position[num_procs].x = 0;
  p4est->global_first_position[num_procs].y = 0;
#ifdef P4_TO_P8
  p4est->global_first_position[num_procs].z = 0;
#endif
  p4est->global_first_position[num_procs].level = P4EST_MAXLEVEL;
  p4est->global_first_position[num_procs].p.which_tree = num_trees;

  if (p4est->mpicomm != MPI_COMM_NULL) {
    SC_BZERO (&input, 1);
    if (first_tree < 0) {
      /* i don't have any quadrants, send negative values */
      P4EST_ASSERT (first_tree == -1 && p4est->last_local_tree == -2);
      input.x = input.y = -1;
#ifdef P4_TO_P8
      input.z = -1;
#endif
      input.level = P4EST_MAXLEVEL;
      input.p.which_tree = -1;
    }
    else {
      /* send values corresponding to my first quadrant */
      tree = sc_array_index (p4est->trees, first_tree);
      quadrant = sc_array_index (&tree->quadrants, 0);
      input.x = quadrant->x;
      input.y = quadrant->y;
#ifdef P4_TO_P8
      input.z = quadrant->z;
#endif
      input.level = P4EST_MAXLEVEL;
      input.p.which_tree = first_tree;
    }
    mpiret = MPI_Allgather (&input, (int) sizeof (p4est_quadrant_t), MPI_BYTE,
                            p4est->global_first_position,
                            (int) sizeof (p4est_quadrant_t), MPI_BYTE,
                            p4est->mpicomm);
    SC_CHECK_MPI (mpiret);

    /* correct for processors that don't have any quadrants */
    for (i = num_procs - 1; i >= 0; --i) {
      pi = &p4est->global_first_position[i];
      if (pi->p.which_tree < 0) {
        P4EST_ASSERT (pi->x == -1 && pi->y == -1);
#ifdef P4_TO_P8
        P4EST_ASSERT (pi->z == -1);
#endif
        memcpy (pi, pi + 1, sizeof (p4est_quadrant_t));
      }
      P4EST_ASSERT (pi->x >= 0 && pi->y >= 0);
#ifdef P4_TO_P8
      P4EST_ASSERT (pi->z >= 0);
#endif
      P4EST_ASSERT (pi->p.which_tree >= 0 && pi->level == P4EST_MAXLEVEL);
    }
  }
}

int
p4est_comm_find_owner (p4est_t * p4est, p4est_locidx_t which_tree,
                       const p4est_quadrant_t * q, int guess)
{
  const int           num_procs = p4est->mpisize;
  const p4est_quadrant_t *global_first_position =
    p4est->global_first_position;
  int                 proc_low, proc_high;
  p4est_topidx_t      ctree;
  p4est_quadrant_t    cur;

  P4EST_ASSERT (0 <= which_tree &&
                which_tree < p4est->connectivity->num_trees);
  P4EST_ASSERT (p4est_quadrant_is_valid (q));

  proc_low = 0;
  proc_high = num_procs - 1;
  cur.level = P4EST_MAXLEVEL;

  for (;;) {
    P4EST_ASSERT (proc_low <= proc_high);
    P4EST_ASSERT (0 <= proc_low && proc_low < num_procs);
    P4EST_ASSERT (0 <= proc_high && proc_high < num_procs);
    P4EST_ASSERT (proc_low <= guess && guess <= proc_high);

    /* check if q is on a lower processor than guess */
    ctree = global_first_position[guess].p.which_tree;
    cur.x = global_first_position[guess].x;
    cur.y = global_first_position[guess].y;
#ifdef P4_TO_P8
    cur.z = global_first_position[guess].z;
#endif
    if (which_tree < ctree ||
        (which_tree == ctree &&
         (p4est_quadrant_compare (q, &cur) < 0 &&
          (q->x != cur.x || q->y != cur.y
#ifdef P4_TO_P8
           || q->z != cur.z
#endif
          )))) {
      proc_high = guess - 1;
      guess = (proc_low + proc_high + 1) / 2;
      continue;
    }

    /* check if q is on a higher processor than guess */
    ctree = global_first_position[guess + 1].p.which_tree;
    cur.x = global_first_position[guess + 1].x;
    cur.y = global_first_position[guess + 1].y;
#ifdef P4_TO_P8
    cur.z = global_first_position[guess + 1].z;
#endif
    if (which_tree > ctree ||
        (which_tree == ctree &&
         (p4est_quadrant_compare (&cur, q) <= 0 ||
          (q->x == cur.x && q->y == cur.y
#ifdef P4_TO_P8
           && q->z == cur.z
#endif
          )))) {
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
                        sizeof (p4est_quadrant_t)) != 0);
  return guess;
}

void
p4est_comm_tree_info (p4est_t * p4est, p4est_locidx_t which_tree,
                      bool full_tree[], const p4est_quadrant_t ** pfirst_pos,
                      const p4est_quadrant_t ** pnext_pos)
{
  const p4est_quadrant_t *first_pos, *next_pos;

  P4EST_ASSERT (p4est->first_local_tree <= which_tree);
  P4EST_ASSERT (which_tree <= p4est->last_local_tree);

  first_pos = &p4est->global_first_position[p4est->mpirank];
  P4EST_ASSERT (first_pos->level == P4EST_MAXLEVEL);
  full_tree[0] = (which_tree > p4est->first_local_tree ||
                  (first_pos->x == 0 && first_pos->y == 0
#ifdef P4_TO_P8
                   && first_pos->z == 0
#endif
                  ));

  next_pos = &p4est->global_first_position[p4est->mpirank + 1];
  P4EST_ASSERT (next_pos->level == P4EST_MAXLEVEL);
  full_tree[1] = (which_tree < p4est->last_local_tree ||
                  (next_pos->x == 0 && next_pos->y == 0
#ifdef P4_TO_P8
                   && next_pos->z == 0
#endif
                  ));

  if (pfirst_pos != NULL) {
    *pfirst_pos = first_pos;
  }
  if (pnext_pos != NULL) {
    *pnext_pos = next_pos;
  }
}

bool
p4est_comm_sync_flag (p4est_t * p4est, bool flag, MPI_Op operation)
{
  int8_t              lbyte, gbyte;
  int                 mpiret;

  P4EST_ASSERT (operation == MPI_BAND || operation == MPI_BOR);

  if (p4est->mpicomm == MPI_COMM_NULL) {
    return flag;
  }

  lbyte = (int8_t) (flag ? 1 : 0);
  mpiret = MPI_Allreduce (&lbyte, &gbyte, 1, MPI_BYTE, operation,
                          p4est->mpicomm);
  SC_CHECK_MPI (mpiret);

  return (bool) gbyte;
}

/* EOF p4est_communication.c */
