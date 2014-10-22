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

#ifdef P4_TO_P8
#include <p8est_communication.h>
#include <p8est_bits.h>
#else
#include <p4est_communication.h>
#include <p4est_bits.h>
#endif /* !P4_TO_P8 */
#ifdef P4EST_HAVE_ZLIB
#include <zlib.h>
#endif

void
p4est_comm_count_quadrants (p4est_t * p4est)
{
  int                 mpiret;
  p4est_gloidx_t      qlocal = p4est->local_num_quadrants;
  p4est_gloidx_t     *global_first_quadrant = p4est->global_first_quadrant;
  int                 i;
  const int           num_procs = p4est->mpisize;

  global_first_quadrant[0] = 0;
  mpiret = sc_MPI_Allgather (&qlocal, 1, P4EST_MPI_GLOIDX,
                             global_first_quadrant + 1, 1, P4EST_MPI_GLOIDX,
                             p4est->mpicomm);
  SC_CHECK_MPI (mpiret);

  for (i = 0; i < num_procs; ++i) {
    global_first_quadrant[i + 1] += global_first_quadrant[i];
  }
  p4est->global_num_quadrants = global_first_quadrant[num_procs];
}

void
p4est_comm_global_partition (p4est_t * p4est, p4est_quadrant_t * first_quad)
{
  const int           num_procs = p4est->mpisize;
  const p4est_topidx_t num_trees = p4est->connectivity->num_trees;
  int                 i;
  int                 mpiret;
  const p4est_topidx_t first_tree = p4est->first_local_tree;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *quadrant;
  p4est_quadrant_t   *pi, input;

  SC_BZERO (&p4est->global_first_position[num_procs], 1);
  p4est->global_first_position[num_procs].level = P4EST_QMAXLEVEL;
  p4est->global_first_position[num_procs].p.which_tree = num_trees;

  SC_BZERO (&input, 1);
  if (first_tree < 0) {
    /* i don't have any quadrants, send negative values */
    P4EST_ASSERT (first_tree == -1 && p4est->last_local_tree == -2);
    input.x = -1;
    input.y = -1;
#ifdef P4_TO_P8
    input.z = -1;
#endif
  }
  else {
    /* send values corresponding to my first quadrant */
    if (first_quad != NULL) {
      tree = NULL;
      quadrant = first_quad;
    }
    else {
      tree = p4est_tree_array_index (p4est->trees, first_tree);
      quadrant = p4est_quadrant_array_index (&tree->quadrants, 0);
    }
    input.x = quadrant->x;
    input.y = quadrant->y;
#ifdef P4_TO_P8
    input.z = quadrant->z;
#endif
  }
  input.level = P4EST_QMAXLEVEL;
  input.p.which_tree = first_tree;
  mpiret = sc_MPI_Allgather (&input, (int) sizeof (p4est_quadrant_t),
                             sc_MPI_BYTE, p4est->global_first_position,
                             (int) sizeof (p4est_quadrant_t), sc_MPI_BYTE,
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
    P4EST_ASSERT (pi->p.which_tree >= 0 && pi->level == P4EST_QMAXLEVEL);
  }
}

void
p4est_comm_count_pertree (p4est_t * p4est, p4est_gloidx_t * pertree)
{
  const int           num_procs = p4est->mpisize;
  const int           rank = p4est->mpirank;
  const p4est_gloidx_t *gfq = p4est->global_first_quadrant;
  const p4est_quadrant_t *gfp = p4est->global_first_position;
  const p4est_topidx_t num_trees = p4est->connectivity->num_trees;
  int                 mpiret;
  int                 p;
  int                 mycount, c, addtomytree;
  int                *treecount, *treeoffset;
  p4est_topidx_t      t;
  p4est_locidx_t      recvbuf, sendbuf;
  p4est_gloidx_t     *mypertree;
  sc_MPI_Request      req_recv, req_send;
  sc_MPI_Status       status;
  p4est_tree_t       *tree;
#ifdef P4EST_ENABLE_DEBUG
  const p4est_quadrant_t *q;
#endif

  /* Tip off valgrind in case input array is too small */
  pertree[num_trees] = 0;

  /*
   * Determine which trees each rank will be counting.
   * A tree is counted by the processor that starts on its first quadrant,
   * even if this processor is empty.
   */
  treecount = P4EST_ALLOC (int, num_procs + 1);
  treeoffset = P4EST_ALLOC (int, num_procs + 1);
  p = 0;
  t = 0;
  treecount[0] = 1;
  treeoffset[0] = 0;
  for (;;) {
    /* Invariant: Rank p is the first that mentions tree t in gfp[p]
       and the ownership of t has been assigned to p or p - 1 */
    P4EST_ASSERT (gfp[p].p.which_tree == t);
    P4EST_ASSERT (p == 0 || gfp[p - 1].p.which_tree < t);
    do {
      treecount[++p] = 0;
    }
    while (gfp[p].p.which_tree == t);
    /* Assign the trees before the next first quadrant */
    P4EST_ASSERT (t < gfp[p].p.which_tree);
    for (++t; t < gfp[p].p.which_tree; ++t) {
      ++treecount[p - 1];
    }
    if (t < num_trees) {
      P4EST_ASSERT (p < num_procs);
      /* Check if the processor has the beginning of the tree */
      if (gfp[p].x == 0 && gfp[p].y == 0
#ifdef P4_TO_P8
          && gfp[p].z == 0
#endif
        ) {
        ++treecount[p];
      }
      else {
        ++treecount[p - 1];
      }
    }
    else {
      while (p < num_procs) {
        treecount[++p] = 0;
      }
      break;
    }
  }
  P4EST_ASSERT (p == num_procs);
  P4EST_ASSERT (t == num_trees);
  P4EST_ASSERT (treecount[num_procs] == 0);
  for (p = 0; p < num_procs; ++p) {
    treeoffset[p + 1] = treeoffset[p] + treecount[p];
  }
  P4EST_ASSERT ((p4est_topidx_t) treeoffset[num_procs] == num_trees);
  mycount = treecount[rank];
#ifdef P4EST_ENABLE_DEBUG
  P4EST_ASSERT (p4est->first_local_tree <= treeoffset[rank]);
  P4EST_ASSERT (gfq[rank + 1] - gfq[rank] ==
                (p4est_gloidx_t) p4est->local_num_quadrants);
  for (c = 0; c < mycount; ++c) {
    t = (p4est_topidx_t) (treeoffset[rank] + c);
    P4EST_ASSERT (rank == 0 || gfp[rank - 1].p.which_tree < t);
    if (p4est->local_num_quadrants > 0) {
      tree = p4est_tree_array_index (p4est->trees, t);
      q = p4est_quadrant_array_index (&tree->quadrants, 0);
    }
    else {
      q = gfp + rank;
    }
    P4EST_ASSERT (q->x == 0 && q->y == 0);
#ifdef P4_TO_P8
    P4EST_ASSERT (q->z == 0);
#endif
  }
#endif

  /* Go through trees this rank is responsible for and collect information */
  recvbuf = sendbuf = -1;
  addtomytree = -1;
  mypertree = P4EST_ALLOC (p4est_gloidx_t, mycount);
  for (c = 0; c < mycount; ++c) {
    /* Rank owns at least the first quadrant on this tree */
    t = (p4est_topidx_t) (treeoffset[rank] + c);
    tree = p4est_tree_array_index (p4est->trees, t);
    mypertree[c] = (p4est_gloidx_t) tree->quadrants.elem_count;
    if (c == mycount - 1) {
      /* Only the last tree in the counted list may be partially owned */
      for (p = rank + 1; p < num_procs && treecount[p] == 0; ++p) {
        P4EST_ASSERT (p < num_procs);
      }
      mypertree[c] += gfq[p] - gfq[rank + 1];
      if (gfp[p].p.which_tree == t) {
        P4EST_ASSERT (p < num_procs);
        /* Processor p has part of this tree too and needs to tell me */
        mpiret = sc_MPI_Irecv (&recvbuf, 1, P4EST_MPI_LOCIDX, p,
                               P4EST_COMM_COUNT_PERTREE, p4est->mpicomm,
                               &req_recv);
        SC_CHECK_MPI (mpiret);
        addtomytree = c;
      }
      else {
        P4EST_ASSERT (p <= num_procs);
        P4EST_ASSERT (gfp[p].p.which_tree == t + 1);
      }
    }
  }
  if (mycount > 0 && (t = gfp[rank].p.which_tree) < treeoffset[rank]) {
    /* Send information to processor that counts my first local quadrants */
    P4EST_ASSERT (rank > 0 && p4est->first_local_tree == t);
    tree = p4est_tree_array_index (p4est->trees, t);
    /* Always below the 32bit limit since this is processor-local data */
    sendbuf = (p4est_locidx_t) tree->quadrants.elem_count;
    for (p = rank - 1; treecount[p] == 0; --p) {
      P4EST_ASSERT (p > 0);
    }
    mpiret = sc_MPI_Isend (&sendbuf, 1, P4EST_MPI_LOCIDX, p,
                           P4EST_COMM_COUNT_PERTREE, p4est->mpicomm,
                           &req_send);
    SC_CHECK_MPI (mpiret);
  }

  /* Complete MPI operations and cumulative count */
  if (addtomytree >= 0) {
    mpiret = sc_MPI_Wait (&req_recv, &status);
    SC_CHECK_MPI (mpiret);
    mypertree[addtomytree] += (p4est_gloidx_t) recvbuf;
  }
  pertree[0] = 0;
  mpiret = sc_MPI_Allgatherv (mypertree, mycount, P4EST_MPI_GLOIDX,
                              pertree + 1, treecount, treeoffset,
                              P4EST_MPI_GLOIDX, p4est->mpicomm);
  SC_CHECK_MPI (mpiret);
  for (c = 0; c < (int) num_trees; ++c) {
    pertree[c + 1] += pertree[c];
  }
  if (sendbuf >= 0) {
    mpiret = sc_MPI_Wait (&req_send, &status);
    SC_CHECK_MPI (mpiret);
  }

  /* Clean up */
  P4EST_FREE (treecount);
  P4EST_FREE (treeoffset);
  P4EST_FREE (mypertree);
}

int
p4est_comm_is_owner (p4est_t * p4est, p4est_locidx_t which_tree,
                     const p4est_quadrant_t * q, int rank)
{
  p4est_topidx_t      ctree;
  p4est_quadrant_t    cur;
  const p4est_quadrant_t *global_first_position =
    p4est->global_first_position;

  cur.level = P4EST_QMAXLEVEL;
  P4EST_ASSERT (0 <= which_tree &&
                which_tree < p4est->connectivity->num_trees);
  P4EST_ASSERT (0 <= rank && rank < p4est->mpisize);
  P4EST_ASSERT (p4est_quadrant_is_node (q, 1) || p4est_quadrant_is_valid (q));

  /* check if q is on a lower processor than guess */
  ctree = global_first_position[rank].p.which_tree;
  cur.x = global_first_position[rank].x;
  cur.y = global_first_position[rank].y;
#ifdef P4_TO_P8
  cur.z = global_first_position[rank].z;
#endif
  if (which_tree < ctree ||
      (which_tree == ctree &&
       (p4est_quadrant_compare (q, &cur) < 0 &&
        (q->x != cur.x || q->y != cur.y
#ifdef P4_TO_P8
         || q->z != cur.z
#endif
        )))) {
    return 0;
  }

  /* check if q is on a higher processor than guess */
  ctree = global_first_position[rank + 1].p.which_tree;
  cur.x = global_first_position[rank + 1].x;
  cur.y = global_first_position[rank + 1].y;
#ifdef P4_TO_P8
  cur.z = global_first_position[rank + 1].z;
#endif
  if (which_tree > ctree ||
      (which_tree == ctree &&
       (p4est_quadrant_compare (&cur, q) <= 0 ||
        (q->x == cur.x && q->y == cur.y
#ifdef P4_TO_P8
         && q->z == cur.z
#endif
        )))) {
    return 0;
  }

  return 1;
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
  P4EST_ASSERT (p4est_quadrant_is_node (q, 1) || p4est_quadrant_is_valid (q));

  proc_low = 0;
  proc_high = num_procs - 1;
  cur.level = P4EST_QMAXLEVEL;

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
                      int full_tree[], int tree_contact[],
                      const p4est_quadrant_t ** pfirst_pos,
                      const p4est_quadrant_t ** pnext_pos)
{
  const p4est_quadrant_t *first_pos, *next_pos;
  p4est_connectivity_t *conn = p4est->connectivity;
  int                 face;

  P4EST_ASSERT (p4est->first_local_tree <= which_tree);
  P4EST_ASSERT (which_tree <= p4est->last_local_tree);

  first_pos = &p4est->global_first_position[p4est->mpirank];
  P4EST_ASSERT (first_pos->level == P4EST_QMAXLEVEL);
  full_tree[0] = (which_tree > p4est->first_local_tree ||
                  (first_pos->x == 0 && first_pos->y == 0
#ifdef P4_TO_P8
                   && first_pos->z == 0
#endif
                  ));

  next_pos = &p4est->global_first_position[p4est->mpirank + 1];
  P4EST_ASSERT (next_pos->level == P4EST_QMAXLEVEL);
  full_tree[1] = (which_tree < p4est->last_local_tree ||
                  (next_pos->x == 0 && next_pos->y == 0
#ifdef P4_TO_P8
                   && next_pos->z == 0
#endif
                  ));

  if (tree_contact != NULL) {
    for (face = 0; face < P4EST_FACES; ++face) {
      tree_contact[face] =
        (conn->tree_to_tree[P4EST_FACES * which_tree + face] != which_tree
         || (int) conn->tree_to_face[P4EST_FACES * which_tree + face] !=
         face);
    }
  }

  if (pfirst_pos != NULL) {
    *pfirst_pos = first_pos;
  }
  if (pnext_pos != NULL) {
    *pnext_pos = next_pos;
  }
}

int
p4est_comm_neighborhood_owned (p4est_t * p4est, p4est_locidx_t which_tree,
                               int full_tree[], int tree_contact[],
                               p4est_quadrant_t * q)
{
  const p4est_qcoord_t qh = P4EST_QUADRANT_LEN (q->level);
  const int           rank = p4est->mpirank;
  int                 n0_proc, n1_proc;
  p4est_quadrant_t    n0, n1;

  if (full_tree[0] && full_tree[1]) {
    /* need only to consider boundary quadrants */
    if (!((tree_contact[0] && q->x == 0) ||
          (tree_contact[1] && q->x == P4EST_ROOT_LEN - qh) ||
          (tree_contact[2] && q->y == 0) ||
          (tree_contact[3] && q->y == P4EST_ROOT_LEN - qh) ||
#ifdef P4_TO_P8
          (tree_contact[4] && q->z == 0) ||
          (tree_contact[5] && q->z == P4EST_ROOT_LEN - qh) ||
#endif
          0)) {
      return 1;
    }
  }
  else {
    /* test lowest and highest neighbor first */
    n0.x = q->x - qh;
    n0.y = q->y - qh;
#ifdef P4_TO_P8
    n0.z = q->z - qh;
#endif
    n0.level = q->level;
    if (n0.x >= 0 && n0.y >= 0
#ifdef P4_TO_P8
        && n0.z >= 0
#endif
      ) {
      n1.x = q->x + qh;
      n1.y = q->y + qh;
#ifdef P4_TO_P8
      n1.z = q->z + qh;
#endif
      n1.level = q->level;
      if (n1.x < P4EST_ROOT_LEN && n1.y < P4EST_ROOT_LEN
#ifdef P4_TO_P8
          && n1.z < P4EST_ROOT_LEN
#endif
        ) {
        n0_proc = p4est_comm_find_owner (p4est, which_tree, &n0, rank);
        if (n0_proc == rank) {
          p4est_quadrant_last_descendant (&n1, &n0, P4EST_QMAXLEVEL);
          n1_proc = p4est_comm_find_owner (p4est, which_tree, &n0, rank);
          if (n1_proc == rank) {
            return 1;
          }
        }
      }
    }
  }

  return 0;
}

int
p4est_comm_sync_flag (p4est_t * p4est, int flag, sc_MPI_Op operation)
{
  int8_t              lbyte, gbyte;
  int                 mpiret;

  P4EST_ASSERT (operation == sc_MPI_BAND || operation == sc_MPI_BOR);

  lbyte = (int8_t) (flag ? 1 : 0);
  mpiret = sc_MPI_Allreduce (&lbyte, &gbyte, 1, sc_MPI_BYTE, operation,
                             p4est->mpicomm);
  SC_CHECK_MPI (mpiret);

  return (int) gbyte;
}

unsigned
p4est_comm_checksum (p4est_t * p4est, unsigned local_crc, size_t local_bytes)
{
#ifdef P4EST_HAVE_ZLIB
  uLong               crc = (uLong) local_crc;

#ifdef P4EST_ENABLE_MPI
  int                 mpiret;
  int                 p;
  uint64_t            send[2];
  uint64_t           *gather;

  send[0] = (uint64_t) local_crc;
  send[1] = (uint64_t) local_bytes;
  gather = NULL;
  if (p4est->mpirank == 0) {
    gather = P4EST_ALLOC (uint64_t, 2 * p4est->mpisize);
  }
  mpiret = MPI_Gather (send, 2, MPI_LONG_LONG_INT,
                       gather, 2, MPI_LONG_LONG_INT, 0, p4est->mpicomm);
  SC_CHECK_MPI (mpiret);

  if (p4est->mpirank == 0) {
    for (p = 1; p < p4est->mpisize; ++p) {
      crc = adler32_combine (crc, (uLong) gather[2 * p + 0],
                             (z_off_t) gather[2 * p + 1]);
    }
    P4EST_FREE (gather);
  }
  else {
    crc = 0;
  }
#endif /* P4EST_ENABLE_MPI */

  return (unsigned) crc;
#else
  sc_abort_collective
    ("Configure did not find a recent enough zlib.  Abort.\n");

  return 0;
#endif /* !P4EST_HAVE_ZLIB */
}
