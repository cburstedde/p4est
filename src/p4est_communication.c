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
#include <p8est_algorithms.h>
#include <p8est_communication.h>
#include <p8est_bits.h>
#else
#include <p4est_algorithms.h>
#include <p4est_communication.h>
#include <p4est_bits.h>
#endif /* !P4_TO_P8 */
#include <sc_search.h>
#ifdef P4EST_HAVE_ZLIB
#include <zlib.h>
#endif

void
p4est_comm_parallel_env_create (p4est_t * p4est, sc_MPI_Comm mpicomm)
{
  int                 mpiret;

  /* duplicate MPI communicator */
  mpiret = sc_MPI_Comm_dup (mpicomm, &(p4est->mpicomm));
  SC_CHECK_MPI (mpiret);
  p4est->mpicomm_owned = 1;

  /* retrieve MPI information */
  mpiret = sc_MPI_Comm_size (mpicomm, &(p4est->mpisize));
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (mpicomm, &(p4est->mpirank));
  SC_CHECK_MPI (mpiret);
}

void
p4est_comm_parallel_env_free (p4est_t * p4est)
{
  int                 mpiret;

  /* free MPI communicator if it's owned */
  if (p4est->mpicomm_owned) {
    mpiret = sc_MPI_Comm_free (&(p4est->mpicomm));
    SC_CHECK_MPI (mpiret);
  }
  p4est->mpicomm = sc_MPI_COMM_NULL;
  p4est->mpicomm_owned = 0;

  /* set MPI information */
  p4est->mpisize = 0;
  p4est->mpirank = sc_MPI_UNDEFINED;
}

int
p4est_comm_parallel_env_is_null (p4est_t * p4est)
{
  return (p4est->mpicomm == sc_MPI_COMM_NULL);
}

void
p4est_comm_parallel_env_assign (p4est_t * p4est, sc_MPI_Comm mpicomm)
{
  int                 mpiret;

  /* check if input MPI communicator has same size and same rank order */
#ifdef P4EST_DEBUG
  {
    int                 result;

    mpiret = sc_MPI_Comm_compare (p4est->mpicomm, mpicomm, &result);
    SC_CHECK_MPI (mpiret);

    P4EST_ASSERT (result == sc_MPI_IDENT || result == sc_MPI_CONGRUENT);
  }
#endif

  /* free the current parallel environment */
  p4est_comm_parallel_env_free (p4est);

  /* assign MPI communicator of input, it is therefore not owned */
  p4est->mpicomm = mpicomm;
  p4est->mpicomm_owned = 0;

  /* retrieve MPI information */
  mpiret = sc_MPI_Comm_size (mpicomm, &(p4est->mpisize));
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (mpicomm, &(p4est->mpirank));
  SC_CHECK_MPI (mpiret);
}

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
p4est_comm_is_empty (p4est_t * p4est, int p)
{
  const p4est_gloidx_t *gfq;

  P4EST_ASSERT (p4est != NULL);
  P4EST_ASSERT (0 <= p && p < p4est->mpisize);

  gfq = p4est->global_first_quadrant;
  P4EST_ASSERT (gfq != NULL);

  return gfq[p] == gfq[p + 1];
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

#ifndef P4_TO_P8

void
p4est_transfer_fixed (p4est_t * dest, p4est_t * src,
                      p4est_transfer_comm_t which_comm, sc_MPI_Comm mpicomm,
                      int tag, void *dest_data,
                      const void *src_data, size_t data_size)
{
  p4est_transfer_context_t *tc;

  tc = p4est_transfer_fixed_begin (dest, src, which_comm, mpicomm, tag,
                                   dest_data, src_data, data_size);
  p4est_transfer_fixed_end (tc);
}

static void
p4est_transfer_determine_comm (p4est_transfer_context_t * tc,
                               p4est_t * dest, p4est_t * src,
                               p4est_transfer_comm_t which_comm,
                               sc_MPI_Comm mpicomm)
{
  int                 mpiret;
#ifdef P4EST_ENABLE_DEBUG
  int                 mpisize, mpirank;
#endif

  P4EST_ASSERT (tc != NULL);
  P4EST_ASSERT (dest != NULL && src != NULL);
  P4EST_ASSERT (p4est_is_valid (dest));
  P4EST_ASSERT (p4est_is_valid (src));
  P4EST_ASSERT (dest->mpisize == src->mpisize);
  P4EST_ASSERT (dest->mpirank == src->mpirank);
  P4EST_ASSERT (dest->global_num_quadrants == src->global_num_quadrants);

  tc->dest = dest;
  tc->src = src;
  tc->which_comm = which_comm;
  switch (which_comm) {
  case P4EST_TRANSFER_COMM_SRC:
    tc->mpicomm = src->mpicomm;
    break;
  case P4EST_TRANSFER_COMM_DEST:
    tc->mpicomm = dest->mpicomm;
    break;
  case P4EST_TRANSFER_COMM_SRC_DUP:
    mpiret = sc_MPI_Comm_dup (src->mpicomm, &tc->mpicomm);
    SC_CHECK_MPI (mpiret);
    break;
  case P4EST_TRANSFER_COMM_DEST_DUP:
    mpiret = sc_MPI_Comm_dup (dest->mpicomm, &tc->mpicomm);
    SC_CHECK_MPI (mpiret);
    break;
  case P4EST_TRANSFER_COMM_EXTERNAL:
    P4EST_ASSERT (mpicomm != sc_MPI_COMM_NULL);
#ifdef P4EST_ENABLE_DEBUG
    mpiret = sc_MPI_Comm_size (mpicomm, &mpisize);
    SC_CHECK_MPI (mpiret);
    P4EST_ASSERT (mpisize == src->mpisize);
    mpiret = sc_MPI_Comm_rank (mpicomm, &mpirank);
    SC_CHECK_MPI (mpiret);
    P4EST_ASSERT (mpirank == src->mpirank);
#endif
    tc->mpicomm = mpicomm;
    break;
  default:
    SC_ABORT_NOT_REACHED ();
    break;
  }
}

/** Given target, find index p such that gfq[p] <= target < gfq[p + 1].
 * \param [in] nmemb    Number of entries in array MINUS ONE.
 */
static int
p4est_bsearch_partition (p4est_gloidx_t target,
                         const p4est_gloidx_t * gfq, int nmemb)
{
  size_t              res;

  P4EST_ASSERT (nmemb > 0);
  P4EST_ASSERT (gfq[0] <= target);
  P4EST_ASSERT (target < gfq[nmemb]);

  res = sc_bsearch_range (&target, gfq, (size_t) nmemb,
                          sizeof (p4est_gloidx_t), p4est_gloidx_compare);
  P4EST_ASSERT (res < (size_t) nmemb);

  return (int) res;
}

p4est_transfer_context_t *
p4est_transfer_fixed_begin (p4est_t * dest, p4est_t * src,
                            p4est_transfer_comm_t which_comm,
                            sc_MPI_Comm mpicomm, int tag, void *dest_data,
                            const void *src_data, size_t data_size)
{
  p4est_transfer_context_t *tc;
  const p4est_gloidx_t *dest_gfq;
  const p4est_gloidx_t *src_gfq;
  int                 mpiret;
  int                 mpisize, mpirank;
  int                 q;
  int                 first_sender, last_sender;
  int                 first_receiver, last_receiver;
  char               *rb;
  char               *dest_cp, *src_cp;
  size_t              byte_len, cp_len;
  p4est_gloidx_t      dest_begin, dest_end;
  p4est_gloidx_t      src_begin, src_end;
  p4est_gloidx_t      gbegin, gend;
  sc_MPI_Request     *rq;

  /* setup context structure */
  tc = P4EST_ALLOC_ZERO (p4est_transfer_context_t, 1);
  p4est_transfer_determine_comm (tc, dest, src, which_comm, mpicomm);
  tc->tag = tag;
  tc->dest_data = dest_data;
  tc->src_data = src_data;
  tc->data_size = data_size;
  tc->variable = 0;

  /* there is nothing to do when there is no data */
  if (tc->data_size == 0) {
    return tc;
  }

  /* grab some p4est information */
  mpisize = dest->mpisize;
  mpirank = dest->mpirank;
  dest_gfq = dest->global_first_quadrant;
  dest_begin = dest_gfq[mpirank];
  dest_end = dest_gfq[mpirank + 1];
  src_gfq = src->global_first_quadrant;
  src_begin = src_gfq[mpirank];
  src_end = src_gfq[mpirank + 1];

  /* prepare data copy for local overlap */
  dest_cp = src_cp = NULL;
  cp_len = 0;

  /* figure out subset of processes to receive from */
  if (dest_begin < dest_end) {
    /* our process as the receiver is not empty */
    first_sender = p4est_bsearch_partition (dest_begin, src_gfq, mpisize);
    P4EST_ASSERT (0 <= first_sender && first_sender < mpisize);
    last_sender = p4est_bsearch_partition (dest_end - 1, src_gfq, mpisize);
    P4EST_ASSERT (first_sender <= last_sender && last_sender < mpisize);
    tc->num_senders = last_sender - first_sender + 1;
    P4EST_ASSERT (tc->num_senders > 0);

    /* go through sender processes and post receive calls */
    gend = dest_begin;
    rq = tc->recv_req = P4EST_ALLOC (sc_MPI_Request, tc->num_senders);
    rb = (char *) dest_data;
    for (q = first_sender; q <= last_sender; ++q) {
      /* prepare positions for the sender process q */
      gbegin = gend;
      gend = src_gfq[q + 1];
      if (gend > dest_end) {
        P4EST_ASSERT (q == last_sender);
        gend = dest_end;
      }
      P4EST_ASSERT (q == first_sender || q == last_sender ?
                    gbegin < gend : gbegin <= gend);
      byte_len = (size_t) (gend - gbegin) * tc->data_size;

      /* choose how to treat the sender process */
      if (byte_len == 0) {
        /* the sender process is empty; we need no message */
        P4EST_ASSERT (first_sender < q && q < last_sender);
        *rq++ = sc_MPI_REQUEST_NULL;
      }
      else if (q == mpirank) {
        /* for the same rank we remember pointers to memcpy */
        cp_len = byte_len;
        dest_cp = rb;
        *rq++ = sc_MPI_REQUEST_NULL;
      }
      else {
        mpiret = sc_MPI_Irecv (rb, byte_len, sc_MPI_BYTE, q,
                               tc->tag, tc->mpicomm, rq++);
        SC_CHECK_MPI (mpiret);
      }
      rb += byte_len;
    }
    P4EST_ASSERT (rb - (char *) tc->dest_data ==
                  (ptrdiff_t) ((dest_end - dest_begin) * tc->data_size));
  }

  /* figure out subset of processes to send to */
  if (src_begin < src_end) {
    /* our process as the sender is not empty */
    first_receiver = p4est_bsearch_partition (src_begin, dest_gfq, mpisize);
    P4EST_ASSERT (0 <= first_receiver && first_receiver < mpisize);
    last_receiver = p4est_bsearch_partition (src_end - 1, dest_gfq, mpisize);
    P4EST_ASSERT (first_receiver <= last_receiver && last_receiver < mpisize);
    tc->num_receivers = last_receiver - first_receiver + 1;
    P4EST_ASSERT (tc->num_receivers > 0);

    /* go through receiver processes and post send calls */
    gend = src_begin;
    rq = tc->send_req = P4EST_ALLOC (sc_MPI_Request, tc->num_receivers);
    rb = (char *) src_data;
    for (q = first_receiver; q <= last_receiver; ++q) {
      /* prepare positions for the receiver process q */
      gbegin = gend;
      gend = dest_gfq[q + 1];
      if (gend > src_end) {
        P4EST_ASSERT (q == last_receiver);
        gend = src_end;
      }
      P4EST_ASSERT (q == first_receiver || q == last_receiver ?
                    gbegin < gend : gbegin <= gend);
      byte_len = (size_t) (gend - gbegin) * tc->data_size;

      /* choose how to treat the receiver process */
      if (byte_len == 0) {
        /* the receiver process is empty; we need no message */
        P4EST_ASSERT (first_receiver < q && q < last_receiver);
        *rq++ = sc_MPI_REQUEST_NULL;
      }
      else if (q == mpirank) {
        /* for the same rank we remember pointers to memcpy */
        P4EST_ASSERT (cp_len == byte_len);
        src_cp = rb;
        *rq++ = sc_MPI_REQUEST_NULL;
      }
      else {
        mpiret = sc_MPI_Isend (rb, byte_len, sc_MPI_BYTE, q,
                               tc->tag, tc->mpicomm, rq++);
        SC_CHECK_MPI (mpiret);
      }
      rb += byte_len;
    }
    P4EST_ASSERT (rb - (char *) tc->src_data ==
                  (ptrdiff_t) ((src_end - src_begin) * tc->data_size));
  }

  /* copy the data that remains local */
  P4EST_ASSERT ((dest_cp == NULL) == (src_cp == NULL));
  if (cp_len > 0) {
    P4EST_ASSERT (dest_cp != NULL && src_cp != NULL);
    memcpy (dest_cp, src_cp, cp_len);
  }

  /* the rest goes into the p4est_transfer_fixed_end function */
  return tc;
}

static void
p4est_transfer_end (p4est_transfer_context_t * tc)
{
  int                 mpiret;

  P4EST_ASSERT (tc != NULL);

  if (tc->variable) {
    /* hmm what to do here */
  }

  /* wait for messages to complete and deallocate request buffers */
  if (tc->num_senders > 0) {
    mpiret = sc_MPI_Waitall (tc->num_senders, tc->recv_req,
                             sc_MPI_STATUSES_IGNORE);
    SC_CHECK_MPI (mpiret);
  }
  if (tc->num_receivers > 0) {
    mpiret = sc_MPI_Waitall (tc->num_receivers, tc->send_req,
                             sc_MPI_STATUSES_IGNORE);
    SC_CHECK_MPI (mpiret);
  }
  P4EST_FREE (tc->recv_req);
  P4EST_FREE (tc->send_req);

  /* free communicator if we have created it */
  switch (tc->which_comm) {
  case P4EST_TRANSFER_COMM_SRC_DUP:
  case P4EST_TRANSFER_COMM_DEST_DUP:
    mpiret = sc_MPI_Comm_free (&tc->mpicomm);
    SC_CHECK_MPI (mpiret);
    break;
  default:
    break;
  }

  /* the context must disappear too */
  tc->dest = tc->src = NULL;
  P4EST_FREE (tc);
}

void
p4est_transfer_fixed_end (p4est_transfer_context_t * tc)
{
  P4EST_ASSERT (tc != NULL);
  P4EST_ASSERT (!tc->variable);

  p4est_transfer_end (tc);
}

void
p4est_transfer_custom (p4est_t * dest, p4est_t * src,
                       p4est_transfer_comm_t which_comm, sc_MPI_Comm mpicomm,
                       int tag, void **dest_data, size_t ** dest_sizes,
                       const void *src_data, const size_t * src_sizes)
{
  p4est_transfer_context_t *tc;

  tc = p4est_transfer_custom_begin (dest, src, which_comm, mpicomm, tag,
                                    dest_data, dest_sizes,
                                    src_data, src_sizes);
  p4est_transfer_end (tc);
}

p4est_transfer_context_t *
p4est_transfer_custom_begin (p4est_t * dest, p4est_t * src,
                             p4est_transfer_comm_t which_comm,
                             sc_MPI_Comm mpicomm, int tag,
                             void **dest_data, size_t ** dest_sizes,
                             const void *src_data, const size_t * src_sizes)
{
  p4est_transfer_context_t *tc;

  P4EST_ASSERT (src->local_num_quadrants == 0 || src_sizes != NULL);
  P4EST_ASSERT (dest_data != NULL);
  P4EST_ASSERT (dest_sizes != NULL);

  tc = P4EST_ALLOC_ZERO (p4est_transfer_context_t, 1);
  p4est_transfer_determine_comm (tc, dest, src, which_comm, mpicomm);
  tc->tag = tag;
  tc->pdest_data = dest_data;
  tc->pdest_sizes = dest_sizes;
  tc->src_data = src_data;
  tc->src_sizes = src_sizes;
  tc->data_size = 0;
  tc->variable = 1;

  *dest_data = NULL;
  *dest_sizes = NULL;

  return tc;
}

void
p4est_transfer_custom_end (p4est_transfer_context_t * tc)
{
  P4EST_ASSERT (tc != NULL);
  P4EST_ASSERT (tc->variable);

  p4est_transfer_end (tc);
}

void
p4est_transfer_dest_data_free (p4est_t * dest,
                               void **dest_data, size_t ** dest_sizes)
{
  P4EST_ASSERT (p4est_is_valid (dest));
  P4EST_ASSERT (dest_data != NULL);
  P4EST_ASSERT (dest_sizes != NULL);

  P4EST_FREE (*dest_data);
  *dest_data = NULL;
  P4EST_FREE (*dest_sizes);
  *dest_sizes = NULL;
}

#endif
