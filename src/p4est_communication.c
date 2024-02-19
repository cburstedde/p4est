/*
  This file is part of p4est.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2010 The University of Texas System
  Additional copyright (C) 2011 individual authors
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

int
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

void
p4est_comm_parallel_env_assign (p4est_t * p4est, sc_MPI_Comm mpicomm)
{
  /* set MPI communicator */
  p4est->mpicomm = mpicomm;
  p4est->mpicomm_owned = 0;

  /* retrieve MPI information */
  p4est_comm_parallel_env_get_info (p4est);
}

void
p4est_comm_parallel_env_duplicate (p4est_t * p4est)
{
  sc_MPI_Comm         mpicomm = p4est->mpicomm;
  int                 mpiret;

  /* duplicate MPI communicator */
  mpiret = sc_MPI_Comm_dup (mpicomm, &(p4est->mpicomm));
  SC_CHECK_MPI (mpiret);
  p4est->mpicomm_owned = 1;
}

void
p4est_comm_parallel_env_release (p4est_t * p4est)
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

void
p4est_comm_parallel_env_replace (p4est_t * p4est, sc_MPI_Comm mpicomm)
{
  /* check if input MPI communicator has same size and same rank order */
#ifdef P4EST_ENABLE_DEBUG
  {
    int                 mpiret, result;

    mpiret = sc_MPI_Comm_compare (p4est->mpicomm, mpicomm, &result);
    SC_CHECK_MPI (mpiret);

    P4EST_ASSERT (result == sc_MPI_IDENT || result == sc_MPI_CONGRUENT);
  }
#endif

  /* release the current parallel environment */
  p4est_comm_parallel_env_release (p4est);

  /* assign new MPI communicator */
  p4est_comm_parallel_env_assign (p4est, mpicomm);
}

void
p4est_comm_parallel_env_get_info (p4est_t * p4est)
{
  int                 mpiret;

  mpiret = sc_MPI_Comm_size (p4est->mpicomm, &(p4est->mpisize));
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (p4est->mpicomm, &(p4est->mpirank));
  SC_CHECK_MPI (mpiret);
}

int
p4est_comm_parallel_env_is_null (p4est_t * p4est)
{
  return (p4est->mpicomm == sc_MPI_COMM_NULL);
}

int
p4est_comm_parallel_env_reduce (p4est_t ** p4est_supercomm)
{
  return p4est_comm_parallel_env_reduce_ext (p4est_supercomm,
                                             sc_MPI_GROUP_NULL, 0, NULL);
}

int
p4est_comm_parallel_env_reduce_ext (p4est_t ** p4est_supercomm,
                                    sc_MPI_Group group_add,
                                    int add_to_beginning, int **ranks_subcomm)
{
  const char         *this_fn_name = "comm_parallel_env_reduce";
  p4est_t            *p4est = *p4est_supercomm;
  sc_MPI_Comm         mpicomm = p4est->mpicomm;
  int                 mpisize = p4est->mpisize;
  int                 mpiret;
  p4est_gloidx_t     *global_first_quadrant = p4est->global_first_quadrant;
  p4est_quadrant_t   *global_first_position = p4est->global_first_position;

  p4est_gloidx_t     *n_quadrants;
  int                *include;
  sc_MPI_Group        group, subgroup;
  sc_MPI_Comm         submpicomm;
  int                 submpisize, submpirank;
  int                *ranks, *subranks;
  int                 p;

  /* exit if MPI communicator cannot be reduced */
  if (mpisize == 1) {
    return 1;
  }

  /* create array of non-empty processes that will be included to sub-comm */
  n_quadrants = P4EST_ALLOC (p4est_gloidx_t, mpisize);
  include = P4EST_ALLOC (int, mpisize);
  submpisize = 0;
  for (p = 0; p < mpisize; p++) {
    n_quadrants[p] = global_first_quadrant[p + 1] - global_first_quadrant[p];
    if (global_first_quadrant[p] < global_first_quadrant[p + 1]) {
      include[submpisize++] = p;
    }
  }

  /* exit if reduction not possible */
  if (submpisize == mpisize) {
    P4EST_FREE (n_quadrants);
    P4EST_FREE (include);
    return 1;
  }

  /* create sub-group of non-empty processors */
  mpiret = sc_MPI_Comm_group (mpicomm, &group);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Group_incl (group, submpisize, include, &subgroup);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Group_free (&group);
  SC_CHECK_MPI (mpiret);
  P4EST_FREE (include);

  /* create sub-communicator */
  if (group_add != sc_MPI_GROUP_NULL) {
    sc_MPI_Group        group_union;

    /* create union with optional group */
    if (add_to_beginning) {
      mpiret = sc_MPI_Group_union (group_add, subgroup, &group_union);
    }
    else {
      mpiret = sc_MPI_Group_union (subgroup, group_add, &group_union);
    }
    SC_CHECK_MPI (mpiret);

    /* create sub-communicator */
    mpiret = sc_MPI_Comm_create (mpicomm, group_union, &submpicomm);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Group_free (&group_union);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Group_free (&subgroup);
    SC_CHECK_MPI (mpiret);
  }
  else {
    /* create sub-communicator */
    mpiret = sc_MPI_Comm_create (mpicomm, subgroup, &submpicomm);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Group_free (&subgroup);
    SC_CHECK_MPI (mpiret);
  }

  /* destroy p4est and exit if this rank is empty */
  if (submpicomm == sc_MPI_COMM_NULL) {
    /* destroy */
    P4EST_FREE (n_quadrants);
    p4est_destroy (p4est);
    *p4est_supercomm = NULL;
    if (ranks_subcomm) {
      *ranks_subcomm = NULL;
    }

    /* return that p4est does not exist on this rank */
    return 0;
  }

  /* update parallel environment */
  mpiret = sc_MPI_Comm_size (submpicomm, &submpisize);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (submpicomm, &submpirank);
  SC_CHECK_MPI (mpiret);

  if (submpirank == 0) {
    P4EST_VERBOSEF ("%s: Reduce MPI communicator from %i to %i\n",
                    this_fn_name, mpisize, submpisize);
    /* TODO: There should be a function for printing to stdout that works with
     *       sub-communicators. */
  }

  /* translate MPI ranks */
  ranks = P4EST_ALLOC (int, submpisize);
  subranks = P4EST_ALLOC (int, submpisize);
  for (p = 0; p < submpisize; p++) {
    subranks[p] = p;
  }
  mpiret = sc_MPI_Comm_group (submpicomm, &subgroup);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_group (mpicomm, &group);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Group_translate_ranks (subgroup, submpisize, subranks,
                                         group, ranks);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Group_free (&subgroup);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Group_free (&group);
  SC_CHECK_MPI (mpiret);
  P4EST_FREE (subranks);

  /* allocate and set global quadrant count */
  P4EST_FREE (p4est->global_first_quadrant);
  p4est->global_first_quadrant = P4EST_ALLOC (p4est_gloidx_t, submpisize + 1);
  p4est->global_first_quadrant[0] = 0;
  for (p = 0; p < submpisize; p++) {
    P4EST_ASSERT (ranks[p] != sc_MPI_UNDEFINED);
    P4EST_ASSERT (group_add != sc_MPI_GROUP_NULL
                  || 0 < n_quadrants[ranks[p]]);
    p4est->global_first_quadrant[p + 1] =
      p4est->global_first_quadrant[p] + n_quadrants[ranks[p]];
  }
  P4EST_ASSERT (p4est->global_first_quadrant[submpisize] =
                p4est->global_num_quadrants);
  P4EST_FREE (n_quadrants);

  /* set new parallel environment */
  p4est_comm_parallel_env_release (p4est);
  p4est_comm_parallel_env_assign (p4est, submpicomm);
  p4est_comm_parallel_env_duplicate (p4est);
  mpiret = sc_MPI_Comm_free (&submpicomm);
  SC_CHECK_MPI (mpiret);
  P4EST_ASSERT (p4est->mpisize == submpisize);

  /* allocate and set partition information */
  p4est->global_first_position =
    P4EST_ALLOC (p4est_quadrant_t, submpisize + 1);
  if (group_add != sc_MPI_GROUP_NULL) { /* if communication is required */
    p4est_comm_global_partition (p4est, NULL);
  }
  else {                        /* if we can set partition information communication-free */
    for (p = 0; p < submpisize; p++) {
      P4EST_ASSERT (0 == p || ranks[p - 1] < ranks[p]);
      p4est->global_first_position[p] = global_first_position[ranks[p]];
    }
    p4est->global_first_position[submpisize] = global_first_position[mpisize];
  }
  P4EST_FREE (global_first_position);
  if (ranks_subcomm) {
    *ranks_subcomm = ranks;
  }
  else {
    P4EST_FREE (ranks);
  }

  /* check for valid p4est */
  P4EST_ASSERT (p4est_is_valid (p4est));

  /* return that p4est exists on this rank */
  return 1;
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
p4est_comm_global_first_quadrant (p4est_gloidx_t global_num_quadrants,
                                  int mpisize, p4est_gloidx_t * gfq)
{
  int                 i;

  P4EST_ASSERT (gfq != NULL);
  P4EST_ASSERT (mpisize >= 1);
  P4EST_ASSERT (global_num_quadrants >= 0);

  gfq[0] = 0;
  for (i = 1; i < mpisize; ++i) {
    gfq[i] = p4est_partition_cut_gloidx (global_num_quadrants, i, mpisize);
  }
  gfq[mpisize] = global_num_quadrants;
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
p4est_comm_is_empty (p4est_t *p4est, int p)
{
  P4EST_ASSERT (p4est != NULL);
  return p4est_comm_is_empty_gfq (p4est->global_first_quadrant,
                                  p4est->mpisize, p);
}

int
p4est_comm_is_empty_gfq (const p4est_gloidx_t *gfq, int num_procs, int p)
{
  P4EST_ASSERT (gfq != NULL);
  P4EST_ASSERT (0 <= p && p < num_procs);

  return gfq[p] == gfq[p + 1];
}

int
p4est_comm_is_empty_gfp (const p4est_quadrant_t *gfp, int num_procs, int p)
{
  P4EST_ASSERT (gfp != NULL);
  P4EST_ASSERT (0 <= p && p < num_procs);

  return p4est_quadrant_is_equal_piggy (&gfp[p], &gfp[p + 1]);
}

int
p4est_comm_is_contained (p4est_t * p4est, p4est_locidx_t which_tree,
                         const p4est_quadrant_t * q, int rank)
{
  p4est_topidx_t      ctree;
  p4est_quadrant_t    qlast;
  const p4est_quadrant_t *cur;

  P4EST_ASSERT (p4est != NULL && p4est->connectivity != NULL);
  P4EST_ASSERT (p4est->global_first_position != NULL);
  P4EST_ASSERT (0 <= which_tree &&
                which_tree < p4est->connectivity->num_trees);
  P4EST_ASSERT (q != NULL);
  P4EST_ASSERT (0 <= rank && rank < p4est->mpisize);
  P4EST_ASSERT (p4est_quadrant_is_node (q, 1) || p4est_quadrant_is_valid (q));

  /* check whether q begins on a lower processor than rank */
  cur = &p4est->global_first_position[rank];
  P4EST_ASSERT (cur->level == P4EST_QMAXLEVEL);
  ctree = cur->p.which_tree;
  if (which_tree < ctree ||
      (which_tree == ctree &&
       (p4est_quadrant_compare (q, cur) < 0 &&
        (q->x != cur->x || q->y != cur->y
#ifdef P4_TO_P8
         || q->z != cur->z
#endif
        )))) {
    return 0;
  }

  /* check whether q ends on a higher processor than rank */
  ++cur;
  P4EST_ASSERT (cur == &p4est->global_first_position[rank + 1]);
  P4EST_ASSERT (cur->level == P4EST_QMAXLEVEL);
  ctree = cur->p.which_tree;
  if (which_tree > ctree ||
      (which_tree == ctree &&
       (p4est_quadrant_last_descendant (q, &qlast, P4EST_QMAXLEVEL),
        p4est_quadrant_compare (cur, &qlast) <= 0))) {
    return 0;
  }

  /* the quadrant lies fully in the ownership region of rank */
  return 1;
}

int
p4est_comm_is_owner (p4est_t *p4est, p4est_locidx_t which_tree,
                     const p4est_quadrant_t *q, int rank)
{
  P4EST_ASSERT (p4est != NULL);
  P4EST_ASSERT (p4est->connectivity != NULL);

  return p4est_comm_is_owner_gfp
    (p4est->global_first_position, p4est->mpisize,
     p4est->connectivity->num_trees, which_tree, q, rank);
}

int                 p4est_comm_is_owner_gfp
  (const p4est_quadrant_t *gfp, int num_procs,
   p4est_topidx_t num_trees, p4est_locidx_t which_tree,
   const p4est_quadrant_t *q, int rank)
{
  p4est_topidx_t      ctree;
  const p4est_quadrant_t *cur;

  P4EST_ASSERT (gfp != NULL);
  P4EST_ASSERT (0 <= which_tree && which_tree < num_trees);
  P4EST_ASSERT (q != NULL);
  P4EST_ASSERT (0 <= rank && rank < num_procs);
  P4EST_ASSERT (p4est_quadrant_is_node (q, 1) || p4est_quadrant_is_valid (q));

  /* check whether q begins on a lower processor than rank */
  cur = &gfp[rank];
  P4EST_ASSERT (cur->level == P4EST_QMAXLEVEL);
  ctree = cur->p.which_tree;
  if (which_tree < ctree ||
      (which_tree == ctree &&
       (p4est_quadrant_compare (q, cur) < 0 &&
        (q->x != cur->x || q->y != cur->y
#ifdef P4_TO_P8
         || q->z != cur->z
#endif
        )))) {
    return 0;
  }

  /* check whether q lies fully on a higher processor than rank */
  ++cur;
  P4EST_ASSERT (cur == &gfp[rank + 1]);
  P4EST_ASSERT (cur->level == P4EST_QMAXLEVEL);
  ctree = cur->p.which_tree;
  if (which_tree > ctree ||
      (which_tree == ctree &&
       (p4est_quadrant_compare (cur, q) <= 0 ||
        (q->x == cur->x && q->y == cur->y
#ifdef P4_TO_P8
         && q->z == cur->z
#endif
        )))) {
    return 0;
  }

  /* we have not covered the case that q may end on a higher process */
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
  mpiret = sc_MPI_Gather (send, 2, sc_MPI_LONG_LONG_INT,
                          gather, 2, sc_MPI_LONG_LONG_INT, 0, p4est->mpicomm);
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

void
p4est_transfer_fixed (const p4est_gloidx_t * dest_gfq,
                      const p4est_gloidx_t * src_gfq,
                      sc_MPI_Comm mpicomm, int tag,
                      void *dest_data, const void *src_data, size_t data_size)
{
  p4est_transfer_context_t *tc;

  tc = p4est_transfer_fixed_begin (dest_gfq, src_gfq, mpicomm, tag,
                                   dest_data, src_data, data_size);
  p4est_transfer_fixed_end (tc);
}

static void
p4est_transfer_assign_comm (const p4est_gloidx_t * dest_gfq,
                            const p4est_gloidx_t * src_gfq,
                            sc_MPI_Comm mpicomm, int *mpisize, int *mpirank)
{
  int                 mpiret;

  P4EST_ASSERT (dest_gfq != NULL && src_gfq != NULL);
  P4EST_ASSERT (dest_gfq[0] == 0 && src_gfq[0] == 0);
  P4EST_ASSERT (mpicomm != sc_MPI_COMM_NULL);
  P4EST_ASSERT (mpisize != NULL && mpirank != NULL);

  mpiret = sc_MPI_Comm_size (mpicomm, mpisize);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (mpicomm, mpirank);
  SC_CHECK_MPI (mpiret);

  P4EST_ASSERT (dest_gfq[*mpisize] == src_gfq[*mpisize]);
  P4EST_ASSERT (0 <= dest_gfq[*mpirank] &&
                dest_gfq[*mpirank] <= dest_gfq[*mpirank + 1] &&
                dest_gfq[*mpirank + 1] <= dest_gfq[*mpisize]);
  P4EST_ASSERT (0 <= src_gfq[*mpirank] &&
                src_gfq[*mpirank] <= src_gfq[*mpirank + 1] &&
                src_gfq[*mpirank + 1] <= src_gfq[*mpisize]);
}

p4est_transfer_context_t *
p4est_transfer_fixed_begin (const p4est_gloidx_t * dest_gfq,
                            const p4est_gloidx_t * src_gfq,
                            sc_MPI_Comm mpicomm, int tag, void *dest_data,
                            const void *src_data, size_t data_size)
{
  p4est_transfer_context_t *tc;
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
  tc->variable = 0;

  /* there is nothing to do when there is no data */
  if (data_size == 0) {
    return tc;
  }

  /* grab local partition information */
  p4est_transfer_assign_comm (dest_gfq, src_gfq, mpicomm, &mpisize, &mpirank);
  dest_begin = dest_gfq[mpirank];
  dest_end = dest_gfq[mpirank + 1];
  src_begin = src_gfq[mpirank];
  src_end = src_gfq[mpirank + 1];

  /* prepare data copy for local overlap */
  dest_cp = src_cp = NULL;
  cp_len = 0;

  /* figure out subset of processes to receive from */
  if (dest_begin < dest_end) {
    P4EST_ASSERT (dest_data != NULL);

    /* our process as the receiver is not empty */
    first_sender = p4est_bsearch_partition (dest_begin, src_gfq, mpisize);
    P4EST_ASSERT (0 <= first_sender && first_sender < mpisize);
    last_sender =
      p4est_bsearch_partition (dest_end - 1, &src_gfq[first_sender],
                               mpisize - first_sender) + first_sender;
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

      /* choose how to treat the sender process */
      if (gbegin == gend) {
        /* the sender process is empty; we need no message */
        P4EST_ASSERT (first_sender < q && q < last_sender);
        *rq++ = sc_MPI_REQUEST_NULL;
      }
      else {
        /* nonzero message from this sender */
        byte_len = (size_t) (gend - gbegin) * data_size;
        if (q == mpirank) {
          /* on the same rank we remember pointers for memcpy */
          cp_len = byte_len;
          dest_cp = rb;
          *rq++ = sc_MPI_REQUEST_NULL;
        }
        else {
          /* we receive a proper message */
          mpiret = sc_MPI_Irecv (rb, byte_len, sc_MPI_BYTE, q,
                                 tag, mpicomm, rq++);
          SC_CHECK_MPI (mpiret);
        }
        rb += byte_len;
      }
    }
    P4EST_ASSERT (rb - (char *) dest_data ==
                  (ptrdiff_t) ((dest_end - dest_begin) * data_size));
  }

  /* figure out subset of processes to send to */
  if (src_begin < src_end) {
    P4EST_ASSERT (src_data != NULL);

    /* our process as the sender is not empty */
    first_receiver = p4est_bsearch_partition (src_begin, dest_gfq, mpisize);
    P4EST_ASSERT (0 <= first_receiver && first_receiver < mpisize);
    last_receiver =
      p4est_bsearch_partition (src_end - 1, &dest_gfq[first_receiver],
                               mpisize - first_receiver) + first_receiver;
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

      /* choose how to treat the receiver process */
      if (gbegin == gend) {
        /* the receiver process is empty; we need no message */
        P4EST_ASSERT (first_receiver < q && q < last_receiver);
        *rq++ = sc_MPI_REQUEST_NULL;
      }
      else {
        /* nonzero message for this receiver */
        byte_len = (size_t) (gend - gbegin) * data_size;
        if (q == mpirank) {
          /* on the same rank we remember pointers for memcpy */
          P4EST_ASSERT (cp_len == byte_len);
          src_cp = rb;
          *rq++ = sc_MPI_REQUEST_NULL;
        }
        else {
          /* we send a proper message */
          mpiret = sc_MPI_Isend (rb, byte_len, sc_MPI_BYTE, q,
                                 tag, mpicomm, rq++);
          SC_CHECK_MPI (mpiret);
        }
        rb += byte_len;
      }
    }
    P4EST_ASSERT (rb - (char *) src_data ==
                  (ptrdiff_t) ((src_end - src_begin) * data_size));
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

void
p4est_transfer_end (p4est_transfer_context_t * tc)
{
  int                 mpiret;

  P4EST_ASSERT (tc != NULL);

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

  /* the context must disappear */
  P4EST_FREE (tc);
}

void
p4est_transfer_fixed_end (p4est_transfer_context_t * tc)
{
  P4EST_ASSERT (tc != NULL);
  P4EST_ASSERT (tc->variable == 0);

  p4est_transfer_end (tc);
}

void
p4est_transfer_custom (const p4est_gloidx_t * dest_gfq,
                       const p4est_gloidx_t * src_gfq,
                       sc_MPI_Comm mpicomm, int tag,
                       void *dest_data, const int *dest_sizes,
                       const void *src_data, const int *src_sizes)
{
  p4est_transfer_context_t *tc;

  tc = p4est_transfer_custom_begin (dest_gfq, src_gfq, mpicomm, tag,
                                    dest_data, dest_sizes,
                                    src_data, src_sizes);
  p4est_transfer_custom_end (tc);
}

static p4est_transfer_context_t *
p4est_transfer_begin (const p4est_gloidx_t * dest_gfq,
                      const p4est_gloidx_t * src_gfq,
                      sc_MPI_Comm mpicomm, int tag,
                      void *dest_data, const int *dest_sizes,
                      const void *src_data, const int *src_sizes,
                      size_t item_size, int variable)
{
  p4est_transfer_context_t *tc;
  int                 mpiret;
  int                 mpisize, mpirank;
  int                 q;
  int                 first_sender, last_sender;
#ifdef P4EST_ENABLE_DEBUG
  int                 old_last_sender, old_last_receiver;
#endif
  int                 first_receiver, last_receiver;
  int                 i, ilen;
  const int          *rs;
  char               *rb;
  char               *dest_cp, *src_cp;
  size_t              byte_len, cp_len;
  p4est_gloidx_t      dest_begin, dest_end;
  p4est_gloidx_t      src_begin, src_end;
  p4est_gloidx_t      gbegin, gend;
  sc_MPI_Request     *rq;

  /* consistency of internal helper function */
  P4EST_ASSERT (variable == 1 || variable == 2);
  P4EST_ASSERT ((variable == 1 && item_size == 1) || variable == 2);

  /* setup context structure */
  tc = P4EST_ALLOC_ZERO (p4est_transfer_context_t, 1);
  tc->variable = variable;

  /* there is nothing to do when there is no data */
  if (item_size == 0) {
    return tc;
  }

  /* grab local partition information */
  p4est_transfer_assign_comm (dest_gfq, src_gfq, mpicomm, &mpisize, &mpirank);
  dest_begin = dest_gfq[mpirank];
  dest_end = dest_gfq[mpirank + 1];
  src_begin = src_gfq[mpirank];
  src_end = src_gfq[mpirank + 1];

  /* prepare data copy for local overlap */
  dest_cp = src_cp = NULL;
  cp_len = 0;

  /* figure out subset of processes to receive from */
  if (dest_begin < dest_end) {
    P4EST_ASSERT (dest_sizes != NULL);

    /* our process as the receiver is not empty */
    first_sender = p4est_bsearch_partition (dest_begin, src_gfq, mpisize);
    P4EST_ASSERT (0 <= first_sender && first_sender < mpisize);
    last_sender =
      p4est_bsearch_partition (dest_end - 1, &src_gfq[first_sender],
                               mpisize - first_sender) + first_sender;
#ifdef P4EST_ENABLE_DEBUG
    old_last_sender =
      p4est_bsearch_partition (dest_end - 1, src_gfq, mpisize);
    P4EST_ASSERT (last_sender == old_last_sender);
#endif
    P4EST_ASSERT (first_sender <= last_sender && last_sender < mpisize);
    tc->num_senders = last_sender - first_sender + 1;
    P4EST_ASSERT (tc->num_senders > 0);

    /* go through sender processes and post receive calls */
    gend = dest_begin;
    rq = tc->recv_req = P4EST_ALLOC (sc_MPI_Request, tc->num_senders);
    rb = (char *) dest_data;
    rs = dest_sizes;
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

      /* determine message size for this sender */
      byte_len = 0;
      ilen = (int) (gend - gbegin);
      for (i = 0; i < ilen; ++i) {
        byte_len += item_size * *rs++;
      }

      /* choose how to treat the sender process */
      if (byte_len == 0) {
        /* the sender process or the message is empty; we need no message */
        *rq++ = sc_MPI_REQUEST_NULL;
      }
      else {
        P4EST_ASSERT (dest_data != NULL);
        if (q == mpirank) {
          /* on the same rank we remember pointers for memcpy */
          cp_len = byte_len;
          dest_cp = rb;
          *rq++ = sc_MPI_REQUEST_NULL;
        }
        else {
          /* we receive a proper message */
          mpiret = sc_MPI_Irecv (rb, byte_len, sc_MPI_BYTE, q,
                                 tag, mpicomm, rq++);
          SC_CHECK_MPI (mpiret);
        }
        rb += byte_len;
      }
    }
    P4EST_ASSERT (rs - dest_sizes == (ptrdiff_t) (dest_end - dest_begin));
  }

  /* figure out subset of processes to send to */
  if (src_begin < src_end) {
    P4EST_ASSERT (src_sizes != NULL);

    /* our process as the sender is not empty */
    first_receiver = p4est_bsearch_partition (src_begin, dest_gfq, mpisize);
    P4EST_ASSERT (0 <= first_receiver && first_receiver < mpisize);
    last_receiver =
      p4est_bsearch_partition (src_end - 1, &dest_gfq[first_receiver],
                               mpisize - first_receiver) + first_receiver;
#ifdef P4EST_ENABLE_DEBUG
    old_last_receiver =
      p4est_bsearch_partition (src_end - 1, dest_gfq, mpisize);
    P4EST_ASSERT (last_receiver == old_last_receiver);
#endif
    P4EST_ASSERT (first_receiver <= last_receiver && last_receiver < mpisize);
    tc->num_receivers = last_receiver - first_receiver + 1;
    P4EST_ASSERT (tc->num_receivers > 0);

    /* go through receiver processes and post send calls */
    gend = src_begin;
    rq = tc->send_req = P4EST_ALLOC (sc_MPI_Request, tc->num_receivers);
    rb = (char *) src_data;
    rs = src_sizes;
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

      /* determine message size for this receiver */
      byte_len = 0;
      ilen = (int) (gend - gbegin);
      for (i = 0; i < ilen; ++i) {
        byte_len += item_size * *rs++;
      }

      /* choose how to treat the receiver process */
      if (byte_len == 0) {
        /* the receiver process or the message is empty; we need no message */
        *rq++ = sc_MPI_REQUEST_NULL;
      }
      else {
        P4EST_ASSERT (src_data != NULL);
        if (q == mpirank) {
          /* on the same rank we remember pointers for memcpy */
          P4EST_ASSERT (cp_len == byte_len);
          src_cp = rb;
          *rq++ = sc_MPI_REQUEST_NULL;
        }
        else {
          /* we send a proper message */
          mpiret = sc_MPI_Isend (rb, byte_len, sc_MPI_BYTE, q,
                                 tag, mpicomm, rq++);
          SC_CHECK_MPI (mpiret);
        }
        rb += byte_len;
      }
    }
    P4EST_ASSERT (rs - src_sizes == (ptrdiff_t) (src_end - src_begin));
  }

  /* copy the data that remains local */
  P4EST_ASSERT ((dest_cp == NULL) == (src_cp == NULL));
  if (cp_len > 0) {
    P4EST_ASSERT (dest_cp != NULL && src_cp != NULL);
    memcpy (dest_cp, src_cp, cp_len);
  }

  /* the rest goes into the p4est_transfer_custom_end function */
  return tc;
}

p4est_transfer_context_t *
p4est_transfer_custom_begin (const p4est_gloidx_t * dest_gfq,
                             const p4est_gloidx_t * src_gfq,
                             sc_MPI_Comm mpicomm, int tag,
                             void *dest_data, const int *dest_sizes,
                             const void *src_data, const int *src_sizes)
{
  return p4est_transfer_begin
    (dest_gfq, src_gfq, mpicomm, tag,
     dest_data, dest_sizes, src_data, src_sizes, 1, 1);
}

p4est_transfer_context_t *
p4est_transfer_items_begin (const p4est_gloidx_t * dest_gfq,
                            const p4est_gloidx_t * src_gfq,
                            sc_MPI_Comm mpicomm, int tag,
                            void *dest_data, const int *dest_counts,
                            const void *src_data, const int *src_counts,
                            size_t item_size)
{
  return p4est_transfer_begin
    (dest_gfq, src_gfq, mpicomm, tag,
     dest_data, dest_counts, src_data, src_counts, item_size, 2);
}

void
p4est_transfer_custom_end (p4est_transfer_context_t * tc)
{
  P4EST_ASSERT (tc != NULL);
  P4EST_ASSERT (tc->variable == 1);

  p4est_transfer_end (tc);
}

void
p4est_transfer_items (const p4est_gloidx_t * dest_gfq,
                      const p4est_gloidx_t * src_gfq,
                      sc_MPI_Comm mpicomm, int tag,
                      void *dest_data, const int *dest_sizes,
                      const void *src_data, const int *src_sizes,
                      size_t item_size)
{
  p4est_transfer_context_t *tc;

  tc = p4est_transfer_items_begin (dest_gfq, src_gfq, mpicomm, tag,
                                   dest_data, dest_sizes,
                                   src_data, src_sizes, item_size);
  p4est_transfer_items_end (tc);
}

void
p4est_transfer_items_end (p4est_transfer_context_t * tc)
{
  P4EST_ASSERT (tc != NULL);
  P4EST_ASSERT (tc->variable == 2);

  p4est_transfer_end (tc);
}
