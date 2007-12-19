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

#include <p4est.h>
#include <p4est_algorithms.h>
#include <p4est_communication.h>
#include <p4est_base.h>

typedef struct
{
  int8_t              will_receive, have_count, have_load;
  int32_t             expect_count;
  p4est_array_t       send_first, recv_first;
  p4est_array_t       send_second, recv_second;
#ifdef HAVE_MPI
  MPI_Request         request_send_first_count, request_send_first_load;
  MPI_Request         request_send_second[2], request_recv_second[2];
#endif
}
p4est_balance_peer_t;

static const int64_t initial_quadrants_per_processor = 15;
static const int    number_toread_quadrants = 32;

p4est_t            *
p4est_new (MPI_Comm mpicomm, FILE * nout, p4est_connectivity_t * connectivity,
           int data_size, p4est_init_t init_fn)
{
#ifdef HAVE_MPI
  int                 mpiret;
#endif
  int                 i, must_remove_last_quadrant;
  int8_t              level;
  int32_t             j, num_trees;
  int64_t             tree_num_quadrants, global_num_quadrants;
  int64_t             first_tree, first_quadrant, first_tree_quadrant;
  int64_t             last_tree, last_quadrant, last_tree_quadrant;
  int64_t             quadrant_index;
  p4est_t            *p4est;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *quad;
  p4est_quadrant_t    a;
  p4est_quadrant_t    b;

  p4est = P4EST_ALLOC_ZERO (p4est_t, 1);
  P4EST_CHECK_ALLOC (p4est);

  /* assign some data members */
  p4est->nout = nout;
  p4est->data_size = data_size;
  p4est->connectivity = connectivity;
  num_trees = connectivity->num_trees;

  p4est->mpicomm = mpicomm;
  p4est->mpisize = 1;
  p4est->mpirank = 0;
#ifdef HAVE_MPI
  if (p4est->mpicomm != MPI_COMM_NULL) {
    mpiret = MPI_Comm_size (p4est->mpicomm, &p4est->mpisize);
    P4EST_CHECK_MPI (mpiret);

    mpiret = MPI_Comm_rank (p4est->mpicomm, &p4est->mpirank);
    P4EST_CHECK_MPI (mpiret);
  }
#endif

  /* allocate memory pools */
  if (p4est->data_size > 0) {
    p4est->user_data_pool = p4est_mempool_new (p4est->data_size);
  }
  else {
    p4est->user_data_pool = NULL;
  }
  p4est->quadrant_pool = p4est_mempool_new (sizeof (p4est_quadrant_t));

  /* determine uniform level of initial tree */
  tree_num_quadrants = 1;
  for (level = 0; level < 16; ++level) {
    if (tree_num_quadrants >=
        (p4est->mpisize * initial_quadrants_per_processor) / num_trees) {
      break;
    }
    tree_num_quadrants *= 4;
  }
  P4EST_ASSERT (level < 16 && tree_num_quadrants <= INT32_MAX);

  /* compute index of first tree for this processor */
  global_num_quadrants = tree_num_quadrants * num_trees;
  first_quadrant = (global_num_quadrants * p4est->mpirank) / p4est->mpisize;
  first_tree = first_quadrant / tree_num_quadrants;
  first_tree_quadrant = first_quadrant - first_tree * tree_num_quadrants;
  last_quadrant =
    (global_num_quadrants * (p4est->mpirank + 1)) / p4est->mpisize - 1;
  last_tree = last_quadrant / tree_num_quadrants;
  last_tree_quadrant = last_quadrant - last_tree * tree_num_quadrants;

  /* check ranges of various integers to be 32bit compatible */
  P4EST_ASSERT (first_tree <= last_tree && last_tree < num_trees);
  P4EST_ASSERT (0 <= first_tree_quadrant && 0 <= last_tree_quadrant);
  P4EST_ASSERT (first_tree_quadrant < tree_num_quadrants);
  P4EST_ASSERT (last_tree_quadrant < tree_num_quadrants);
  if (first_tree == last_tree) {
    P4EST_ASSERT (first_tree_quadrant < last_tree_quadrant);
  }

  /* print some diagnostics */
  if (p4est->nout != NULL) {
    if (p4est->mpirank == 0) {
      fprintf (p4est->nout, "New forest: %d trees %d processors\n",
               num_trees, p4est->mpisize);
      fprintf (p4est->nout,
               "   initial level %d potential global quadrants %lld per tree %lld\n",
               level, (long long) global_num_quadrants,
               (long long) tree_num_quadrants);
    }
    fprintf (p4est->nout, "[%d] first tree %lld first quadrant %lld"
             " global quadrant %lld\n", p4est->mpirank,
             (long long) first_tree, (long long) first_tree_quadrant,
             (long long) first_quadrant);
    fprintf (p4est->nout, "[%d] last tree %lld last quadrant %lld"
             " global quadrant %lld\n", p4est->mpirank,
             (long long) last_tree, (long long) last_tree_quadrant,
             (long long) last_quadrant);
  }

  /* allocate trees and quadrants */
  p4est->trees = p4est_array_new (sizeof (p4est_tree_t));
  p4est_array_resize (p4est->trees, num_trees);
  for (j = 0; j < num_trees; ++j) {
    tree = p4est_array_index (p4est->trees, j);
    tree->quadrants = p4est_array_new (sizeof (p4est_quadrant_t));
    for (i = 0; i <= P4EST_MAXLEVEL; ++i) {
      tree->quadrants_per_level[i] = 0;
    }
    tree->maxlevel = 0;
  }
  p4est->local_num_quadrants = 0;
  p4est->global_num_quadrants = 0;

  /* for every locally non-empty tree fill first and last quadrant */
  for (j = first_tree; j <= last_tree; ++j) {
    tree = p4est_array_index (p4est->trees, j);
    must_remove_last_quadrant = 0;

    /* set morton id of first quadrant and initialize user data */
    if (j == first_tree) {
      p4est_quadrant_set_morton (&a, level, first_tree_quadrant);
    }
    else {
      p4est_quadrant_set_morton (&a, level, 0);
    }
    if (p4est->nout != NULL) {
      fprintf (p4est->nout, "[%d] tree %d first morton 0x%x 0x%x\n",
               p4est->mpirank, j, a.x, a.y);
    }
    p4est_quadrant_init_data (p4est, j, &a, init_fn);

    /* set morton id of last quadrant */
    if (tree_num_quadrants == 1 ||
        (j == first_tree && first_tree_quadrant == tree_num_quadrants - 1)) {
      /* There is only a in the tree */
      p4est_array_resize (tree->quadrants, 1);
      quad = p4est_array_index (tree->quadrants, 0);
      *quad = a;
      tree->maxlevel = a.level;
      ++tree->quadrants_per_level[a.level];
    }
    else {
      if (j == last_tree) {
        if (last_tree_quadrant == tree_num_quadrants - 1) {
          quadrant_index = last_tree_quadrant;
        }
        else {
          quadrant_index = last_tree_quadrant + 1;
          must_remove_last_quadrant = 1;
        }
        p4est_quadrant_set_morton (&b, level, quadrant_index);
      }
      else {
        p4est_quadrant_set_morton (&b, level, tree_num_quadrants - 1);
      }
      if (p4est->nout != NULL) {
        fprintf (p4est->nout, "[%d] tree %d last morton 0x%x 0x%x\n",
                 p4est->mpirank, j, b.x, b.y);
      }
      if (!must_remove_last_quadrant) {
        p4est_quadrant_init_data (p4est, j, &b, init_fn);
      }

      /* now run algorithm CompleteRegion (&tree->quadrants) here */
      p4est_complete_region (p4est, &a, 1, &b, !must_remove_last_quadrant,
                             tree, j, init_fn);
    }
    if (p4est->nout != NULL) {
      fprintf (p4est->nout, "[%d] tree %d quadrants %d\n",
               p4est->mpirank, j, tree->quadrants->elem_count);
    }
    p4est->local_num_quadrants += tree->quadrants->elem_count;
  }

  /* compute some member variables */
  p4est->first_local_tree = first_tree;
  p4est->last_local_tree = last_tree;
  p4est->local_num_trees = last_tree - first_tree + 1;
  p4est_comm_count_quadrants (p4est);

  /* compute global partition information */
  p4est->global_first_indices = P4EST_ALLOC (int32_t,
                                             3 * (p4est->mpisize + 1));
  P4EST_CHECK_ALLOC (p4est->global_first_indices);
  p4est_comm_global_partition (p4est);

  /* print more statistics */
  if (p4est->nout != NULL) {
    fprintf (p4est->nout, "[%d] total local quadrants %d\n",
             p4est->mpirank, p4est->local_num_quadrants);
    if (p4est->mpirank == 0) {
      fprintf (p4est->nout, "   total global quadrants %lld\n",
               (long long int) p4est->global_num_quadrants);
    }
  }

  P4EST_ASSERT (p4est_is_valid (p4est));

  return p4est;
}

void
p4est_destroy (p4est_t * p4est)
{
  int32_t             j;
  p4est_tree_t       *tree;

  for (j = 0; j < p4est->connectivity->num_trees; ++j) {
    tree = p4est_array_index (p4est->trees, j);
    p4est_array_destroy (tree->quadrants);
  }
  p4est_array_destroy (p4est->trees);

  if (p4est->user_data_pool != NULL) {
    p4est_mempool_destroy (p4est->user_data_pool);
  }
  p4est_mempool_destroy (p4est->quadrant_pool);

  P4EST_FREE (p4est->global_first_indices);

  P4EST_FREE (p4est);
}

void
p4est_refine (p4est_t * p4est, p4est_refine_t refine_fn, p4est_init_t init_fn)
{
  int                 quadrant_pool_size, data_pool_size;
  int                 dorefine;
  int8_t              i, maxlevel;
  int32_t             j, movecount;
  int32_t             current, restpos, incount;
  p4est_list_t       *list;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *q, *qalloc, *qpop;
  p4est_quadrant_t   *c0, *c1, *c2, *c3;
  int                *key;

  P4EST_ASSERT (p4est_is_valid (p4est));

  /*
     q points to a quadrant that is an array member
     qalloc is a quadrant that has been allocated through quadrant_pool
     qpop is a quadrant that has been allocated through quadrant_pool
     never mix these two types of quadrant pointers
   */
  key = &dorefine;              /* use this to create a unique user_data pointer */
  list = p4est_list_new (NULL);
  p4est->local_num_quadrants = 0;

  /* loop over all local trees */
  for (j = p4est->first_local_tree; j <= p4est->last_local_tree; ++j) {
    tree = p4est_array_index (p4est->trees, j);
    quadrant_pool_size = p4est->quadrant_pool->elem_count;
    if (p4est->user_data_pool != NULL) {
      data_pool_size = p4est->user_data_pool->elem_count;
    }

    /* initial log message for this tree */
    if (p4est->nout != NULL) {
      fprintf (p4est->nout, "[%d] Into refine tree %d with %d\n",
               p4est->mpirank, j, tree->quadrants->elem_count);
    }

    /* reset the quadrant counters */
    maxlevel = 0;
    for (i = 0; i <= P4EST_MAXLEVEL; ++i) {
      tree->quadrants_per_level[i] = 0;
    }

    /* run through the array to find first quadrant to be refined */
    q = NULL;
    dorefine = 0;
    incount = tree->quadrants->elem_count;
    for (current = 0; current < incount; ++current) {
      q = p4est_array_index (tree->quadrants, current);
      dorefine = ((q->level < P4EST_MAXLEVEL) && refine_fn (p4est, j, q));
      if (dorefine) {
        break;
      }
      maxlevel = (int8_t) P4EST_MAX (maxlevel, q->level);
      ++tree->quadrants_per_level[q->level];
    }
    if (!dorefine) {
      p4est->local_num_quadrants += incount;
      continue;
    }

    /* now we have a quadrant to refine, prepend it to the list */
    qalloc = p4est_mempool_alloc (p4est->quadrant_pool);
    *qalloc = *q;               /* never prepend array members */
    p4est_list_prepend (list, qalloc);  /* only newly allocated quadrants */

    /*
       current points to the next array member to write
       restpos points to the next array member to read
     */
    restpos = current + 1;

    /* run through the list and refine recursively */
    while (list->elem_count > 0) {
      qpop = p4est_list_pop (list);
      if (dorefine ||
          ((qpop->level < P4EST_MAXLEVEL) && refine_fn (p4est, j, qpop))) {
        dorefine = 0;           /* a marker so that refine_fn is never called twice */
        p4est_array_resize (tree->quadrants, tree->quadrants->elem_count + 3);

        /* compute children and prepend them to the list */
        if (qpop->user_data != key) {
          p4est_quadrant_free_data (p4est, qpop);
        }
        c0 = qpop;
        c1 = p4est_mempool_alloc (p4est->quadrant_pool);
        c2 = p4est_mempool_alloc (p4est->quadrant_pool);
        c3 = p4est_mempool_alloc (p4est->quadrant_pool);
        p4est_quadrant_children (qpop, c0, c1, c2, c3);
        c0->user_data = key;
        c1->user_data = key;
        c2->user_data = key;
        c3->user_data = key;
        p4est_list_prepend (list, c3);
        p4est_list_prepend (list, c2);
        p4est_list_prepend (list, c1);
        p4est_list_prepend (list, c0);
      }
      else {
        /* need to make room in the array to store this new quadrant */
        if (restpos < incount && current == restpos) {
          movecount = P4EST_MIN (incount - restpos, number_toread_quadrants);
          while (movecount > 0) {
            q = p4est_array_index (tree->quadrants, restpos);
            qalloc = p4est_mempool_alloc (p4est->quadrant_pool);
            *qalloc = *q;       /* never append array members */
            p4est_list_append (list, qalloc);   /* only newly allocated quadrants */
            --movecount;
            ++restpos;
          }
        }

        /* store new quadrant and update counters */
        if (qpop->user_data == key) {
          p4est_quadrant_init_data (p4est, j, qpop, init_fn);
        }
        q = p4est_array_index (tree->quadrants, current);
        *q = *qpop;
        maxlevel = (int8_t) P4EST_MAX (maxlevel, qpop->level);
        ++tree->quadrants_per_level[qpop->level];
        ++current;
        p4est_mempool_free (p4est->quadrant_pool, qpop);
      }
    }
    tree->maxlevel = maxlevel;
    p4est->local_num_quadrants += tree->quadrants->elem_count;

    P4EST_ASSERT (restpos == incount);
    P4EST_ASSERT (current == tree->quadrants->elem_count);
    P4EST_ASSERT (list->first == NULL && list->last == NULL);
    P4EST_ASSERT (quadrant_pool_size == p4est->quadrant_pool->elem_count);
    if (p4est->user_data_pool != NULL) {
      P4EST_ASSERT (data_pool_size + tree->quadrants->elem_count ==
                    p4est->user_data_pool->elem_count + incount);
    }
    P4EST_ASSERT (p4est_tree_is_sorted (tree));
    P4EST_ASSERT (p4est_tree_is_complete (tree));

    /* final log message for this tree */
    if (p4est->nout != NULL) {
      fprintf (p4est->nout, "[%d] Done refine tree %d now %d\n",
               p4est->mpirank, j, tree->quadrants->elem_count);
    }
  }

  p4est_list_destroy (list);

  /* compute global number of quadrants */
  p4est_comm_count_quadrants (p4est);

  P4EST_ASSERT (p4est_is_valid (p4est));
}

void
p4est_coarsen (p4est_t * p4est, p4est_coarsen_t coarsen_fn,
               p4est_init_t init_fn)
{
  int                 k, couldbegood, data_pool_size;
  int                 incount, removed, num_quadrants;
  int                 first, last, rest, before;
  int8_t              i, maxlevel;
  int32_t             j;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *c[4];
  p4est_quadrant_t   *cfirst, *clast;

  P4EST_ASSERT (p4est_is_valid (p4est));

  /* loop over all local trees */
  for (j = p4est->first_local_tree; j <= p4est->last_local_tree; ++j) {
    tree = p4est_array_index (p4est->trees, j);
    if (p4est->user_data_pool != NULL) {
      data_pool_size = p4est->user_data_pool->elem_count;
    }
    removed = 0;

    /* initial log message for this tree */
    if (p4est->nout != NULL) {
      fprintf (p4est->nout, "[%d] Into coarsen tree %d with %d\n",
               p4est->mpirank, j, tree->quadrants->elem_count);
    }

    /* Initialize array indices.
       If children are coarsened, the array will have an empty window.
       first   index of the first child to be considered
       last    index of the last child before the hole in the array
       before  number of children before the hole in the array
       rest    index of the first child after the hole in the array
     */
    first = last = 0;
    before = rest = 1;

    /* run through the array and coarsen recursively */
    incount = tree->quadrants->elem_count;
    while (rest + 3 - before < incount) {
      couldbegood = 1;
      for (k = 0; k < 4; ++k) {
        if (k < before) {
          c[k] = p4est_array_index (tree->quadrants, first + k);
          if (k != p4est_quadrant_child_id (c[k])) {
            couldbegood = 0;
            break;
          }
        }
        else {
          c[k] = p4est_array_index (tree->quadrants, rest + k - before);
        }
      }
      if (couldbegood &&
          p4est_quadrant_is_family (c[0], c[1], c[2], c[3]) &&
          coarsen_fn (p4est, j, c[0], c[1], c[2], c[3])) {
        /* coarsen now */
        for (k = 0; k < 4; ++k) {
          p4est_quadrant_free_data (p4est, c[k]);
        }
        tree->quadrants_per_level[c[0]->level] -= 4;
        cfirst = c[0];
        p4est_quadrant_parent (c[0], cfirst);
        p4est_quadrant_init_data (p4est, j, cfirst, init_fn);
        tree->quadrants_per_level[cfirst->level] += 1;
        p4est->local_num_quadrants -= 3;
        removed += 3;

        rest += 4 - before;
        last = first;
        first -= p4est_quadrant_child_id (cfirst);
        first = P4EST_MAX (first, 0);
      }
      else {
        /* do nothing, just move the counters and the hole */
        ++first;
        if (first > last) {
          if (first != rest) {
            cfirst = p4est_array_index (tree->quadrants, first);
            clast = p4est_array_index (tree->quadrants, rest);
            *cfirst = *clast;
          }
          last = first;
          ++rest;
        }
      }
      before = last - first + 1;
    }

    /* adjust final array size */
    first = last;
    if (first + 1 < rest) {
      while (rest < incount) {
        ++first;
        cfirst = p4est_array_index (tree->quadrants, first);
        clast = p4est_array_index (tree->quadrants, rest);
        *cfirst = *clast;
        ++rest;
      }
      p4est_array_resize (tree->quadrants, first + 1);
    }

    /* compute maximum level */
    maxlevel = 0;
    num_quadrants = 0;
    for (i = 0; i <= P4EST_MAXLEVEL; ++i) {
      P4EST_ASSERT (tree->quadrants_per_level[i] >= 0);
      num_quadrants += tree->quadrants_per_level[i];
      if (tree->quadrants_per_level[i] > 0) {
        maxlevel = i;
      }
    }
    tree->maxlevel = maxlevel;

    /* do some sanity checks */
    P4EST_ASSERT (num_quadrants == tree->quadrants->elem_count);
    P4EST_ASSERT (tree->quadrants->elem_count == incount - removed);
    if (p4est->user_data_pool != NULL) {
      P4EST_ASSERT (data_pool_size - removed ==
                    p4est->user_data_pool->elem_count);
    }
    P4EST_ASSERT (p4est_tree_is_sorted (tree));
    P4EST_ASSERT (p4est_tree_is_complete (tree));

    /* final log message for this tree */
    if (p4est->nout != NULL) {
      fprintf (p4est->nout, "[%d] Done coarsen tree %d now %d\n",
               p4est->mpirank, j, tree->quadrants->elem_count);
    }
  }

  /* compute global number of quadrants */
  p4est_comm_count_quadrants (p4est);

  P4EST_ASSERT (p4est_is_valid (p4est));
}

void
p4est_balance (p4est_t * p4est, p4est_init_t init_fn)
{
  const int           rank = p4est->mpirank;
  int32_t             j;
  int32_t             rh;
  int32_t             qcount, treecount;
  int32_t             first_tree, last_tree, next_tree;
  int32_t             first_peer, last_peer;
  p4est_tree_t       *tree;
  p4est_balance_peer_t *peer;
  p4est_quadrant_t    mylow, nextlow;
#ifdef HAVE_MPI
  int                 mpiret, qbytes;
  int                 first_index, last_index;
  int                 k, l;
  int                 which, comp, scount;
  int                 request_count, outcount;
  int                *wait_indices;
  int32_t             i, qtree;
  int32_t             qh;
  int32_t             owner, first_owner, last_owner;
  int32_t             mypb[2], *peer_boundaries;
  p4est_quadrant_t    ld, ins[9];       /* insulation layer including its center */
  p4est_quadrant_t   *q, *s;
  p4est_array_t      *peers, *qarray;
  MPI_Request        *requests;
  MPI_Status         *statuses;
#endif

  P4EST_ASSERT (p4est_is_valid (p4est));

#ifdef HAVE_MPI
  /* allocate temporary storage */
  peers = p4est_array_new (sizeof (p4est_balance_peer_t));
  p4est_array_resize (peers, p4est->mpisize);
  for (i = 0; i < p4est->mpisize; ++i) {
    peer = p4est_array_index (peers, i);
    p4est_array_init (&peer->send_first, sizeof (p4est_quadrant_t));
    p4est_array_init (&peer->recv_first, sizeof (p4est_quadrant_t));
    p4est_array_init (&peer->send_second, sizeof (p4est_quadrant_t));
    p4est_array_init (&peer->recv_second, sizeof (p4est_quadrant_t));
  }
  /* will contain first and last peer (inclusive) for each processor */
  peer_boundaries = P4EST_ALLOC (int32_t, 2 * p4est->mpisize);
  P4EST_CHECK_ALLOC (peer_boundaries);
  requests = P4EST_ALLOC (MPI_Request, p4est->mpisize);
  P4EST_CHECK_ALLOC (requests);
  statuses = P4EST_ALLOC (MPI_Status, p4est->mpisize);
  P4EST_CHECK_ALLOC (statuses);
  wait_indices = P4EST_ALLOC (int, p4est->mpisize);
  P4EST_CHECK_ALLOC (wait_indices);
#endif

  /* compute first quadrants on finest level for comparison for me and next */
  first_peer = last_peer = rank;
  first_tree = p4est->first_local_tree;
  last_tree = p4est->last_local_tree;
  P4EST_ASSERT (p4est->global_first_indices[3 * rank + 0] == first_tree);
  mylow.x = p4est->global_first_indices[3 * rank + 1];
  mylow.y = p4est->global_first_indices[3 * rank + 2];
  mylow.level = P4EST_MAXLEVEL;
  next_tree = p4est->global_first_indices[3 * (rank + 1) + 0];
  P4EST_ASSERT (next_tree == last_tree || next_tree == last_tree + 1);
  nextlow.x = p4est->global_first_indices[3 * (rank + 1) + 1];
  nextlow.y = p4est->global_first_indices[3 * (rank + 1) + 2];
  nextlow.level = P4EST_MAXLEVEL;
  rh = (1 << P4EST_MAXLEVEL);

  /* loop over all local trees to assemble first send list */
  p4est->local_num_quadrants = 0;
  for (j = first_tree; j <= last_tree; ++j) {
    tree = p4est_array_index (p4est->trees, j);

    /* initial log message for this tree */
    if (p4est->nout != NULL) {
      fprintf (p4est->nout, "[%d] Into balance tree %d with %d\n",
               rank, j, tree->quadrants->elem_count);
    }

    /* local balance first pass */
    p4est_balance_subtree (p4est, tree, j, init_fn);
    p4est->local_num_quadrants += (qcount = tree->quadrants->elem_count);

#ifdef HAVE_MPI
    /* check if this tree is not shared with other processors */
    if (p4est->mpicomm == MPI_COMM_NULL) {
      continue;
    }
    if ((j > first_tree || (mylow.x == 0 && mylow.y == 0)) &&
        (j < last_tree || (nextlow.x == 0 && nextlow.y == 0))) {
      continue;
    }

    /* identify boundary quadrants and prepare them to be sent */
    for (i = 0; i < qcount; ++i) {
      /* this quadrant may be on the boundary with a range of processors */
      first_owner = last_owner = rank;
      q = p4est_array_index (tree->quadrants, i);
      qh = (1 << (P4EST_MAXLEVEL - q->level));
      for (k = 0; k < 3; ++k) {
        for (l = 0; l < 3; ++l) {
          which = k * 3 + l;    /* 0..8 */
          /* exclude myself from the queries */
          if (which == 4) {
            continue;
          }
          s = &ins[which];
          *s = *q;
          s->x += (l - 1) * qh;
          s->y += (k - 1) * qh;
          if ((s->x < 0 || s->x >= rh) || (s->y < 0 || s->y >= rh)) {
            /* this quadrant is relevant for inter-tree balancing */
            continue;
          }
          comp = p4est_quadrant_compare (s, q);
          P4EST_ASSERT (comp != 0);
          if (comp < 0) {
            /* querying s is equivalent to querying first descendent */
            owner = p4est_comm_find_owner (p4est, j, s, rank);
            P4EST_ASSERT (owner <= rank);
            first_owner = P4EST_MIN (owner, first_owner);
          }
          else {
            p4est_quadrant_last_descendent (s, &ld, P4EST_MAXLEVEL);
            owner = p4est_comm_find_owner (p4est, j, &ld, rank);
            P4EST_ASSERT (owner >= rank);
            last_owner = P4EST_MAX (owner, last_owner);
          }
          if (owner != rank) {
            if (p4est->nout != NULL) {
              fprintf (p4est->nout, "[%d] Tree %d 0x%x 0x%x %d owner %d\n",
                       rank, j,
                       ins[which].x, ins[which].y, ins[which].level, owner);
            }
          }
        }
      }
      /*
       * send q to all processors on its insulation layer
       * rely on the space filling curve to have not too many owners
       */
      for (owner = first_owner; owner <= last_owner; ++owner) {
        if (owner == rank) {
          continue;
        }
        peer = p4est_array_index (peers, owner);
        if (p4est_array_bsearch (&peer->send_first, q,
                                 p4est_quadrant_compare) != NULL) {
          continue;
        }
        scount = peer->send_first.elem_count;
        p4est_array_resize (&peer->send_first, scount + 1);
        s = p4est_array_index (&peer->send_first, scount);
        *s = *q;
        s->user_data = (void *) j;      /* piggy back tree id with quadrant */
      }
      first_peer = P4EST_MIN (first_owner, first_peer);
      last_peer = P4EST_MAX (last_owner, last_peer);
    }
#endif /* HAVE_MPI */
  }

#ifdef HAVE_MPI
  /* distribute information about first and last peers to inform receivers */
  if (p4est->mpicomm != MPI_COMM_NULL) {
    mypb[0] = first_peer;
    mypb[1] = last_peer;
    mpiret = MPI_Allgather (mypb, 2, MPI_INT, peer_boundaries, 2, MPI_INT,
                            p4est->mpicomm);
    P4EST_CHECK_MPI (mpiret);
  }

  /*
   * loop over all peers and send first round of quadrants
   * for intra-tree balancing, each load is contained in one tree
   */
  for (j = first_peer; j <= last_peer; ++j) {
    if (j == rank) {
      continue;
    }
    peer = p4est_array_index (peers, j);

    /* first send number of quadrants to be expected */
    qcount = peer->send_first.elem_count;
    if (p4est->nout != NULL) {
      fprintf (p4est->nout, "[%d] Balance A send %d quadrants to %d\n",
               rank, qcount, j);
    }
    mpiret = MPI_Isend (&qcount, 1, MPI_INT,
                        j, P4EST_COMM_BALANCE_FIRST_COUNT,
                        p4est->mpicomm, &peer->request_send_first_count);
    P4EST_CHECK_MPI (mpiret);

    /* then send the actual quadrants */
    if (qcount > 0) {
      qbytes = qcount * sizeof (p4est_quadrant_t);
      mpiret = MPI_Isend (peer->send_first.array, qbytes, MPI_CHAR,
                          j, P4EST_COMM_BALANCE_FIRST_LOAD,
                          p4est->mpicomm, &peer->request_send_first_load);
      P4EST_CHECK_MPI (mpiret);
    }
  }

  /* find out who is sending to me and receive quadrant counts */
  request_count = 0;
  for (j = 0; j < p4est->mpisize; ++j) {
    peer = p4est_array_index (peers, j);
    peer->will_receive = 0;
    peer->have_count = 0;
    peer->have_load = 0;
    peer->expect_count = 0;
    requests[j] = MPI_REQUEST_NULL;
    wait_indices[j] = -1;
    if (j == rank) {
      continue;
    }
    first_peer = peer_boundaries[2 * j + 0];
    last_peer = peer_boundaries[2 * j + 1];
    if (rank < first_peer || rank > last_peer) {
      continue;
    }
    peer = p4est_array_index (peers, j);
    peer->will_receive = 1;
    ++request_count;

    /* processor j is sending to me */
    mpiret = MPI_Irecv (&peer->expect_count, 1, MPI_INT,
                        j, P4EST_COMM_BALANCE_FIRST_COUNT,
                        p4est->mpicomm, &requests[j]);
    P4EST_CHECK_MPI (mpiret);
  }

  /* wait for quadrant counts and post receive and send for quadrants */
  while (request_count > 0) {
    mpiret = MPI_Waitsome (p4est->mpisize, requests,
                           &outcount, wait_indices, statuses);
    P4EST_CHECK_MPI (mpiret);
    P4EST_ASSERT (outcount > 0);
    for (i = 0; i < outcount; ++i) {
      /* retrieve sender's rank */
      j = wait_indices[i];
      wait_indices[i] = -1;
      P4EST_ASSERT (0 <= j && j < p4est->mpisize);

      /* check if we are in receiving count or load */
      peer = p4est_array_index (peers, j);
      P4EST_ASSERT (peer->will_receive);
      P4EST_ASSERT (!peer->have_load);
      if (!peer->have_count) {
        peer->have_count = 1;
        qcount = peer->expect_count;
        if (p4est->nout != NULL) {
          fprintf (p4est->nout,
                   "[%d] Balance A recv %d quadrants from %d\n",
                   rank, qcount, j);
        }
        if (qcount > 0) {
          /* received nonzero count, post receive for load */
          p4est_array_resize (&peer->recv_first, qcount);
          qbytes = qcount * sizeof (p4est_quadrant_t);
          mpiret = MPI_Irecv (peer->recv_first.array, qbytes, MPI_CHAR,
                              j, P4EST_COMM_BALANCE_FIRST_LOAD,
                              p4est->mpicomm, &requests[j]);
          P4EST_CHECK_MPI (mpiret);
        }
        else {
          /* will not receive load, close this request */
          requests[j] = MPI_REQUEST_NULL;
          --request_count;
        }
      }
      else {
        /* received load, close this request */
        P4EST_ASSERT (peer->expect_count > 0);
        peer->have_load = 1;
        requests[j] = MPI_REQUEST_NULL;
        --request_count;

        /* process incoming quadrants to interleave with communication */
        qarray = &peer->recv_first;
        qcount = qarray->elem_count;
        P4EST_ASSERT (peer->expect_count == qcount);
        q = p4est_array_index (qarray, 0);
        qtree = (int) q->user_data;
        P4EST_ASSERT (first_tree <= qtree && qtree <= last_tree);
        for (k = 1; k < qcount; ++k) {
          s = p4est_array_index (qarray, k);
          P4EST_ASSERT ((int) s->user_data == qtree);
        }
        tree = p4est_array_index (p4est->trees, qtree);
        P4EST_ASSERT (tree->quadrants->elem_count > 0);
        p4est_tree_compute_overlap (tree, qarray, &peer->send_second);

        for (k = 0; k < qarray->elem_count; ++k) {
          s = p4est_array_index (qarray, k);
          printf ("[%d] Tree %d inq 0x%x 0x%x %d\n", rank, qtree,
                  s->x, s->y, s->level);
        }
        for (k = 0; k < peer->send_second.elem_count; ++k) {
          s = p4est_array_index (&peer->send_second, k);
          printf ("[%d] Tree %d s2q 0x%x 0x%x %d\n", rank, qtree,
                  s->x, s->y, s->level);
        }

        /* remove what has already been in first send and send overlap */
      }
    }
  }

  /* receive second round */

  /* merge received quadrants */
  for (j = 0; j < p4est->mpisize; ++j) {
    peer = p4est_array_index (peers, j);
    if (!peer->will_receive || peer->expect_count == 0) {
      continue;
    }
    qarray = &peer->recv_first;
    qcount = qarray->elem_count;
    P4EST_ASSERT (qcount = peer->expect_count);
    q = p4est_array_index (qarray, 0);
    qtree = (int) q->user_data;
    P4EST_ASSERT (first_tree <= qtree && qtree <= last_tree);
    tree = p4est_array_index (p4est->trees, qtree);
    treecount = tree->quadrants->elem_count;
    p4est_array_resize (tree->quadrants, treecount + qcount);
    for (k = 0; k < qcount; ++k) {
      s = p4est_array_index (qarray, k);
      P4EST_ASSERT ((int) s->user_data == qtree);
      q = p4est_array_index (tree->quadrants, treecount + k);
      *q = *s;
      ++tree->quadrants_per_level[q->level];
      tree->maxlevel = P4EST_MAX (q->level, tree->maxlevel);
      ++p4est->local_num_quadrants;
      p4est_quadrant_init_data (p4est, qtree, q, init_fn);
      /*
        printf ("[%d] MR %d quad 0x%x 0x%x %d\n", rank, qtree,
        q->x, q->y, q->level);
      */
    }
    p4est_array_sort (tree->quadrants, p4est_quadrant_compare);
    /*
    printf ("[%d] tree %d MR %d %d %d\n", rank, qtree,
            treecount, qcount, tree->quadrants->elem_count);
    */
  }
  
  /* rebalance and clamp result back to original tree boundaries */
  p4est->local_num_quadrants = 0;
  for (j = first_tree; j <= last_tree; ++j) {
    tree = p4est_array_index (p4est->trees, j);
    treecount = tree->quadrants->elem_count;
    if ((j > first_tree || (mylow.x == 0 && mylow.y == 0)) &&
        (j < last_tree || (nextlow.x == 0 && nextlow.y == 0))) {
      p4est->local_num_quadrants += treecount;
      continue;
    }
    p4est_balance_subtree (p4est, tree, j, init_fn);
    treecount = tree->quadrants->elem_count;
    /*
    printf ("[%d] tree %d BL %d %d\n", rank, j,
            treecount, tree->quadrants->elem_count);
    */
    /* figure out the new elements outside the original tree */
    first_index = 0;
    last_index = treecount - 1;
    if (j == first_tree) {
      for (first_index = 0; first_index < treecount; ++first_index) {
        q = p4est_array_index (tree->quadrants, first_index);
        if (p4est_quadrant_compare (q, &mylow) >= 0 ||
            (q->x == mylow.x && q->y == mylow.y)) {
          break;
        }
      }
    }
    if (j == next_tree) {
      for (last_index = treecount - 1; last_index >= 0; --last_index) {
        q = p4est_array_index (tree->quadrants, last_index);
        /*
        printf ("[%d] LI %d quad 0x%x 0x%x %d\n", rank,
                last_index, q->x, q->y, q->level);
        */
        if (p4est_quadrant_compare (q, &nextlow) < 0) {
          break;
        }
      }
    }
    /*
    printf ("TC %d FI %d LI %d\n", treecount, first_index, last_index);
    */
    P4EST_ASSERT (first_index <= last_index);

    /* remove first part of tree */
    if (first_index > 0) {
      k = 0;
      while (first_index + k <= last_index) {
        q = p4est_array_index (tree->quadrants, k);
        s = p4est_array_index (tree->quadrants, first_index + k);
        p4est_quadrant_free_data (p4est, q);
        /*
        printf ("[%d] tree %d copy 0x%x 0x%x %d\n",
                rank, j, s->x, s->y, s->level);
        */
        *q = *s;
        ++k;
      }
    }
    /* remove last part of tree */
    qcount = last_index - first_index + 1;
    for (k = last_index + 1; k < treecount; ++k) {
      q = p4est_array_index (tree->quadrants, k);
      p4est_quadrant_free_data (p4est, q);
    }
    /*
      printf ("Resize %d %d\n", treecount, qcount);
    */
    p4est_array_resize (tree->quadrants, qcount);
    for (l = 0; l <= P4EST_MAXLEVEL; ++l) {
      tree->quadrants_per_level[l] = 0;
    }
    tree->maxlevel = 0;
    for (k = 0; k < qcount; ++k) {
      q = p4est_array_index (tree->quadrants, k);
      ++tree->quadrants_per_level[q->level];
      tree->maxlevel = P4EST_MAX (tree->maxlevel, q->level);
    }
    p4est->local_num_quadrants += qcount;
  }
#endif /* HAVE_MPI */

  /* loop over all local trees to finalize balance */
  for (j = first_tree; j <= last_tree; ++j) {
    tree = p4est_array_index (p4est->trees, j);

    /* final log message for this tree */
    if (p4est->nout != NULL) {
      fprintf (p4est->nout, "[%d] Done balance tree %d now %d\n",
               rank, j, tree->quadrants->elem_count);
    }
  }

#ifdef HAVE_MPI
  /* cleanup temporary storage */
  for (i = 0; i < p4est->mpisize; ++i) {
    peer = p4est_array_index (peers, i);
    p4est_array_reset (&peer->send_first);
    p4est_array_reset (&peer->recv_first);
    p4est_array_reset (&peer->send_second);
    p4est_array_reset (&peer->recv_second);
  }
  p4est_array_destroy (peers);
  P4EST_FREE (peer_boundaries);
  P4EST_FREE (requests);
  P4EST_FREE (statuses);
  P4EST_FREE (wait_indices);
#endif

  /* compute global number of quadrants */
  p4est_comm_count_quadrants (p4est);

  P4EST_ASSERT (p4est_is_valid (p4est));
}

/* EOF p4est.c */
