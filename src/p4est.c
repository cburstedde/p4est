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

#include <p4est.h>
#include <p4est_algorithms.h>
#include <p4est_communication.h>
#include <p4est_base.h>

/* require zlib header for adler32 checksums */
#ifdef HAVE_ZLIB_H
#include <zlib.h>
#endif

typedef struct
{
  int8_t              have_first_count, have_first_load;
  int8_t              have_second_count, have_second_load;
  int32_t             first_count, second_count;
  p4est_array_t       send_first, send_second, recv_both;
}
p4est_balance_peer_t;

static const int64_t initial_quadrants_per_processor = 15;
static const int    number_toread_quadrants = 32;
static const int    number_peer_windows = 5;

static const int8_t fully_owned_flag = 0x01;
static const int8_t any_face_flag = 0x02;

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

  P4EST_ASSERT (p4est_connectivity_is_valid (connectivity));

  P4EST_QUADRANT_INIT (&a);
  P4EST_QUADRANT_INIT (&b);

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
  p4est->global_last_quad_index = P4EST_ALLOC (int64_t, p4est->mpisize);
  P4EST_CHECK_ALLOC (p4est->global_last_quad_index);
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

#ifdef P4EST_HAVE_DEBUG
    int                 q;
    for (q = 0; q < tree->quadrants->elem_count; ++q) {
      p4est_quadrant_t   *quad = p4est_array_index (tree->quadrants, q);
      p4est_quadrant_free_data (p4est, quad);
    }
#endif

    p4est_array_destroy (tree->quadrants);
  }
  p4est_array_destroy (p4est->trees);

  if (p4est->user_data_pool != NULL) {
    p4est_mempool_destroy (p4est->user_data_pool);
  }
  p4est_mempool_destroy (p4est->quadrant_pool);

  P4EST_FREE (p4est->global_first_indices);
  P4EST_FREE (p4est->global_last_quad_index);

  P4EST_FREE (p4est);
}

p4est_t            *
p4est_copy (p4est_t * input, int copy_data)
{
  const int32_t       num_trees = input->connectivity->num_trees;
  const int32_t       first_tree = input->first_local_tree;
  const int32_t       last_tree = input->last_local_tree;
  int                 icount;
  int32_t             j, k;
  p4est_t            *p4est;
  p4est_tree_t       *itree, *ptree;
  p4est_quadrant_t   *iq, *pq;

  /* create a shallow copy and zero out dependent fields */
  p4est = P4EST_ALLOC (p4est_t, 1);
  P4EST_CHECK_ALLOC (p4est);
  memcpy (p4est, input, sizeof (p4est_t));
  p4est->global_last_quad_index = NULL;
  p4est->global_first_indices = NULL;
  p4est->trees = NULL;
  p4est->user_data_pool = NULL;
  p4est->quadrant_pool = NULL;

  /* allocate a user data pool if necessary and a quadrant pool */
  if (copy_data && p4est->data_size > 0) {
    p4est->user_data_pool = p4est_mempool_new (p4est->data_size);
  }
  else {
    p4est->data_size = 0;
  }
  p4est->quadrant_pool = p4est_mempool_new (sizeof (p4est_quadrant_t));

  /* copy quadrants for each tree */
  p4est->trees = p4est_array_new (sizeof (p4est_tree_t));
  p4est_array_resize (p4est->trees, num_trees);
  for (j = 0; j < num_trees; ++j) {
    itree = p4est_array_index (input->trees, j);
    ptree = p4est_array_index (p4est->trees, j);
    memcpy (ptree, itree, sizeof (p4est_tree_t));
    ptree->quadrants = p4est_array_new (sizeof (p4est_quadrant_t));
  }
  for (j = first_tree; j <= last_tree; ++j) {
    itree = p4est_array_index (input->trees, j);
    icount = itree->quadrants->elem_count;
    ptree = p4est_array_index (p4est->trees, j);
    p4est_array_resize (ptree->quadrants, icount);
    memcpy (ptree->quadrants->array, itree->quadrants->array,
            icount * sizeof (p4est_quadrant_t));
    if (p4est->data_size > 0) {
      for (k = 0; k < icount; ++k) {
        iq = p4est_array_index (itree->quadrants, k);
        pq = p4est_array_index (ptree->quadrants, k);
        pq->user_data = p4est_mempool_alloc (p4est->user_data_pool);
        memcpy (pq->user_data, iq->user_data, p4est->data_size);
      }
    }
    else {
      for (k = 0; k < icount; ++k) {
        pq = p4est_array_index (ptree->quadrants, k);
        pq->user_data = NULL;
      }
    }
  }

  /* allocate and copy global quadrant count */
  p4est->global_last_quad_index = P4EST_ALLOC (int64_t, p4est->mpisize);
  P4EST_CHECK_ALLOC (p4est->global_last_quad_index);
  memcpy (p4est->global_last_quad_index, input->global_last_quad_index,
          p4est->mpisize * sizeof (int64_t));

  /* allocate and copy global partition information */
  p4est->global_first_indices = P4EST_ALLOC (int32_t,
                                             3 * (p4est->mpisize + 1));
  P4EST_CHECK_ALLOC (p4est->global_first_indices);
  memcpy (p4est->global_first_indices, input->global_first_indices,
          3 * (p4est->mpisize + 1) * sizeof (int32_t));

  /* check for valid p4est and return */
  P4EST_ASSERT (p4est_is_valid (p4est));

  return p4est;
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

/** Check if the insulation layer of a quadrant overlaps anybody.
 * If yes, the quadrant itself is scheduled for sending.
 * Both quadrants are in the receiving tree's coordinates.
 * \param [in]  qtree       Tree id of the receiving tree.
 * \param [in]  inter_tree  Boolean flag to specify inter-tree communication.
 * \param [in]  q           The quadrant to be sent if there is overlap.
 * \param [in]  insul       An insulation quadrant of \a q.
 * \param [in,out]  first_peer  Lowest peer, will be updated.
 * \param [in,out]  last_peer   Highest peer, will be updated.
 */
static void
p4est_balance_schedule (p4est_t * p4est, p4est_array_t * peers,
                        int32_t qtree, int inter_tree,
                        const p4est_quadrant_t * q,
                        const p4est_quadrant_t * insul,
                        int *first_peer, int *last_peer)
{
  const int           rank = p4est->mpirank;
  int                 found, back, pos, scount;
  int32_t             owner, first_owner, last_owner;
  p4est_quadrant_t    ld, *s;
  p4est_balance_peer_t *peer;

  P4EST_QUADRANT_INIT (&ld);

  /* querying insul is equivalent to querying first descendent */
  first_owner = p4est_comm_find_owner (p4est, qtree, insul, rank);
  /* querying last descendent */
  p4est_quadrant_last_descendent (insul, &ld, P4EST_MAXLEVEL);
  last_owner = p4est_comm_find_owner (p4est, qtree, &ld, rank);

  /* send to all processors possibly intersecting insulation */
  for (owner = first_owner; owner <= last_owner; ++owner) {
    if (owner == rank && !inter_tree) {
      continue;
    }
    peer = p4est_array_index (peers, owner);
    /* avoid duplicates in the send array */
    found = 0;
    for (back = 0; back < 8; ++back) {
      pos = peer->send_first.elem_count - back - 1;
      if (pos < 0) {
        break;
      }
      s = p4est_array_index (&peer->send_first, pos);
      if (p4est_quadrant_is_equal (s, q) && (int32_t) s->user_data == qtree) {
        found = 1;
        break;
      }
    }
    if (found) {
      continue;
    }

    /* copy quadrant into shipping list */
    scount = peer->send_first.elem_count;
    p4est_array_resize (&peer->send_first, scount + 1);
    s = p4est_array_index (&peer->send_first, scount);
    *s = *q;
    s->user_data = (void *) qtree;      /* piggy back tree id  with quadrant */

    /* update lowest and highest peer */
    *first_peer = P4EST_MIN (owner, *first_peer);
    *last_peer = P4EST_MAX (owner, *last_peer);
  }
}

static void
p4est_balance_response (p4est_t * p4est, int32_t peer_id,
                        p4est_balance_peer_t * peer)
{
#ifdef P4EST_HAVE_DEBUG
  const int32_t       first_tree = p4est->first_local_tree;
  const int32_t       last_tree = p4est->last_local_tree;
#endif /* P4EST_HAVE_DEBUG */
  int32_t             k, qcount, qtree;
  int32_t             prev, num_receive_trees, nt;
  int32_t            *pi;
  p4est_array_t      *qarray, tree_array;
  p4est_quadrant_t   *q;

  qarray = &peer->recv_both;
  qcount = qarray->elem_count;
  P4EST_ASSERT (peer->first_count == qcount);

  /* build list of received tree ids */
  prev = -1;
  num_receive_trees = 0;
  p4est_array_init (&tree_array, sizeof (int32_t));
  for (k = 0; k < qcount; ++k) {
    q = p4est_array_index (qarray, k);
    qtree = (int32_t) q->user_data;
    P4EST_ASSERT (first_tree <= qtree && qtree <= last_tree);
    P4EST_ASSERT (qtree >= prev);
    if (qtree > prev) {
      p4est_array_resize (&tree_array, num_receive_trees + 1);
      pi = p4est_array_index (&tree_array, num_receive_trees);
      *pi = qtree;
      ++num_receive_trees;
      prev = qtree;
    }
  }
  if (p4est->nout != NULL) {
    fprintf (p4est->nout, "[%d] first load from %d into %d trees\n",
             p4est->mpirank, peer_id, num_receive_trees);
  }
  P4EST_ASSERT (num_receive_trees == tree_array.elem_count);

  /* loop to the trees to receive into and update overlap quadrants */
  for (nt = 0; nt < num_receive_trees; ++nt) {
    pi = p4est_array_index (&tree_array, nt);
    qtree = *pi;

    /* compute overlap quadrants */
    p4est_tree_compute_overlap (p4est, qtree, qarray, &peer->send_second);
  }
  p4est_tree_uniqify_overlap (&peer->send_first, &peer->send_second);
  p4est_array_reset (&tree_array);
}

void
p4est_balance (p4est_t * p4est, p4est_init_t init_fn)
{
  const int           rank = p4est->mpirank;
  int                 data_pool_size, all_incount, all_outcount;
  int                 k, l, ctree;
  int                 any_face, face_contact[4];
  int                 any_quad, quad_contact[4];
  int                 tree_fully_owned, transform;
  int                 first_index, last_index;
  int                 which;
  int8_t              face, corner;
  int8_t             *tree_flags;
  int32_t             i, j;
  int32_t             qtree;
  int32_t             qh, rh;
  int32_t             treecount, qcount, qbytes, offset, obytes;
  int32_t             first_tree, last_tree, next_tree;
  int32_t             first_peer, last_peer, rank_in_peers, over_peer_count;
  p4est_array_t      *peers, *qarray, corner_info;
  p4est_balance_peer_t *peer;
  p4est_tree_t       *tree;
  p4est_quadrant_t    mylow, nextlow;
  p4est_quadrant_t    tosend, insulq, tempq;
  p4est_quadrant_t   *q, *s;
  p4est_connectivity_t *conn = p4est->connectivity;
  p4est_corner_info_t *ci;
#ifdef HAVE_MPI
#ifdef P4EST_HAVE_DEBUG
  unsigned            checksum;
  p4est_array_t       checkarray;
  int32_t             ltotal[2], gtotal[2];
#endif /* P4EST_HAVE_DEBUG */
  const int           twopeerw = 2 * number_peer_windows;
  int                 mpiret;
  int                 first_bound, last_bound;
  int                 request_first_count, request_second_count, outcount;
  int                 request_send_count, total_send_count, total_recv_count;
  int                 lastw, nwin;
  int                 send_zero[2], send_load[2];
  int                 recv_zero[2], recv_load[2];
  int                *wait_indices;
  int32_t             prev, start, end;
  int32_t             length, shortest_window, shortest_length;
  int32_t             peer_windows[twopeerw];
  int32_t            *peer_boundaries;
  MPI_Request        *requests_first, *requests_second;
  MPI_Request        *send_requests_first_count, *send_requests_first_load;
  MPI_Request        *send_requests_second_count, *send_requests_second_load;
#endif /* HAVE_MPI */

  P4EST_ASSERT (p4est_is_valid (p4est));

  /* prepare sanity checks */
  if (p4est->user_data_pool != NULL) {
    data_pool_size = p4est->user_data_pool->elem_count;
  }

  P4EST_QUADRANT_INIT (&mylow);
  P4EST_QUADRANT_INIT (&nextlow);
  P4EST_QUADRANT_INIT (&tosend);
  P4EST_QUADRANT_INIT (&insulq);
  P4EST_QUADRANT_INIT (&tempq);

  /* tree status flags (max 8 per tree) */
  tree_flags = P4EST_ALLOC (int8_t, conn->num_trees);
  P4EST_CHECK_ALLOC (tree_flags);
  for (i = 0; i < conn->num_trees; ++i) {
    tree_flags[i] = 0x00;
  }

#ifdef HAVE_MPI
  /* will contain first and last peer (inclusive) for each processor */
  peer_boundaries = P4EST_ALLOC (int32_t, twopeerw * p4est->mpisize);
  P4EST_CHECK_ALLOC (peer_boundaries);
  /* request and status buffers for receive operations */
  requests_first = P4EST_ALLOC (MPI_Request, 6 * p4est->mpisize);
  P4EST_CHECK_ALLOC (requests_first);
  requests_second = requests_first + 1 * p4est->mpisize;
  send_requests_first_count = requests_first + 2 * p4est->mpisize;
  send_requests_first_load = requests_first + 3 * p4est->mpisize;
  send_requests_second_count = requests_first + 4 * p4est->mpisize;
  send_requests_second_load = requests_first + 5 * p4est->mpisize;
  for (i = 0; i < p4est->mpisize; ++i) {
    requests_first[i] = requests_second[i] = MPI_REQUEST_NULL;
    send_requests_first_count[i] = MPI_REQUEST_NULL;
    send_requests_first_load[i] = MPI_REQUEST_NULL;
    send_requests_second_count[i] = MPI_REQUEST_NULL;
    send_requests_second_load[i] = MPI_REQUEST_NULL;
  }
  /* index buffer for call to waitsome */
  wait_indices = P4EST_ALLOC (int, 4 * p4est->mpisize);
  P4EST_CHECK_ALLOC (wait_indices);
#ifdef P4EST_HAVE_DEBUG
  p4est_array_init (&checkarray, 4);
#endif /* P4EST_HAVE_DEBUG */
#endif /* HAVE_MPI */

  /* allocate per peer storage and initialize requests */
  peers = p4est_array_new (sizeof (p4est_balance_peer_t));
  p4est_array_resize (peers, p4est->mpisize);
  for (i = 0; i < p4est->mpisize; ++i) {
    peer = p4est_array_index (peers, i);
    p4est_array_init (&peer->send_first, sizeof (p4est_quadrant_t));
    p4est_array_init (&peer->send_second, sizeof (p4est_quadrant_t));
    p4est_array_init (&peer->recv_both, sizeof (p4est_quadrant_t));
    peer->first_count = peer->second_count = 0;
    peer->have_first_count = peer->have_first_load = 0;
    peer->have_second_count = peer->have_second_load = 0;
  }
  p4est_array_init (&corner_info, sizeof (p4est_corner_info_t));

  /* compute first quadrants on finest level for comparison for me and next */
  first_peer = p4est->mpisize;
  last_peer = -1;
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
  all_incount = 0;
  for (j = first_tree; j <= last_tree; ++j) {
    any_face = 0;
    for (face = 0; face < 4; ++face) {
      face_contact[face] = (conn->tree_to_tree[4 * j + face] != j);
      any_face = any_face || face_contact[face];
    }
    if (any_face) {
      tree_flags[j] |= any_face_flag;
    }
    tree = p4est_array_index (p4est->trees, j);
    all_incount += tree->quadrants->elem_count;

    /* initial log message for this tree */
    if (p4est->nout != NULL) {
      fprintf (p4est->nout, "[%d] Into balance tree %d with %d\n",
               rank, j, tree->quadrants->elem_count);
    }

    /* local balance first pass */
    p4est_balance_subtree (p4est, tree, j, init_fn);
    treecount = tree->quadrants->elem_count;
    if (p4est->nout != NULL) {
      fprintf (p4est->nout, "[%d] Balance tree %d A %d\n",
               rank, j, treecount);
    }

    /* check if this tree is not shared with other processors */
    tree_fully_owned = 0;
    if ((j > first_tree || (mylow.x == 0 && mylow.y == 0)) &&
        (j < last_tree || (nextlow.x == 0 && nextlow.y == 0))) {
      /* all quadrants in this tree are owned by me */
      tree_fully_owned = 1;
      tree_flags[j] |= fully_owned_flag;
      if (!any_face) {
        /* this tree is isolated, no balance between trees */
        continue;
      }
    }

    /* identify boundary quadrants and prepare them to be sent */
    for (i = 0; i < treecount; ++i) {
      /* this quadrant may be on the boundary with a range of processors */
      q = p4est_array_index (tree->quadrants, i);
      qh = (1 << (P4EST_MAXLEVEL - q->level));
      if (tree_fully_owned) {
        /* need only to consider boundary quadrants */
        any_quad =
          (face_contact[0] && q->y == 0) ||
          (face_contact[1] && q->x == rh - qh) ||
          (face_contact[2] && q->y == rh - qh) ||
          (face_contact[3] && q->x == 0);
        if (!any_quad) {
          continue;
        }
      }
      for (k = 0; k < 3; ++k) {
        for (l = 0; l < 3; ++l) {
          which = k * 3 + l;    /* 0..8 */
          /* exclude myself from the queries */
          if (which == 4) {
            continue;
          }
          /* may modify insulq below, never modify q itself! */
          insulq = *q;
          insulq.x += (l - 1) * qh;
          insulq.y += (k - 1) * qh;

          /* check boundary status of s */
          quad_contact[0] = (insulq.y < 0);
          quad_contact[1] = (insulq.x >= rh);
          quad_contact[2] = (insulq.y >= rh);
          quad_contact[3] = (insulq.x < 0);
          if (quad_contact[0] || quad_contact[1] ||
              quad_contact[2] || quad_contact[3]) {
            /* this quadrant is relevant for inter-tree balancing */
            if ((quad_contact[0] || quad_contact[2]) &&
                (quad_contact[1] || quad_contact[3])) {
              /* this quadrant goes across a corner */
              for (corner = 0; corner < 4; ++corner) {
                if (quad_contact[(corner + 3) % 4] && quad_contact[corner]) {
                  break;
                }
              }
              p4est_find_corner_info (conn, j, corner, &corner_info);
              for (ctree = 0; ctree < corner_info.elem_count; ++ctree) {
                ci = p4est_array_index (&corner_info, ctree);
                tosend = *q;
                p4est_quadrant_corner (&tosend, ci->ncorner, 0);
                p4est_quadrant_corner (&insulq, ci->ncorner, 1);
                p4est_balance_schedule (p4est, peers, ci->ntree, 1,
                                        &tosend, &insulq,
                                        &first_peer, &last_peer);
              }
            }
            else {
              /* this quadrant goes across a face */
              for (face = 0; face < 4; ++face) {
                if (quad_contact[face] && face_contact[face]) {
                  qtree = conn->tree_to_tree[4 * j + face];
                  P4EST_ASSERT (qtree != j);
                  break;
                }
              }
              if (face == 4) {
                /* this quadrant ran across a face with no neighbor */
                continue;
              }
              /* transform both q and insulq into the neighbor's coordinates */
              transform = p4est_find_face_transform (conn, j, face);
              tempq = *q;
              p4est_quadrant_translate (&tempq, face);
              p4est_quadrant_transform (&tempq, &tosend, transform);
              p4est_quadrant_translate (&insulq, face);
              p4est_quadrant_transform (&insulq, &tempq, transform);
              p4est_balance_schedule (p4est, peers, qtree, 1,
                                      &tosend, &tempq,
                                      &first_peer, &last_peer);
            }
          }
          else {
            /* no inter-tree contact */
            p4est_balance_schedule (p4est, peers, j, 0,
                                    q, &insulq, &first_peer, &last_peer);
          }
        }
      }
    }
  }
  if (last_peer < 0) {
    P4EST_ASSERT (first_peer == p4est->mpisize);
    first_peer = last_peer = rank;
  }
  P4EST_ASSERT (first_peer <= last_peer);
  P4EST_ASSERT (0 <= first_peer && last_peer < p4est->mpisize);
  rank_in_peers = (first_peer <= rank && rank <= last_peer) ? 1 : 0;
  over_peer_count = last_peer - first_peer + 1 - rank_in_peers;
  P4EST_ASSERT (0 <= over_peer_count && over_peer_count < p4est->mpisize);
  if (p4est->mpisize == 1) {
    P4EST_ASSERT (first_peer == rank && last_peer == rank);
  }

#ifdef HAVE_MPI
  /* compute information about peers to inform receivers */
  for (i = 0; i < twopeerw; ++i) {
    peer_windows[i] = -1;
  }
  /* find a maximum of npw-1 empty ranges with (start, end) */
  lastw = number_peer_windows - 1;
  prev = -1;
  nwin = 0;
  for (j = 0; j < p4est->mpisize; ++j) {
    peer = p4est_array_index (peers, j);
    if (peer->send_first.elem_count == 0 || j == rank) {
      continue;
    }
    if (prev == -1) {
      prev = j;
      continue;
    }
    if (prev < j - 1) {
      length = j - 1 - prev;
#ifdef P4EST_HAVE_DEBUG
      if (p4est->nout != NULL) {
        fprintf (p4est->nout,
                 "[%d] found empty range prev %d j %d length %d\n",
                 rank, prev, j, length);
      }
#endif /* P4EST_HAVE_DEBUG */
      /* claim unused window */
      for (i = 0; i < number_peer_windows; ++i) {
        if (peer_windows[2 * i] == -1) {
          peer_windows[2 * i] = prev + 1;
          peer_windows[2 * i + 1] = j - 1;
          break;
        }
      }
      P4EST_ASSERT (i < number_peer_windows);
      nwin = i + 1;
      /* if all ranges are used, remove the shortest */
      if (nwin == number_peer_windows) {
        nwin = lastw;
        shortest_window = -1;
        shortest_length = p4est->mpisize + 1;
        for (i = 0; i < number_peer_windows; ++i) {
          length = peer_windows[2 * i + 1] - peer_windows[2 * i] + 1;
          if (length < shortest_length) {
            shortest_window = i;
            shortest_length = length;
          }
        }
        P4EST_ASSERT (shortest_window >= 0 && shortest_window <= lastw);
        if (shortest_window < lastw) {
          peer_windows[2 * shortest_window] = peer_windows[2 * lastw];
          peer_windows[2 * shortest_window + 1] = peer_windows[2 * lastw + 1];
        }
        peer_windows[2 * lastw] = -1;
        peer_windows[2 * lastw + 1] = -1;
      }
    }
    prev = j;
  }
  P4EST_ASSERT (nwin >= 0 && nwin < number_peer_windows);
  /* bubble sort empty ranges by start rank */
  for (i = nwin - 1; i >= 0; --i) {
    for (j = 0; j < i; ++j) {
      if (peer_windows[2 * j] > peer_windows[2 * (j + 1)]) {
        start = peer_windows[2 * j];
        end = peer_windows[2 * j + 1];
        peer_windows[2 * j] = peer_windows[2 * (j + 1)];
        peer_windows[2 * j + 1] = peer_windows[2 * (j + 1) + 1];
        peer_windows[2 * (j + 1)] = start;
        peer_windows[2 * (j + 1) + 1] = end;
      }
    }
  }
#ifdef P4EST_HAVE_DEBUG
  for (i = 0; i < nwin; ++i) {
    P4EST_ASSERT (peer_windows[2 * i] <= peer_windows[2 * i + 1]);
    if (i < nwin - 1) {
      P4EST_ASSERT (peer_windows[2 * i + 1] < peer_windows[2 * (i + 1)] - 1);
    }
  }
  for (i = nwin; i < number_peer_windows; ++i) {
    P4EST_ASSERT (peer_windows[2 * i] == -1);
    P4EST_ASSERT (peer_windows[2 * i + 1] == -1);
  }
  if (p4est->nout != NULL) {
    for (i = 0; i < nwin; ++i) {
      fprintf (p4est->nout,
               "[%d] empty range %d from %d to %d\n",
               rank, i, peer_windows[2 * i], peer_windows[2 * i + 1]);
    }
  }
#endif /* P4EST_HAVE_DEBUG */
  /* compute windows from empty ranges */
  peer_windows[2 * nwin + 1] = last_peer;
  for (i = nwin; i > 0; --i) {
    peer_windows[2 * i] = peer_windows[2 * i - 1] + 1;
    peer_windows[2 * i - 1] = peer_windows[2 * (i - 1)] - 1;
  }
  peer_windows[0] = first_peer;
  ++nwin;
#ifdef P4EST_HAVE_DEBUG
  for (i = 0; i < nwin; ++i) {
    P4EST_ASSERT (peer_windows[2 * i] <= peer_windows[2 * i + 1]);
    if (i < nwin - 1) {
      P4EST_ASSERT (peer_windows[2 * i + 1] < peer_windows[2 * (i + 1)] - 1);
    }
  }
  for (i = nwin; i < number_peer_windows; ++i) {
    P4EST_ASSERT (peer_windows[2 * i] == -1);
    P4EST_ASSERT (peer_windows[2 * i + 1] == -1);
  }
  if (p4est->nout != NULL) {
    for (i = 0; i < nwin; ++i) {
      fprintf (p4est->nout, "[%d] window %d from %d to %d\n",
               rank, i, peer_windows[2 * i], peer_windows[2 * i + 1]);
    }
  }
#endif /* P4EST_HAVE_DEBUG */

  /* distribute information about peer windows to inform receivers */
  if (p4est->mpicomm != MPI_COMM_NULL) {
    mpiret = MPI_Allgather (peer_windows, twopeerw, MPI_INT,
                            peer_boundaries, twopeerw, MPI_INT,
                            p4est->mpicomm);
    P4EST_CHECK_MPI (mpiret);
  }
  if (p4est->nout != NULL) {
    fprintf (p4est->nout, "[%d] Peers first %d last %d count %d\n",
             rank, first_peer, last_peer, over_peer_count);
  }

  /*
   * loop over all peers and send first round of quadrants
   * for intra-tree balancing, each load is contained in one tree
   */
  total_send_count = total_recv_count = 0;
  request_first_count = request_second_count = request_send_count = 0;
  send_zero[0] = send_load[0] = recv_zero[0] = recv_load[0] = 0;
  send_zero[1] = send_load[1] = recv_zero[1] = recv_load[1] = 0;
  for (j = first_peer; j <= last_peer; ++j) {
    if (j == rank) {
      continue;
    }
    peer = p4est_array_index (peers, j);
    qcount = peer->send_first.elem_count;

    /* check windows here for now, eventually merge into outer loop */
    for (i = 0; i < nwin - 1; ++i) {
      if (j > peer_windows[2 * i + 1] && j < peer_windows[2 * (i + 1)]) {
        break;
      }
    }
    if (i < nwin - 1) {
      P4EST_ASSERT (qcount == 0);
      continue;
    }

    /* first send number of quadrants to be expected */
    if (qcount > 0) {
      if (p4est->nout != NULL) {
        fprintf (p4est->nout, "[%d] Balance A send %d quadrants to %d\n",
                 rank, qcount, j);
      }
      ++send_load[0];
    }
    else {
      ++send_zero[0];
    }
    mpiret = MPI_Isend (&qcount, 1, MPI_INT,
                        j, P4EST_COMM_BALANCE_FIRST_COUNT,
                        p4est->mpicomm, &send_requests_first_count[j]);
    P4EST_CHECK_MPI (mpiret);
    ++request_send_count;

    /* sort and send the actual quadrants and post receive for reply */
    if (qcount > 0) {
      p4est_array_sort (&peer->send_first, p4est_quadrant_compare_piggy);
#ifdef P4EST_HAVE_DEBUG
      checksum = p4est_quadrant_checksum (&peer->send_first, &checkarray, 0);
      if (p4est->nout != NULL) {
        fprintf (p4est->nout, "[%d] Balance A send checksum %x to %d\n",
                 rank, checksum, j);
      }
#endif /* P4EST_HAVE_DEBUG */
      total_send_count += qcount;
      qbytes = qcount * sizeof (p4est_quadrant_t);
      mpiret = MPI_Isend (peer->send_first.array, qbytes, MPI_CHAR,
                          j, P4EST_COMM_BALANCE_FIRST_LOAD,
                          p4est->mpicomm, &send_requests_first_load[j]);
      P4EST_CHECK_MPI (mpiret);
      ++request_send_count;
      mpiret = MPI_Irecv (&peer->second_count, 1, MPI_INT,
                          j, P4EST_COMM_BALANCE_SECOND_COUNT,
                          p4est->mpicomm, &requests_second[j]);
      P4EST_CHECK_MPI (mpiret);
      ++request_second_count;
    }
  }

  /* find out who is sending to me and receive quadrant counts */
  for (j = 0; j < p4est->mpisize; ++j) {
    if (j == rank) {
      continue;
    }
    peer = p4est_array_index (peers, j);
    for (i = 0; i < number_peer_windows; ++i) {
      first_bound = peer_boundaries[twopeerw * j + 2 * i];
      last_bound = peer_boundaries[twopeerw * j + 2 * i + 1];
      if (first_bound <= rank && rank <= last_bound) {
        break;
      }
    }
    if (i == number_peer_windows) {
      continue;
    }
    peer = p4est_array_index (peers, j);
    ++request_first_count;

    /* processor j is sending to me */
    mpiret = MPI_Irecv (&peer->first_count, 1, MPI_INT,
                        j, P4EST_COMM_BALANCE_FIRST_COUNT,
                        p4est->mpicomm, &requests_first[j]);
    P4EST_CHECK_MPI (mpiret);
  }

  /* wait for quadrant counts and post receive and send for quadrants */
  while (request_first_count > 0) {
    mpiret = MPI_Waitsome (p4est->mpisize, requests_first,
                           &outcount, wait_indices, MPI_STATUSES_IGNORE);
    P4EST_CHECK_MPI (mpiret);
    P4EST_ASSERT (outcount != MPI_UNDEFINED);
    P4EST_ASSERT (outcount > 0);
    for (i = 0; i < outcount; ++i) {
      /* retrieve sender's rank */
      j = wait_indices[i];
      wait_indices[i] = -1;
      P4EST_ASSERT (0 <= j && j < p4est->mpisize);

      /* check if we are in receiving count or load */
      peer = p4est_array_index (peers, j);
      P4EST_ASSERT (!peer->have_first_load);
      if (!peer->have_first_count) {
        peer->have_first_count = 1;
        qcount = peer->first_count;
        if (qcount > 0) {
          /* received nonzero count, post receive for load */
          if (p4est->nout != NULL) {
            fprintf (p4est->nout,
                     "[%d] Balance A recv %d quadrants from %d\n",
                     rank, qcount, j);
          }
          P4EST_ASSERT (peer->recv_both.elem_count == 0);
          p4est_array_resize (&peer->recv_both, qcount);
          total_recv_count += qcount;
          qbytes = qcount * sizeof (p4est_quadrant_t);
          mpiret = MPI_Irecv (peer->recv_both.array, qbytes, MPI_CHAR,
                              j, P4EST_COMM_BALANCE_FIRST_LOAD,
                              p4est->mpicomm, &requests_first[j]);
          P4EST_CHECK_MPI (mpiret);
          ++recv_load[0];
        }
        else {
          /* will not receive load, close this request */
          requests_first[j] = MPI_REQUEST_NULL;
          --request_first_count;
          ++recv_zero[0];
        }
      }
      else {
        /* received load, close this request */
        P4EST_ASSERT (peer->first_count > 0);
        peer->have_first_load = 1;
        requests_first[j] = MPI_REQUEST_NULL;
        --request_first_count;
#ifdef P4EST_HAVE_DEBUG
        checksum = p4est_quadrant_checksum (&peer->recv_both, &checkarray, 0);
        if (p4est->nout != NULL) {
          fprintf (p4est->nout, "[%d] Balance A recv checksum %x from %d\n",
                   rank, checksum, j);
        }
#endif /* P4EST_HAVE_DEBUG */

        /* process incoming quadrants to interleave with communication */
        p4est_balance_response (p4est, j, peer);
        qcount = peer->send_second.elem_count;
        if (qcount > 0) {
          if (p4est->nout != NULL) {
            fprintf (p4est->nout, "[%d] Balance B send %d quadrants to %d\n",
                     rank, qcount, j);
          }
          ++send_load[1];
        }
        else {
          ++send_zero[1];
        }
        mpiret = MPI_Isend (&qcount, 1, MPI_INT,
                            j, P4EST_COMM_BALANCE_SECOND_COUNT,
                            p4est->mpicomm, &send_requests_second_count[j]);
        P4EST_CHECK_MPI (mpiret);
        ++request_send_count;
        if (qcount > 0) {
#ifdef P4EST_HAVE_DEBUG
          checksum =
            p4est_quadrant_checksum (&peer->send_second, &checkarray, 0);
          if (p4est->nout != NULL) {
            fprintf (p4est->nout, "[%d] Balance B send checksum %x to %d\n",
                     rank, checksum, j);
          }
#endif /* P4EST_HAVE_DEBUG */
          total_send_count += qcount;
          qbytes = qcount * sizeof (p4est_quadrant_t);
          mpiret = MPI_Isend (peer->send_second.array, qbytes, MPI_CHAR,
                              j, P4EST_COMM_BALANCE_SECOND_LOAD,
                              p4est->mpicomm, &send_requests_second_load[j]);
          P4EST_CHECK_MPI (mpiret);
          ++request_send_count;
        }
      }
    }
  }
#endif /* HAVE_MPI */

  /* simulate send and receive with myself across tree boundaries */
  peer = p4est_array_index (peers, rank);
  p4est_array_sort (&peer->send_first, p4est_quadrant_compare_piggy);
  offset = peer->first_count = peer->send_first.elem_count;
  obytes = offset * sizeof (p4est_quadrant_t);
  qarray = &peer->recv_both;
  p4est_array_resize (qarray, offset);
  memcpy (qarray->array, peer->send_first.array, obytes);
  p4est_balance_response (p4est, rank, peer);
  qcount = peer->second_count = peer->send_second.elem_count;
  qbytes = qcount * sizeof (p4est_quadrant_t);
  p4est_array_resize (qarray, offset + qcount);
  memcpy (qarray->array + obytes, peer->send_second.array, qbytes);

#ifdef HAVE_MPI
  /* receive second round appending to the same receive buffer */
  while (request_second_count > 0) {
    mpiret = MPI_Waitsome (p4est->mpisize, requests_second,
                           &outcount, wait_indices, MPI_STATUSES_IGNORE);
    P4EST_CHECK_MPI (mpiret);
    P4EST_ASSERT (outcount != MPI_UNDEFINED);
    P4EST_ASSERT (outcount > 0);
    for (i = 0; i < outcount; ++i) {
      /* retrieve sender's rank */
      j = wait_indices[i];
      wait_indices[i] = -1;
      P4EST_ASSERT (0 <= j && j < p4est->mpisize);

      /* check if we are in receiving count or load */
      peer = p4est_array_index (peers, j);
      P4EST_ASSERT (!peer->have_second_load);
      if (!peer->have_second_count) {
        peer->have_second_count = 1;
        qcount = peer->second_count;
        if (qcount > 0) {
          /* received nonzero count, post receive for load */
          if (p4est->nout != NULL) {
            fprintf (p4est->nout,
                     "[%d] Balance B recv %d quadrants from %d\n",
                     rank, qcount, j);
          }
          offset = peer->recv_both.elem_count;
          P4EST_ASSERT (offset == peer->first_count);
          obytes = offset * sizeof (p4est_quadrant_t);
          p4est_array_resize (&peer->recv_both, offset + qcount);
          total_recv_count += qcount;
          qbytes = qcount * sizeof (p4est_quadrant_t);
          mpiret = MPI_Irecv (peer->recv_both.array + obytes, qbytes,
                              MPI_CHAR, j, P4EST_COMM_BALANCE_SECOND_LOAD,
                              p4est->mpicomm, &requests_second[j]);
          P4EST_CHECK_MPI (mpiret);
          ++recv_load[1];
        }
        else {
          /* will not receive load, close this request */
          requests_second[j] = MPI_REQUEST_NULL;
          --request_second_count;
          ++recv_zero[1];
        }
      }
      else {
        /* received load, close this request */
        P4EST_ASSERT (peer->second_count > 0);
        peer->have_second_load = 1;
        requests_second[j] = MPI_REQUEST_NULL;
        --request_second_count;
#ifdef P4EST_HAVE_DEBUG
        checksum = p4est_quadrant_checksum (&peer->recv_both, &checkarray,
                                            peer->first_count);
        if (p4est->nout != NULL) {
          fprintf (p4est->nout, "[%d] Balance B recv checksum %x from %d\n",
                   rank, checksum, j);
        }
#endif /* P4EST_HAVE_DEBUG */
      }
    }
  }

  /* print buffer statistics */
  if (p4est->nout != NULL) {
    fprintf (p4est->nout, "[%d] first send Z %d L %d recv Z %d L %d\n",
             rank, send_zero[0], send_load[0], recv_zero[0], recv_load[0]);
    fprintf (p4est->nout, "[%d] second send Z %d L %d recv Z %d L %d\n",
             rank, send_zero[1], send_load[1], recv_zero[1], recv_load[1]);
    fprintf (p4est->nout, "[%d] total send %d recv %d\n",
             rank, total_send_count, total_recv_count);
#ifdef P4EST_HAVE_DEBUG
    for (j = 0; j < p4est->mpisize; ++j) {
      peer = p4est_array_index (peers, j);
      if (peer->send_first.elem_count > 0 || peer->first_count > 0 ||
          peer->send_second.elem_count > 0 || peer->second_count > 0) {
        fprintf (p4est->nout,
                 "[%d] peer %d first S %d R %d second S %d R %d\n", rank, j,
                 peer->send_first.elem_count, peer->first_count,
                 peer->send_second.elem_count, peer->second_count);
      }
    }
#endif /* P4EST_HAVE_DEBUG */
  }
  if (number_peer_windows == 1) {
    P4EST_ASSERT (send_zero[0] + send_load[0] == over_peer_count);
    P4EST_ASSERT (recv_zero[1] + recv_load[1] == over_peer_count);
  }
#endif /* HAVE_MPI */

  /* merge received quadrants */
  for (j = 0; j < p4est->mpisize; ++j) {
    /* access peer information */
    peer = p4est_array_index (peers, j);
    qarray = &peer->recv_both;
    qcount = qarray->elem_count;
    P4EST_ASSERT (qcount == peer->first_count + peer->second_count);
    if (qcount == 0) {
      continue;
    }

    /* merge received quadrants into correct tree */
    for (k = 0; k < qcount; ++k) {
      s = p4est_array_index (qarray, k);
      P4EST_ASSERT (p4est_quadrant_is_extended (s));
      qtree = (int32_t) s->user_data;
      if (qtree < first_tree || qtree > last_tree) {
        /* this is a corner quadrant from the second pass of balance */
        P4EST_ASSERT (k >= peer->first_count);
        P4EST_ASSERT ((s->x < 0 && s->y < 0) || (s->x < 0 && s->y >= rh) ||
                      (s->x >= rh && s->y < 0) || (s->x >= rh && s->y >= rh));
        continue;
      }
      tree = p4est_array_index (p4est->trees, qtree);
      treecount = tree->quadrants->elem_count;
      p4est_array_resize (tree->quadrants, treecount + 1);
      q = p4est_array_index (tree->quadrants, treecount);
      *q = *s;
      ++tree->quadrants_per_level[q->level];
      tree->maxlevel = (int8_t) P4EST_MAX (tree->maxlevel, q->level);
      ++p4est->local_num_quadrants;
      p4est_quadrant_init_data (p4est, qtree, q, init_fn);
    }
  }

  /* rebalance and clamp result back to original tree boundaries */
  p4est->local_num_quadrants = 0;
  for (j = first_tree; j <= last_tree; ++j) {
    /* check if we are the only processor in an isolated tree */
    tree = p4est_array_index (p4est->trees, j);
    if ((tree_flags[j] & fully_owned_flag) &&
        !(tree_flags[j] & any_face_flag)) {
      p4est->local_num_quadrants += tree->quadrants->elem_count;
      continue;
    }

    /* we have most probably received quadrants, run sort and balance */
    p4est_array_sort (tree->quadrants, p4est_quadrant_compare);
    p4est_balance_subtree (p4est, tree, j, init_fn);
    if (p4est->nout != NULL) {
      fprintf (p4est->nout, "[%d] Balance tree %d B %d\n",
               rank, j, tree->quadrants->elem_count);
    }
    treecount = tree->quadrants->elem_count;

    /* figure out the new elements outside the original tree */
    for (first_index = 0; first_index < treecount; ++first_index) {
      q = p4est_array_index (tree->quadrants, first_index);
      P4EST_ASSERT (p4est_quadrant_is_extended (q));
      if (p4est_quadrant_is_valid (q)) {
        break;
      }
    }
    if (j == first_tree) {
      for (; first_index < treecount; ++first_index) {
        q = p4est_array_index (tree->quadrants, first_index);
        if (p4est_quadrant_compare (q, &mylow) >= 0 ||
            (q->x == mylow.x && q->y == mylow.y)) {
          break;
        }
      }
    }
    for (last_index = treecount - 1; last_index >= 0; --last_index) {
      q = p4est_array_index (tree->quadrants, last_index);
      P4EST_ASSERT (p4est_quadrant_is_extended (q));
      if (p4est_quadrant_is_valid (q)) {
        break;
      }
    }
    if (j == next_tree) {
      for (; last_index >= 0; --last_index) {
        q = p4est_array_index (tree->quadrants, last_index);
        if (p4est_quadrant_compare (q, &nextlow) < 0 &&
            (q->x != nextlow.x || q->y != nextlow.y)) {
          break;
        }
      }
    }
    P4EST_ASSERT (first_index <= last_index);

    /* remove first part of tree */
    if (first_index > 0) {
      k = 0;
      while (first_index + k <= last_index) {
        q = p4est_array_index (tree->quadrants, k);
        s = p4est_array_index (tree->quadrants, first_index + k);
        if (k < first_index) {
          p4est_quadrant_free_data (p4est, q);
        }
        *q = *s;
        ++k;
      }
      while (k < first_index) {
        q = p4est_array_index (tree->quadrants, k);
        p4est_quadrant_free_data (p4est, q);
        ++k;
      }
    }

    /* remove last part of tree */
    qcount = last_index - first_index + 1;
    for (k = last_index + 1; k < treecount; ++k) {
      q = p4est_array_index (tree->quadrants, k);
      p4est_quadrant_free_data (p4est, q);
    }

    P4EST_ASSERT (qcount <= treecount);
    p4est_array_resize (tree->quadrants, qcount);
    for (l = 0; l <= P4EST_MAXLEVEL; ++l) {
      tree->quadrants_per_level[l] = 0;
    }
    tree->maxlevel = 0;
    for (k = 0; k < qcount; ++k) {
      q = p4est_array_index (tree->quadrants, k);
      P4EST_ASSERT (p4est_quadrant_is_valid (q));
      ++tree->quadrants_per_level[q->level];
      tree->maxlevel = (int8_t) P4EST_MAX (tree->maxlevel, q->level);
    }
    p4est->local_num_quadrants += qcount;

    if (p4est->nout != NULL) {
      fprintf (p4est->nout, "[%d] Balance tree %d C %d\n",
               rank, j, tree->quadrants->elem_count);
    }
  }

#ifdef HAVE_MPI
  /* compute global sum of send and receive counts */
#ifdef P4EST_HAVE_DEBUG
  ltotal[0] = total_send_count;
  ltotal[1] = total_recv_count;
  if (p4est->mpicomm != MPI_COMM_NULL) {
    mpiret = MPI_Reduce (ltotal, gtotal, 2, MPI_INT,
                         MPI_SUM, 0, p4est->mpicomm);
    P4EST_CHECK_MPI (mpiret);
    if (rank == 0) {
      if (p4est->nout != NULL) {
        fprintf (p4est->nout,
                 "   global number of shipped quadrants %d\n", gtotal[0]);
      }
      P4EST_ASSERT (gtotal[0] == gtotal[1]);
    }
  }
  else {
    P4EST_ASSERT (ltotal[0] == 0 && ltotal[1] == 0);
  }
#endif /* P4EST_HAVE_DEBUG */

  /* wait for all send operations */
  if (request_send_count > 0) {
    mpiret = MPI_Waitall (4 * p4est->mpisize,
                          send_requests_first_count, MPI_STATUSES_IGNORE);
    P4EST_CHECK_MPI (mpiret);
  }
#endif /* HAVE_MPI */

  /* loop over all local trees to finalize balance */
  all_outcount = 0;
  for (j = first_tree; j <= last_tree; ++j) {
    tree = p4est_array_index (p4est->trees, j);
    all_outcount += tree->quadrants->elem_count;

    /* final log message for this tree */
    if (p4est->nout != NULL) {
      fprintf (p4est->nout, "[%d] Done balance tree %d now %d\n",
               rank, j, tree->quadrants->elem_count);
    }
  }

  /* cleanup temporary storage */
  P4EST_FREE (tree_flags);
  for (i = 0; i < p4est->mpisize; ++i) {
    peer = p4est_array_index (peers, i);
    p4est_array_reset (&peer->send_first);
    p4est_array_reset (&peer->send_second);
    p4est_array_reset (&peer->recv_both);
  }
  p4est_array_destroy (peers);
  p4est_array_reset (&corner_info);
#ifdef HAVE_MPI
  P4EST_FREE (peer_boundaries);
  P4EST_FREE (requests_first);  /* includes allocation for requests_second */
  P4EST_FREE (wait_indices);
#ifdef P4EST_HAVE_DEBUG
  p4est_array_reset (&checkarray);
#endif /* P4EST_HAVE_DEBUG */
#endif /* HAVE_MPI */

  /* compute global number of quadrants */
  p4est_comm_count_quadrants (p4est);

  /* some sanity checks */
  P4EST_ASSERT (all_outcount == p4est->local_num_quadrants);
  if (p4est->user_data_pool != NULL) {
    P4EST_ASSERT (data_pool_size + (all_outcount - all_incount) ==
                  p4est->user_data_pool->elem_count);
  }
  P4EST_ASSERT (p4est_is_valid (p4est));
}

/** Find the lowest position k in a sorted array such that array[k] >= target.
 * \return  Returns the matching position, or -1 if not found.
 */
static int
int64_find_lower_bound (int64_t target, const int64_t * array,
                        int count, int guess)
{
  int                 k_low, k_high;
  int64_t             cur;

  k_low = 0;
  k_high = count - 1;
  for (;;) {
    P4EST_ASSERT (k_low <= k_high);
    P4EST_ASSERT (0 <= k_low && k_low < count);
    P4EST_ASSERT (0 <= k_high && k_high < count);
    P4EST_ASSERT (k_low <= guess && guess <= k_high);

    /* compare two quadrants */
    cur = array[guess];

    /* check if guess is higher or equal target and there's room below it */
    if (target <= cur && (guess > 0 && target <= array[guess - 1])) {
      k_high = guess - 1;
      guess = (k_low + k_high + 1) / 2;
      continue;
    }

    /* check if guess is lower than target */
    if (target > cur) {
      k_low = guess + 1;
      if (k_low > k_high) {
        return -1;
      }
      guess = (k_low + k_high) / 2;
      continue;
    }

    /* otherwise guess is the correct position */
    break;
  }

  return guess;
}

void
p4est_partition (p4est_t * p4est, p4est_weight_t weight_fn)
{
#ifdef HAVE_MPI
  int                 mpiret;
  const int           num_procs = p4est->mpisize;
  const int           rank = p4est->mpirank;
  const int32_t       first_tree = p4est->first_local_tree;
  const int32_t       last_tree = p4est->last_local_tree;
  const int32_t       local_num_quadrants = p4est->local_num_quadrants;
  const int64_t       global_num_quadrants = p4est->global_num_quadrants;
  int                 i, p;
  int                 send_lowest, send_highest;
  int                 num_sends;
  int32_t             j, k, qcount;
  int32_t             cut, my_lowcut, my_highcut;
  int32_t            *num_quadrants_in_proc;
  int64_t             prev_quadrant, next_quadrant;
  int64_t             weight, weight_sum;
  int64_t             send_index, recv_low, recv_high;
  int64_t            *local_weights;    /* cumulative weights by quadrant */
  int64_t            *global_weight_sums;
  p4est_quadrant_t   *q;
  p4est_tree_t       *tree;
  MPI_Request        *send_requests, recv_requests[2];

  /* this function does nothing in a serial setup */
  if (p4est->mpicomm == MPI_COMM_NULL || p4est->mpisize == 1) {
    return;
  }

  /* allocate new quadrant distribution counts */
  num_quadrants_in_proc = P4EST_ALLOC (int32_t, num_procs);
  P4EST_CHECK_ALLOC (num_quadrants_in_proc);

  if (weight_fn == NULL) {
    /* Divide up the quadants equally */
    for (p = 0, next_quadrant = 0; p < num_procs; ++p) {
      prev_quadrant = next_quadrant;
      next_quadrant = (global_num_quadrants * (p + 1)) / num_procs;
      num_quadrants_in_proc[p] = (int32_t) (next_quadrant - prev_quadrant);
    }
  }
  else {
    /* do a weighted partition */
    local_weights = P4EST_ALLOC (int64_t, local_num_quadrants + 1);
    P4EST_CHECK_ALLOC (local_weights);
    global_weight_sums = P4EST_ALLOC (int64_t, num_procs + 1);
    P4EST_CHECK_ALLOC (global_weight_sums);

    /* linearly sum weights across all trees */
    k = 0;
    local_weights[0] = 0;
    for (j = first_tree; j <= last_tree; ++j) {
      tree = p4est_array_index (p4est->trees, j);
      for (i = 0; i < tree->quadrants->elem_count; ++i, ++k) {
        q = p4est_array_index (tree->quadrants, i);
        weight = weight_fn (p4est, j, q);
        P4EST_ASSERT (weight >= 0);
        local_weights[k + 1] = local_weights[k] + weight;
      }
    }
    P4EST_ASSERT (k == local_num_quadrants);
    weight_sum = local_weights[local_num_quadrants];

#if(0)
    printf ("[%d] local weight sum %lld\n", rank, weight_sum);
#endif

    /* distribute local weight sums */
    global_weight_sums[0] = 0;
    mpiret = MPI_Allgather (&weight_sum, 1, MPI_LONG_LONG,
                            &global_weight_sums[1], 1, MPI_LONG_LONG,
                            p4est->mpicomm);
    P4EST_CHECK_MPI (mpiret);

    /* adjust all arrays to reflect the global weight */
    for (i = 0; i < num_procs; ++i) {
      global_weight_sums[i + 1] += global_weight_sums[i];
    }
    if (rank > 0) {
      weight_sum = global_weight_sums[rank];
      for (k = 0; k <= local_num_quadrants; ++k) {
        local_weights[k] += weight_sum;
      }
    }
    P4EST_ASSERT (local_weights[0] == global_weight_sums[rank]);
    P4EST_ASSERT (local_weights[local_num_quadrants] ==
                  global_weight_sums[rank + 1]);
    weight_sum = global_weight_sums[num_procs];
    if (weight_sum == 0) {
      /* all quadrants have zero weight, we do nothing */
      P4EST_FREE (local_weights);
      P4EST_FREE (global_weight_sums);
      P4EST_FREE (num_quadrants_in_proc);
      return;
    }

#if(0)
    if (rank == 0) {
      for (i = 0; i <= num_procs; ++i) {
        printf ("global weight sum [%d] %lld\n", i, global_weight_sums[i]);
      }
    }
#endif

    /* determine processor ids to send to */
    send_lowest = num_procs;
    send_highest = 0;
    for (i = 1; i <= num_procs; ++i) {
      cut = (weight_sum * i) / num_procs;
      if (global_weight_sums[rank] < cut &&
          cut <= global_weight_sums[rank + 1]) {
        send_lowest = P4EST_MIN (send_lowest, i);
        send_highest = P4EST_MAX (send_highest, i);
      }
    }
    /*
     * send low cut to send_lowest..send_highest
     * and high cut to send_lowest-1..send_highest-1
     */
#if(0)
    printf ("[%d] send bounds low %d high %d\n",
            rank, send_lowest, send_highest);
#endif
    num_sends = 2 * (send_highest - send_lowest + 1);
    if (num_sends <= 0) {
      num_sends = 0;
      send_requests = NULL;
    }
    else {
      send_requests = P4EST_ALLOC (MPI_Request, num_sends);
      P4EST_CHECK_ALLOC (send_requests);
      k = 0;
      for (i = send_lowest; i <= send_highest; ++i) {
        /* do binary search in the weight array */
        k = int64_find_lower_bound ((weight_sum * i) / num_procs,
                                    local_weights, local_num_quadrants + 1,
                                    k);
        P4EST_ASSERT (k > 0 && k <= local_num_quadrants);
        send_index = k +
          ((rank > 0) ? (p4est->global_last_quad_index[rank - 1] + 1) : -1);

        /* and send the quadrant index as high and low bounds */
        mpiret = MPI_Isend (&send_index, 1, MPI_LONG_LONG, i - 1,
                            P4EST_COMM_PARTITION_WEIGHTED_HIGH,
                            p4est->mpicomm,
                            &send_requests[2 * (i - send_lowest)]);
        P4EST_CHECK_MPI (mpiret);
        if (i < num_procs) {
          mpiret = MPI_Isend (&send_index, 1, MPI_LONG_LONG, i,
                              P4EST_COMM_PARTITION_WEIGHTED_LOW,
                              p4est->mpicomm,
                              &send_requests[2 * (i - send_lowest) + 1]);
        }
        else {
          send_requests[2 * (i - send_lowest) + 1] = MPI_REQUEST_NULL;
        }
        P4EST_CHECK_MPI (mpiret);
      }
    }

    /* determine processor ids to receive from and post irecv */
    i = 0;
    my_lowcut = (weight_sum * rank) / num_procs;
    if (my_lowcut == 0) {
      recv_low = 0;
      recv_requests[0] = MPI_REQUEST_NULL;
    }
    else {
      for (; i < num_procs; ++i) {
        if (global_weight_sums[i] < my_lowcut &&
            my_lowcut <= global_weight_sums[i + 1]) {
          mpiret = MPI_Irecv (&recv_low, 1, MPI_LONG_LONG, i,
                              P4EST_COMM_PARTITION_WEIGHTED_LOW,
                              p4est->mpicomm, &recv_requests[0]);
          P4EST_CHECK_MPI (mpiret);
          break;
        }
      }
    }
    my_highcut = (weight_sum * (rank + 1)) / num_procs;
    if (my_highcut == 0) {
      recv_high = 0;
      recv_requests[1] = MPI_REQUEST_NULL;
    }
    else {
      for (; i < num_procs; ++i) {
        if (global_weight_sums[i] < my_highcut &&
            my_highcut <= global_weight_sums[i + 1]) {
          mpiret = MPI_Irecv (&recv_high, 1, MPI_LONG_LONG, i,
                              P4EST_COMM_PARTITION_WEIGHTED_HIGH,
                              p4est->mpicomm, &recv_requests[1]);
          P4EST_CHECK_MPI (mpiret);
          break;
        }
      }
    }

    /* free temporary memory */
    P4EST_FREE (local_weights);
    P4EST_FREE (global_weight_sums);

    /* wait for sends and receives to complete */
    if (num_sends > 0) {
      mpiret = MPI_Waitall (num_sends, send_requests, MPI_STATUSES_IGNORE);
      P4EST_CHECK_MPI (mpiret);
      P4EST_FREE (send_requests);
    }
    mpiret = MPI_Waitall (2, recv_requests, MPI_STATUSES_IGNORE);
    P4EST_CHECK_MPI (mpiret);

    /* and allgather the quadrant ranges */
    qcount = recv_high - recv_low;
    if (p4est->nout != NULL) {
      fprintf (p4est->nout, "[%d] weighted partition count %d\n",
               rank, qcount);
    }
    mpiret = MPI_Allgather (&qcount, 1, MPI_INT,
                            num_quadrants_in_proc, 1, MPI_INT,
                            p4est->mpicomm);
    P4EST_CHECK_MPI (mpiret);
  }

  p4est_partition_given (p4est, num_quadrants_in_proc);

  P4EST_FREE (num_quadrants_in_proc);
#endif /* HAVE_MPI */
}

unsigned
p4est_checksum (p4est_t * p4est)
{
  unsigned            treecrc, crc;
  int32_t             j;
  p4est_tree_t       *tree;
  p4est_array_t       checkarray;
#ifdef HAVE_MPI
  int                 mpiret;
  uint32_t            send[2];
  uint32_t           *gather;
#endif

  P4EST_ASSERT (p4est_is_valid (p4est));

  p4est_array_init (&checkarray, 4);
  crc = 0;
  for (j = p4est->first_local_tree; j <= p4est->last_local_tree; ++j) {
    tree = p4est_array_index (p4est->trees, j);
    treecrc = p4est_quadrant_checksum (tree->quadrants, &checkarray, 0);
    if (j == p4est->first_local_tree) {
      crc = treecrc;
    }
    else {
      crc = adler32_combine (crc, treecrc, checkarray.elem_count * 4);
    }
  }
  p4est_array_reset (&checkarray);

#ifdef HAVE_MPI
  if (p4est->mpirank == 0) {
    gather = P4EST_ALLOC (uint32_t, 2 * p4est->mpisize);
    P4EST_CHECK_ALLOC (gather);
  }
  else {
    gather = NULL;
  }
  send[0] = crc;
  send[1] = p4est->local_num_quadrants * 12;
  mpiret = MPI_Gather (send, 2, MPI_UNSIGNED, gather, 2, MPI_UNSIGNED,
                       0, p4est->mpicomm);
  P4EST_CHECK_MPI (mpiret);

  crc = 0;
  if (p4est->mpirank == 0) {
    crc = gather[0];
    for (j = 1; j < p4est->mpisize; ++j) {
      crc = adler32_combine (crc, gather[2 * j + 0], gather[2 * j + 1]);
    }
    P4EST_FREE (gather);
  }
#endif

  return crc;
}

/* EOF p4est.c */
