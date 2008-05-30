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

#include <p4est_algorithms.h>
#include <p4est_communication.h>
#include <sc_zlib.h>

typedef struct
{
  int8_t              have_first_count, have_first_load;
  int8_t              have_second_count, have_second_load;
  int                 recv_first_count, recv_second_count;
  int                 send_first_count, send_second_count;
  sc_array_t          send_first, send_second, recv_both;
}
p4est_balance_peer_t;

const int           p4est_corner_to_zorder[5] = { 0, 1, 3, 2, 4 };

static const p4est_gloidx_t initial_quadrants_per_processor = 15;
static const int    number_toread_quadrants = 32;
#ifdef P4EST_MPI
static const int    number_peer_windows = 5;
#endif

static const int8_t fully_owned_flag = 0x01;
static const int8_t any_face_flag = 0x02;

p4est_t            *
p4est_new (MPI_Comm mpicomm, p4est_connectivity_t * connectivity,
           size_t data_size, p4est_init_t init_fn)
{
#ifdef P4EST_MPI
  int                 mpiret;
#endif
  int                 num_procs, rank;
  int                 i, must_remove_last_quadrant;
  int                 level;
  p4est_topidx_t      jl, num_trees;
  p4est_gloidx_t      tree_num_quadrants, global_num_quadrants;
  p4est_gloidx_t      first_tree, first_quadrant, first_tree_quadrant;
  p4est_gloidx_t      last_tree, last_quadrant, last_tree_quadrant;
  p4est_gloidx_t      quadrant_index;
  p4est_t            *p4est;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *quad;
  p4est_quadrant_t    a, b, c;
  p4est_position_t   *global_first_position;

  P4EST_GLOBAL_PRODUCTION ("Into p4est_new\n");
  P4EST_ASSERT (p4est_connectivity_is_valid (connectivity));

  P4EST_QUADRANT_INIT (&a);
  P4EST_QUADRANT_INIT (&b);
  P4EST_QUADRANT_INIT (&c);

  p4est = P4EST_ALLOC_ZERO (p4est_t, 1);

  /* assign some data members */
  p4est->data_size = data_size;
  p4est->user_global_pointer = NULL;
  p4est->connectivity = connectivity;
  num_trees = connectivity->num_trees;

  p4est->mpicomm = mpicomm;
  p4est->mpisize = 1;
  p4est->mpirank = 0;
#ifdef P4EST_MPI
  if (p4est->mpicomm != MPI_COMM_NULL) {
    mpiret = MPI_Comm_size (p4est->mpicomm, &p4est->mpisize);
    SC_CHECK_MPI (mpiret);

    mpiret = MPI_Comm_rank (p4est->mpicomm, &p4est->mpirank);
    SC_CHECK_MPI (mpiret);
  }
#endif
  num_procs = p4est->mpisize;
  rank = p4est->mpirank;

  /* allocate memory pools */
  if (p4est->data_size > 0) {
    p4est->user_data_pool = sc_mempool_new (p4est->data_size);
  }
  else {
    p4est->user_data_pool = NULL;
  }
  p4est->quadrant_pool = sc_mempool_new (sizeof (p4est_quadrant_t));

  /* determine uniform level of initial tree */
  tree_num_quadrants = 1;
  for (level = 0; level < 16; ++level) {
    if (tree_num_quadrants >=
        (num_procs * initial_quadrants_per_processor) / num_trees) {
      break;
    }
    tree_num_quadrants *= 4;
  }
  P4EST_ASSERT (level < 16 && tree_num_quadrants <= INT32_MAX);

  /* compute index of first tree for this processor */
  global_num_quadrants = tree_num_quadrants * num_trees;
  first_quadrant = (global_num_quadrants * rank) / num_procs;
  first_tree = first_quadrant / tree_num_quadrants;
  first_tree_quadrant = first_quadrant - first_tree * tree_num_quadrants;
  last_quadrant = (global_num_quadrants * (rank + 1)) / num_procs - 1;
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
  P4EST_GLOBAL_PRODUCTIONF ("New p4est with %d trees on %d processors\n",
                            num_trees, num_procs);
  P4EST_GLOBAL_INFOF ("Initial level %d potential global quadrants"
                      " %lld per tree %lld\n",
                      level, (long long) global_num_quadrants,
                      (long long) tree_num_quadrants);
  P4EST_VERBOSEF
    ("first tree %lld first quadrant %lld global quadrant %lld\n",
     (long long) first_tree, (long long) first_tree_quadrant,
     (long long) first_quadrant);
  P4EST_VERBOSEF ("last tree %lld last quadrant %lld global quadrant %lld\n",
                  (long long) last_tree, (long long) last_tree_quadrant,
                  (long long) last_quadrant);

  /* allocate trees and quadrants */
  p4est->trees = sc_array_new (sizeof (p4est_tree_t));
  sc_array_resize (p4est->trees, num_trees);
  for (jl = 0; jl < num_trees; ++jl) {
    tree = sc_array_index (p4est->trees, jl);
    sc_array_init (&tree->quadrants, sizeof (p4est_quadrant_t));
    for (i = 0; i <= P4EST_MAXLEVEL; ++i) {
      tree->quadrants_per_level[i] = 0;
    }
    tree->maxlevel = 0;
  }
  p4est->local_num_quadrants = 0;
  p4est->global_num_quadrants = 0;

  /* for every locally non-empty tree fill first and last quadrant */
  for (jl = first_tree; jl <= last_tree; ++jl) {
    tree = sc_array_index (p4est->trees, jl);
    must_remove_last_quadrant = 0;

    /* set morton id of first quadrant and initialize user data */
    if (jl == first_tree) {
      p4est_quadrant_set_morton (&a, level, first_tree_quadrant);
    }
    else {
      p4est_quadrant_set_morton (&a, level, 0);
    }
    P4EST_LDEBUGF ("tree %lld first morton 0x%llx 0x%llx\n",
                   (long long) jl, (long long) a.x, (long long) a.y);
    p4est_quadrant_init_data (p4est, jl, &a, init_fn);

    /* set morton id of last quadrant */
    if (tree_num_quadrants == 1 ||
        (jl == first_tree && first_tree_quadrant == tree_num_quadrants - 1)) {
      /* There is only a in the tree */
      sc_array_resize (&tree->quadrants, 1);
      quad = sc_array_index (&tree->quadrants, 0);
      *quad = a;
      tree->maxlevel = a.level;
      ++tree->quadrants_per_level[a.level];
    }
    else {
      if (jl == last_tree) {
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
      P4EST_LDEBUGF ("tree %lld last morton 0x%llx 0x%llx\n",
                     (long long) jl, (long long) b.x, (long long) b.y);
      if (!must_remove_last_quadrant) {
        p4est_quadrant_init_data (p4est, jl, &b, init_fn);
      }

      /* now run algorithm CompleteRegion (&tree->quadrants) here */
      p4est_complete_region (p4est, &a, 1, &b, !must_remove_last_quadrant,
                             tree, jl, init_fn);
    }
    P4EST_INFOF ("tree %lld quadrants %lu\n",
                 (long long) jl, (unsigned long) tree->quadrants.elem_count);
    p4est->local_num_quadrants += tree->quadrants.elem_count;
  }

  /* compute some member variables */
  p4est->first_local_tree = first_tree;
  p4est->last_local_tree = last_tree;
  p4est->global_last_quad_index = P4EST_ALLOC (p4est_gloidx_t, num_procs);
  p4est_comm_count_quadrants (p4est);

  /* fill in global partition information */
  global_first_position = P4EST_ALLOC_ZERO (p4est_position_t, num_procs + 1);
  for (i = 0; i <= num_procs; ++i) {
    first_quadrant = (global_num_quadrants * i) / num_procs;
    first_tree = first_quadrant / tree_num_quadrants;
    first_tree_quadrant = first_quadrant - first_tree * tree_num_quadrants;
    p4est_quadrant_set_morton (&c, level, first_tree_quadrant);
    global_first_position[i].which_tree = first_tree;
    global_first_position[i].x = c.x;
    global_first_position[i].y = c.y;
  }
  p4est->global_first_position = global_first_position;

  /* print more statistics */
  P4EST_INFOF ("total local quadrants %d\n", p4est->local_num_quadrants);

  P4EST_ASSERT (p4est_is_valid (p4est));
  P4EST_GLOBAL_PRODUCTIONF ("Done p4est_new with %lld total quadrants\n",
                            (long long) p4est->global_num_quadrants);
  return p4est;
}

void
p4est_destroy (p4est_t * p4est)
{
#ifdef P4EST_DEBUG
  int                 q;
#endif
  int32_t             j;
  p4est_tree_t       *tree;

  for (j = 0; j < p4est->connectivity->num_trees; ++j) {
    tree = sc_array_index (p4est->trees, j);

#ifdef P4EST_DEBUG
    for (q = 0; q < tree->quadrants.elem_count; ++q) {
      p4est_quadrant_t   *quad = sc_array_index (&tree->quadrants, q);
      p4est_quadrant_free_data (p4est, quad);
    }
#endif

    sc_array_reset (&tree->quadrants);
  }
  sc_array_destroy (p4est->trees);

  if (p4est->user_data_pool != NULL) {
    sc_mempool_destroy (p4est->user_data_pool);
  }
  sc_mempool_destroy (p4est->quadrant_pool);

  P4EST_FREE (p4est->global_first_position);
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
  sc_array_t         *iquadrants, *pquadrants;

  /* create a shallow copy and zero out dependent fields */
  p4est = P4EST_ALLOC (p4est_t, 1);
  memcpy (p4est, input, sizeof (p4est_t));
  p4est->global_last_quad_index = NULL;
  p4est->global_first_position = NULL;
  p4est->trees = NULL;
  p4est->user_data_pool = NULL;
  p4est->quadrant_pool = NULL;

  /* allocate a user data pool if necessary and a quadrant pool */
  if (copy_data && p4est->data_size > 0) {
    p4est->user_data_pool = sc_mempool_new (p4est->data_size);
  }
  else {
    p4est->data_size = 0;
  }
  p4est->quadrant_pool = sc_mempool_new (sizeof (p4est_quadrant_t));

  /* copy quadrants for each tree */
  p4est->trees = sc_array_new (sizeof (p4est_tree_t));
  sc_array_resize (p4est->trees, num_trees);
  for (j = 0; j < num_trees; ++j) {
    itree = sc_array_index (input->trees, j);
    ptree = sc_array_index (p4est->trees, j);
    memcpy (ptree, itree, sizeof (p4est_tree_t));
    sc_array_init (&ptree->quadrants, sizeof (p4est_quadrant_t));
  }
  for (j = first_tree; j <= last_tree; ++j) {
    itree = sc_array_index (input->trees, j);
    iquadrants = &itree->quadrants;
    icount = iquadrants->elem_count;
    ptree = sc_array_index (p4est->trees, j);
    pquadrants = &ptree->quadrants;
    sc_array_resize (pquadrants, icount);
    memcpy (pquadrants->array, iquadrants->array,
            icount * sizeof (p4est_quadrant_t));
    if (p4est->data_size > 0) {
      for (k = 0; k < icount; ++k) {
        iq = sc_array_index (iquadrants, k);
        pq = sc_array_index (pquadrants, k);
        pq->p.user_data = sc_mempool_alloc (p4est->user_data_pool);
        memcpy (pq->p.user_data, iq->p.user_data, p4est->data_size);
      }
    }
    else {
      for (k = 0; k < icount; ++k) {
        pq = sc_array_index (pquadrants, k);
        pq->p.user_data = NULL;
      }
    }
  }

  /* allocate and copy global quadrant count */
  p4est->global_last_quad_index =
    P4EST_ALLOC (p4est_gloidx_t, p4est->mpisize);
  memcpy (p4est->global_last_quad_index, input->global_last_quad_index,
          p4est->mpisize * sizeof (p4est_gloidx_t));

  /* allocate and copy global partition information */
  p4est->global_first_position = P4EST_ALLOC (p4est_position_t,
                                              p4est->mpisize + 1);
  memcpy (p4est->global_first_position, input->global_first_position,
          (p4est->mpisize + 1) * sizeof (p4est_position_t));

  /* check for valid p4est and return */
  P4EST_ASSERT (p4est_is_valid (p4est));

  return p4est;
}

void
p4est_refine (p4est_t * p4est, p4est_refine_t refine_fn, p4est_init_t init_fn)
{
  size_t              quadrant_pool_size, data_pool_size = 0;
  int                 dorefine;
  int                *key;
  int                 i, maxlevel;
  int32_t             j, movecount;
  int32_t             current, restpos, incount;
  sc_list_t          *list;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *q, *qalloc, *qpop;
  p4est_quadrant_t   *c0, *c1, *c2, *c3;
  sc_array_t         *tquadrants;

  P4EST_GLOBAL_PRODUCTIONF ("Into p4est_refine with %lld total quadrants\n",
                            (long long) p4est->global_num_quadrants);
  P4EST_ASSERT (p4est_is_valid (p4est));

  /*
     q points to a quadrant that is an array member
     qalloc is a quadrant that has been allocated through quadrant_pool
     qpop is a quadrant that has been allocated through quadrant_pool
     never mix these two types of quadrant pointers
   */
  key = &dorefine;              /* use this to create a unique user_data pointer */
  list = sc_list_new (NULL);
  p4est->local_num_quadrants = 0;

  /* loop over all local trees */
  for (j = p4est->first_local_tree; j <= p4est->last_local_tree; ++j) {
    tree = sc_array_index (p4est->trees, j);
    tquadrants = &tree->quadrants;
    quadrant_pool_size = p4est->quadrant_pool->elem_count;
    if (p4est->user_data_pool != NULL) {
      data_pool_size = p4est->user_data_pool->elem_count;
    }

    /* initial log message for this tree */
    P4EST_INFOF ("Into refine tree %d with %lu\n",
                 j, (unsigned long) tquadrants->elem_count);

    /* reset the quadrant counters */
    maxlevel = 0;
    for (i = 0; i <= P4EST_MAXLEVEL; ++i) {
      tree->quadrants_per_level[i] = 0;
    }

    /* run through the array to find first quadrant to be refined */
    q = NULL;
    dorefine = 0;
    incount = tquadrants->elem_count;
    for (current = 0; current < incount; ++current) {
      q = sc_array_index (tquadrants, current);
      dorefine = ((q->level < P4EST_MAXLEVEL) && refine_fn (p4est, j, q));
      if (dorefine) {
        break;
      }
      maxlevel = SC_MAX (maxlevel, q->level);
      ++tree->quadrants_per_level[q->level];
    }
    if (!dorefine) {
      p4est->local_num_quadrants += incount;
      continue;
    }

    /* now we have a quadrant to refine, prepend it to the list */
    qalloc = sc_mempool_alloc (p4est->quadrant_pool);
    *qalloc = *q;               /* never prepend array members */
    sc_list_prepend (list, qalloc);     /* only newly allocated quadrants */

    /*
       current points to the next array member to write
       restpos points to the next array member to read
     */
    restpos = current + 1;

    /* run through the list and refine recursively */
    while (list->elem_count > 0) {
      qpop = sc_list_pop (list);
      if (dorefine ||
          ((qpop->level < P4EST_MAXLEVEL) && refine_fn (p4est, j, qpop))) {
        dorefine = 0;           /* a marker so that refine_fn is never called twice */
        sc_array_resize (tquadrants, tquadrants->elem_count + 3);

        /* compute children and prepend them to the list */
        if (qpop->p.user_data != key) {
          p4est_quadrant_free_data (p4est, qpop);
        }
        c0 = qpop;
        c1 = sc_mempool_alloc (p4est->quadrant_pool);
        c2 = sc_mempool_alloc (p4est->quadrant_pool);
        c3 = sc_mempool_alloc (p4est->quadrant_pool);
        p4est_quadrant_children (qpop, c0, c1, c2, c3);
        c0->p.user_data = key;
        c1->p.user_data = key;
        c2->p.user_data = key;
        c3->p.user_data = key;
        sc_list_prepend (list, c3);
        sc_list_prepend (list, c2);
        sc_list_prepend (list, c1);
        sc_list_prepend (list, c0);
      }
      else {
        /* need to make room in the array to store this new quadrant */
        if (restpos < incount && current == restpos) {
          movecount = SC_MIN (incount - restpos, number_toread_quadrants);
          while (movecount > 0) {
            q = sc_array_index (tquadrants, restpos);
            qalloc = sc_mempool_alloc (p4est->quadrant_pool);
            *qalloc = *q;       /* never append array members */
            sc_list_append (list, qalloc);      /* only newly allocated quadrants */
            --movecount;
            ++restpos;
          }
        }

        /* store new quadrant and update counters */
        if (qpop->p.user_data == key) {
          p4est_quadrant_init_data (p4est, j, qpop, init_fn);
        }
        q = sc_array_index (tquadrants, current);
        *q = *qpop;
        maxlevel = SC_MAX (maxlevel, qpop->level);
        ++tree->quadrants_per_level[qpop->level];
        ++current;
        sc_mempool_free (p4est->quadrant_pool, qpop);
      }
    }
    tree->maxlevel = (int8_t) maxlevel;
    p4est->local_num_quadrants += tquadrants->elem_count;

    P4EST_ASSERT (restpos == incount);
    P4EST_ASSERT (current == tquadrants->elem_count);
    P4EST_ASSERT (list->first == NULL && list->last == NULL);
    P4EST_ASSERT (quadrant_pool_size == p4est->quadrant_pool->elem_count);
    if (p4est->user_data_pool != NULL) {
      P4EST_ASSERT (data_pool_size + tquadrants->elem_count ==
                    p4est->user_data_pool->elem_count + incount);
    }
    P4EST_ASSERT (p4est_tree_is_sorted (tree));
    P4EST_ASSERT (p4est_tree_is_complete (tree));

    /* final log message for this tree */
    P4EST_INFOF ("Done refine tree %d now %lu\n",
                 j, (unsigned long) tquadrants->elem_count);
  }

  sc_list_destroy (list);

  /* compute global number of quadrants */
  p4est_comm_count_quadrants (p4est);

  P4EST_ASSERT (p4est_is_valid (p4est));
  P4EST_GLOBAL_PRODUCTIONF ("Done p4est_refine with %lld total quadrants\n",
                            (long long) p4est->global_num_quadrants);
}

void
p4est_coarsen (p4est_t * p4est, p4est_coarsen_t coarsen_fn,
               p4est_init_t init_fn)
{
  int                 k, couldbegood;
  size_t              data_pool_size = 0;
  int                 incount, removed, num_quadrants;
  int                 first, last, rest, before;
  int                 i, maxlevel;
  int32_t             j;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *c[4];
  p4est_quadrant_t   *cfirst, *clast;
  sc_array_t         *tquadrants;

  P4EST_GLOBAL_PRODUCTIONF ("Into p4est_coarsen with %lld total quadrants\n",
                            (long long) p4est->global_num_quadrants);
  P4EST_ASSERT (p4est_is_valid (p4est));

  /* loop over all local trees */
  for (j = p4est->first_local_tree; j <= p4est->last_local_tree; ++j) {
    tree = sc_array_index (p4est->trees, j);
    tquadrants = &tree->quadrants;
    if (p4est->user_data_pool != NULL) {
      data_pool_size = p4est->user_data_pool->elem_count;
    }
    removed = 0;

    /* initial log message for this tree */
    P4EST_INFOF ("Into coarsen tree %d with %lu\n",
                 j, (unsigned long) tquadrants->elem_count);

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
    incount = tquadrants->elem_count;
    while (rest + 3 - before < incount) {
      couldbegood = 1;
      for (k = 0; k < 4; ++k) {
        if (k < before) {
          c[k] = sc_array_index (tquadrants, first + k);
          if (k != p4est_quadrant_child_id (c[k])) {
            couldbegood = 0;
            break;
          }
        }
        else {
          c[k] = sc_array_index (tquadrants, rest + k - before);
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
        first = SC_MAX (first, 0);
      }
      else {
        /* do nothing, just move the counters and the hole */
        ++first;
        if (first > last) {
          if (first != rest) {
            cfirst = sc_array_index (tquadrants, first);
            clast = sc_array_index (tquadrants, rest);
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
        cfirst = sc_array_index (tquadrants, first);
        clast = sc_array_index (tquadrants, rest);
        *cfirst = *clast;
        ++rest;
      }
      sc_array_resize (tquadrants, first + 1);
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
    tree->maxlevel = (int8_t) maxlevel;

    /* do some sanity checks */
    P4EST_ASSERT (num_quadrants == tquadrants->elem_count);
    P4EST_ASSERT (tquadrants->elem_count == incount - removed);
    if (p4est->user_data_pool != NULL) {
      P4EST_ASSERT (data_pool_size - removed ==
                    p4est->user_data_pool->elem_count);
    }
    P4EST_ASSERT (p4est_tree_is_sorted (tree));
    P4EST_ASSERT (p4est_tree_is_complete (tree));

    /* final log message for this tree */
    P4EST_INFOF ("Done coarsen tree %d now %lu\n",
                 j, (unsigned long) tquadrants->elem_count);
  }

  /* compute global number of quadrants */
  p4est_comm_count_quadrants (p4est);

  P4EST_ASSERT (p4est_is_valid (p4est));
  P4EST_GLOBAL_PRODUCTIONF ("Done p4est_coarsen with %lld total quadrants\n",
                            (long long) p4est->global_num_quadrants);
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
p4est_balance_schedule (p4est_t * p4est, sc_array_t * peers,
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
    peer = sc_array_index (peers, owner);
    /* avoid duplicates in the send array */
    found = 0;
    for (back = 0; back < 8; ++back) {
      pos = peer->send_first.elem_count - back - 1;
      if (pos < 0) {
        break;
      }
      s = sc_array_index (&peer->send_first, pos);
      if (p4est_quadrant_is_equal (s, q) && s->p.piggy.which_tree == qtree) {
        found = 1;
        break;
      }
    }
    if (found) {
      continue;
    }

    /* copy quadrant into shipping list */
    scount = peer->send_first.elem_count;
    sc_array_resize (&peer->send_first, scount + 1);
    s = sc_array_index (&peer->send_first, scount);
    *s = *q;
    s->p.piggy.which_tree = qtree;      /* piggy back tree id */

    /* update lowest and highest peer */
    *first_peer = SC_MIN (owner, *first_peer);
    *last_peer = SC_MAX (owner, *last_peer);
  }
}

static void
p4est_balance_response (p4est_t * p4est, int32_t peer_id,
                        p4est_balance_peer_t * peer)
{
#ifdef P4EST_DEBUG
  const int32_t       first_tree = p4est->first_local_tree;
  const int32_t       last_tree = p4est->last_local_tree;
#endif /* P4EST_DEBUG */
  int32_t             k, qcount, qtree;
  int32_t             prev, num_receive_trees, nt;
  int32_t            *pi;
  sc_array_t         *qarray, tree_array;
  p4est_quadrant_t   *q;

  qarray = &peer->recv_both;
  qcount = qarray->elem_count;
  P4EST_ASSERT (peer->recv_first_count == qcount);

  /* build list of received tree ids */
  prev = -1;
  num_receive_trees = 0;
  sc_array_init (&tree_array, sizeof (int32_t));
  for (k = 0; k < qcount; ++k) {
    q = sc_array_index (qarray, k);
    qtree = q->p.piggy.which_tree;
    P4EST_ASSERT (first_tree <= qtree && qtree <= last_tree);
    P4EST_ASSERT (qtree >= prev);
    if (qtree > prev) {
      sc_array_resize (&tree_array, num_receive_trees + 1);
      pi = sc_array_index (&tree_array, num_receive_trees);
      *pi = qtree;
      ++num_receive_trees;
      prev = qtree;
    }
  }
  P4EST_LDEBUGF ("first load from %d into %d trees\n", peer_id,
                 num_receive_trees);
  P4EST_ASSERT (num_receive_trees == tree_array.elem_count);

  /* loop to the trees to receive into and update overlap quadrants */
  for (nt = 0; nt < num_receive_trees; ++nt) {
    pi = sc_array_index (&tree_array, nt);
    qtree = *pi;

    /* compute overlap quadrants */
    p4est_tree_compute_overlap (p4est, qtree, qarray, &peer->send_second);
  }
  p4est_tree_uniqify_overlap (&peer->send_first, &peer->send_second);
  sc_array_reset (&tree_array);
}

void
p4est_balance (p4est_t * p4est, p4est_init_t init_fn)
{
  const int           rank = p4est->mpirank;
  const int           num_procs = p4est->mpisize;
  size_t              data_pool_size = 0;
  size_t              all_incount, all_outcount;
  int                 k, l, ctree;
  int                 any_face, face_contact[4];
  int                 any_quad, quad_contact[4];
  int                 tree_fully_owned, transform;
  int                 first_index, last_index;
  int                 which;
  int                 face, corner, zcorner;
#ifdef P4EST_MPI
  int                 rcount;
#endif
  int8_t             *tree_flags;
  int32_t             i, j;
  int32_t             qtree;
  p4est_qcoord_t      qh;
  const p4est_qcoord_t rh = P4EST_ROOT_LEN;
  int32_t             treecount, qcount, qbytes, offset, obytes;
  int32_t             first_tree, last_tree, next_tree;
  int32_t             first_peer, last_peer;
  int32_t             over_peer_count;
  sc_array_t         *peers, *qarray, *tquadrants, corner_info;
  p4est_balance_peer_t *peer;
  p4est_tree_t       *tree;
  p4est_quadrant_t    mylow, nextlow;
  p4est_quadrant_t    tosend, insulq, tempq;
  p4est_quadrant_t   *q, *s;
  p4est_connectivity_t *conn = p4est->connectivity;
  p4est_corner_info_t *ci;
#ifdef P4EST_MPI
#ifdef P4EST_DEBUG
  unsigned            checksum;
  sc_array_t          checkarray;
  int32_t             ltotal[2], gtotal[2];
#endif /* P4EST_DEBUG */
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
  int32_t             rip, eff_peer_count;
  int32_t             length, shortest_window, shortest_length;
  int32_t             peer_windows[twopeerw];
  int32_t            *peer_boundaries;
  MPI_Request        *requests_first, *requests_second;
  MPI_Request        *send_requests_first_count, *send_requests_first_load;
  MPI_Request        *send_requests_second_count, *send_requests_second_load;
  MPI_Status         *recv_statuses, *jstatus;
#endif /* P4EST_MPI */

  P4EST_GLOBAL_PRODUCTIONF ("Into p4est_balance with %lld total quadrants\n",
                            (long long) p4est->global_num_quadrants);
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
  for (i = 0; i < conn->num_trees; ++i) {
    tree_flags[i] = 0x00;
  }

#ifdef P4EST_MPI
  /* will contain first and last peer (inclusive) for each processor */
  peer_boundaries = P4EST_ALLOC (int32_t, twopeerw * num_procs);
  /* request and status buffers for receive operations */
  requests_first = P4EST_ALLOC (MPI_Request, 6 * num_procs);
  requests_second = requests_first + 1 * num_procs;
  send_requests_first_count = requests_first + 2 * num_procs;
  send_requests_first_load = requests_first + 3 * num_procs;
  send_requests_second_count = requests_first + 4 * num_procs;
  send_requests_second_load = requests_first + 5 * num_procs;
  recv_statuses = P4EST_ALLOC (MPI_Status, num_procs);
  for (i = 0; i < num_procs; ++i) {
    requests_first[i] = requests_second[i] = MPI_REQUEST_NULL;
    send_requests_first_count[i] = MPI_REQUEST_NULL;
    send_requests_first_load[i] = MPI_REQUEST_NULL;
    send_requests_second_count[i] = MPI_REQUEST_NULL;
    send_requests_second_load[i] = MPI_REQUEST_NULL;
  }
  /* index buffer for call to waitsome */
  wait_indices = P4EST_ALLOC (int, num_procs);
#ifdef P4EST_DEBUG
  sc_array_init (&checkarray, 4);
#endif /* P4EST_DEBUG */
#endif /* P4EST_MPI */

  /* allocate per peer storage and initialize requests */
  peers = sc_array_new (sizeof (p4est_balance_peer_t));
  sc_array_resize (peers, num_procs);
  for (i = 0; i < num_procs; ++i) {
    peer = sc_array_index (peers, i);
    sc_array_init (&peer->send_first, sizeof (p4est_quadrant_t));
    sc_array_init (&peer->send_second, sizeof (p4est_quadrant_t));
    sc_array_init (&peer->recv_both, sizeof (p4est_quadrant_t));
    peer->send_first_count = peer->send_second_count = 0;
    peer->recv_first_count = peer->recv_second_count = 0;
    peer->have_first_count = peer->have_first_load = 0;
    peer->have_second_count = peer->have_second_load = 0;
  }
  sc_array_init (&corner_info, sizeof (p4est_corner_info_t));

  /* compute first quadrants on finest level for comparison for me and next */
  first_peer = num_procs;
  last_peer = -1;
  first_tree = p4est->first_local_tree;
  last_tree = p4est->last_local_tree;
  if (first_tree < 0) {
    P4EST_ASSERT (first_tree == -1 && last_tree == -2);
  }
  else {
    P4EST_ASSERT (p4est->global_first_position[rank].which_tree ==
                  first_tree);
  }
  mylow.x = p4est->global_first_position[rank].x;
  mylow.y = p4est->global_first_position[rank].y;
  mylow.level = P4EST_MAXLEVEL;
  next_tree = p4est->global_first_position[rank + 1].which_tree;
  if (last_tree < 0) {
    P4EST_ASSERT (first_tree == -1 && last_tree == -2);
  }
  else {
    P4EST_ASSERT (next_tree == last_tree || next_tree == last_tree + 1);
  }
  nextlow.x = p4est->global_first_position[rank + 1].x;
  nextlow.y = p4est->global_first_position[rank + 1].y;
  nextlow.level = P4EST_MAXLEVEL;

  /* loop over all local trees to assemble first send list */
  all_incount = 0;
  for (j = first_tree; j <= last_tree; ++j) {
    any_face = 0;
    for (face = 0; face < 4; ++face) {
      face_contact[face] = (conn->tree_to_tree[4 * j + face] != j ||
                            conn->tree_to_face[4 * j + face] != face);
      any_face = any_face || face_contact[face];
    }
    if (any_face) {
      tree_flags[j] |= any_face_flag;
    }
    tree = sc_array_index (p4est->trees, j);
    tquadrants = &tree->quadrants;
    all_incount += tquadrants->elem_count;

    /* initial log message for this tree */
    P4EST_INFOF ("Into balance tree %d with %lu\n",
                 j, (unsigned long) tquadrants->elem_count);

    /* local balance first pass */
    p4est_balance_subtree (p4est, tree, j, init_fn);
    treecount = tquadrants->elem_count;
    P4EST_VERBOSEF ("Balance tree %d A %d\n", j, treecount);

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
      q = sc_array_index (tquadrants, i);
      qh = P4EST_QUADRANT_LEN (q->level);
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
                ci = sc_array_index (&corner_info, ctree);
                tosend = *q;
                zcorner = p4est_corner_to_zorder[ci->ncorner];
                p4est_quadrant_corner (&tosend, zcorner, 0);
                p4est_quadrant_corner (&insulq, zcorner, 1);
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
    tquadrants = NULL;          /* safeguard */
  }
  if (last_peer < 0) {
    P4EST_ASSERT (first_peer == num_procs);
    first_peer = last_peer = rank;
  }
  P4EST_ASSERT (first_peer <= last_peer);
  P4EST_ASSERT (0 <= first_peer && last_peer < num_procs);
  over_peer_count = last_peer - first_peer + 1;
  P4EST_ASSERT (0 <= over_peer_count && over_peer_count <= num_procs);
  if (num_procs == 1) {
    P4EST_ASSERT (first_peer == rank && last_peer == rank);
  }

#ifdef P4EST_MPI
  /* compute information about peers to inform receivers */
  rip = (first_peer <= rank && rank <= last_peer) ? 1 : 0;
  for (i = 0; i < twopeerw; ++i) {
    peer_windows[i] = -1;
  }

  /* find a maximum of npw-1 empty ranges with (start, end) */
  lastw = number_peer_windows - 1;
  prev = -1;
  nwin = 0;
  for (j = 0; j < num_procs; ++j) {
    peer = sc_array_index (peers, j);
    if (peer->send_first.elem_count == 0 || j == rank) {
      continue;
    }
    if (prev == -1) {
      prev = j;
      continue;
    }
    if (prev < j - 1) {
      length = j - 1 - prev;
      P4EST_LDEBUGF ("found empty range prev %d j %d length %d\n",
                     prev, j, length);

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
        shortest_length = num_procs + 1;
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

#ifdef P4EST_DEBUG
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
  for (i = 0; i < nwin; ++i) {
    P4EST_LDEBUGF ("empty range %d from %d to %d\n",
                   i, peer_windows[2 * i], peer_windows[2 * i + 1]);
  }
#endif /* P4EST_DEBUG */

  /* compute windows from empty ranges */
  eff_peer_count = 0;
  peer_windows[2 * nwin + 1] = last_peer;
  for (i = nwin; i > 0; --i) {
    peer_windows[2 * i] = peer_windows[2 * i - 1] + 1;
    peer_windows[2 * i - 1] = peer_windows[2 * (i - 1)] - 1;
  }
  peer_windows[0] = first_peer;
  ++nwin;
  for (i = 0; i < nwin; ++i) {
    P4EST_ASSERT (peer_windows[2 * i] <= peer_windows[2 * i + 1]);
    if (i < nwin - 1) {
      P4EST_ASSERT (peer_windows[2 * i + 1] < peer_windows[2 * (i + 1)] - 1);
    }
    eff_peer_count += peer_windows[2 * i + 1] - peer_windows[2 * i] + 1;
  }

#ifdef P4EST_DEBUG
  for (i = nwin; i < number_peer_windows; ++i) {
    P4EST_ASSERT (peer_windows[2 * i] == -1);
    P4EST_ASSERT (peer_windows[2 * i + 1] == -1);
  }
  for (i = 0; i < nwin; ++i) {
    P4EST_LDEBUGF ("window %d from %d to %d\n",
                   i, peer_windows[2 * i], peer_windows[2 * i + 1]);
  }
#endif /* P4EST_DEBUG */

  /* distribute information about peer windows to inform receivers */
  if (p4est->mpicomm != MPI_COMM_NULL) {
    mpiret = MPI_Allgather (peer_windows, twopeerw, MPI_INT,
                            peer_boundaries, twopeerw, MPI_INT,
                            p4est->mpicomm);
    SC_CHECK_MPI (mpiret);
  }
  P4EST_INFOF ("Peers first %d last %d counts %d %d\n",
               first_peer, last_peer, over_peer_count, eff_peer_count);

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
    peer = sc_array_index (peers, j);
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
      P4EST_LDEBUGF ("Balance A send %d quadrants to %d\n", qcount, j);
      ++send_load[0];
    }
    else {
      ++send_zero[0];
    }
    peer->send_first_count = qcount;
    mpiret = MPI_Isend (&peer->send_first_count, 1, MPI_INT,
                        j, P4EST_COMM_BALANCE_FIRST_COUNT,
                        p4est->mpicomm, &send_requests_first_count[j]);
    SC_CHECK_MPI (mpiret);
    ++request_send_count;

    /* sort and send the actual quadrants and post receive for reply */
    if (qcount > 0) {
      sc_array_sort (&peer->send_first, p4est_quadrant_compare_piggy);

#ifdef P4EST_DEBUG
      checksum = p4est_quadrant_checksum (&peer->send_first, &checkarray, 0);
      P4EST_LDEBUGF ("Balance A send checksum %x to %d\n", checksum, j);
#endif /* P4EST_DEBUG */

      total_send_count += qcount;
      qbytes = qcount * sizeof (p4est_quadrant_t);
      mpiret = MPI_Isend (peer->send_first.array, qbytes, MPI_BYTE,
                          j, P4EST_COMM_BALANCE_FIRST_LOAD,
                          p4est->mpicomm, &send_requests_first_load[j]);
      SC_CHECK_MPI (mpiret);
      ++request_send_count;
      mpiret = MPI_Irecv (&peer->recv_second_count, 1, MPI_INT,
                          j, P4EST_COMM_BALANCE_SECOND_COUNT,
                          p4est->mpicomm, &requests_second[j]);
      SC_CHECK_MPI (mpiret);
      ++request_second_count;
    }
  }

  /* find out who is sending to me and receive quadrant counts */
  for (j = 0; j < num_procs; ++j) {
    if (j == rank) {
      continue;
    }
    peer = sc_array_index (peers, j);
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
    peer = sc_array_index (peers, j);
    ++request_first_count;

    /* processor j is sending to me */
    mpiret = MPI_Irecv (&peer->recv_first_count, 1, MPI_INT,
                        j, P4EST_COMM_BALANCE_FIRST_COUNT,
                        p4est->mpicomm, &requests_first[j]);
    SC_CHECK_MPI (mpiret);
  }

  /* wait for quadrant counts and post receive and send for quadrants */
  while (request_first_count > 0) {
    mpiret = MPI_Waitsome (num_procs, requests_first,
                           &outcount, wait_indices, recv_statuses);
    SC_CHECK_MPI (mpiret);
    P4EST_ASSERT (outcount != MPI_UNDEFINED);
    P4EST_ASSERT (outcount > 0);
    for (i = 0; i < outcount; ++i) {
      /* retrieve sender's rank */
      j = wait_indices[i];
      jstatus = &recv_statuses[i];
      wait_indices[i] = -1;
      P4EST_ASSERT (0 <= j && j < num_procs);
      P4EST_ASSERT (requests_first[j] == MPI_REQUEST_NULL);
      P4EST_ASSERT (jstatus->MPI_SOURCE == j);

      /* check if we are in receiving count or load */
      peer = sc_array_index (peers, j);
      P4EST_ASSERT (!peer->have_first_load);
      if (!peer->have_first_count) {
        /* verify message size */
        P4EST_ASSERT (jstatus->MPI_TAG == P4EST_COMM_BALANCE_FIRST_COUNT);
        mpiret = MPI_Get_count (jstatus, MPI_INT, &rcount);
        SC_CHECK_MPI (mpiret);
        SC_CHECK_ABORTF (rcount == 1, "Receive count mismatch A %d", rcount);

        /* process the count information received */
        peer->have_first_count = 1;
        qcount = peer->recv_first_count;
        if (qcount > 0) {
          /* received nonzero count, post receive for load */
          P4EST_LDEBUGF ("Balance A recv %d quadrants from %d\n", qcount, j);
          P4EST_ASSERT (peer->recv_both.elem_count == 0);
          sc_array_resize (&peer->recv_both, qcount);
          total_recv_count += qcount;
          qbytes = qcount * sizeof (p4est_quadrant_t);
          P4EST_ASSERT (requests_first[j] == MPI_REQUEST_NULL);
          mpiret = MPI_Irecv (peer->recv_both.array, qbytes, MPI_BYTE,
                              j, P4EST_COMM_BALANCE_FIRST_LOAD,
                              p4est->mpicomm, &requests_first[j]);
          SC_CHECK_MPI (mpiret);
          ++recv_load[0];
        }
        else {
          /* will not receive load, close this request */
          P4EST_ASSERT (qcount == 0);
          P4EST_ASSERT (requests_first[j] == MPI_REQUEST_NULL);
          --request_first_count;
          ++recv_zero[0];
        }
      }
      else {
        /* verify received size */
        P4EST_ASSERT (jstatus->MPI_TAG == P4EST_COMM_BALANCE_FIRST_LOAD);
        P4EST_ASSERT (peer->recv_first_count > 0);
        mpiret = MPI_Get_count (jstatus, MPI_BYTE, &rcount);
        SC_CHECK_MPI (mpiret);
        SC_CHECK_ABORTF (rcount ==
                         peer->recv_first_count *
                         sizeof (p4est_quadrant_t),
                         "Receive load mismatch A %d %dx%lld", rcount,
                         peer->recv_first_count,
                         (long long) sizeof (p4est_quadrant_t));

        /* received load, close this request */
        peer->have_first_load = 1;
        P4EST_ASSERT (requests_first[j] == MPI_REQUEST_NULL);
        --request_first_count;

#ifdef P4EST_DEBUG
        checksum = p4est_quadrant_checksum (&peer->recv_both, &checkarray, 0);
        P4EST_LDEBUGF ("Balance A recv checksum %x from %d\n", checksum, j);
#endif /* P4EST_DEBUG */

        /* process incoming quadrants to interleave with communication */
        p4est_balance_response (p4est, j, peer);
        qcount = peer->send_second.elem_count;
        if (qcount > 0) {
          P4EST_LDEBUGF ("Balance B send %d quadrants to %d\n", qcount, j);
          ++send_load[1];
        }
        else {
          ++send_zero[1];
        }
        peer->send_second_count = qcount;
        mpiret = MPI_Isend (&peer->send_second_count, 1, MPI_INT,
                            j, P4EST_COMM_BALANCE_SECOND_COUNT,
                            p4est->mpicomm, &send_requests_second_count[j]);
        SC_CHECK_MPI (mpiret);
        ++request_send_count;
        if (qcount > 0) {

#ifdef P4EST_DEBUG
          checksum =
            p4est_quadrant_checksum (&peer->send_second, &checkarray, 0);
          P4EST_LDEBUGF ("Balance B send checksum %x to %d\n", checksum, j);
#endif /* P4EST_DEBUG */

          total_send_count += qcount;
          qbytes = qcount * sizeof (p4est_quadrant_t);
          mpiret = MPI_Isend (peer->send_second.array, qbytes, MPI_BYTE,
                              j, P4EST_COMM_BALANCE_SECOND_LOAD,
                              p4est->mpicomm, &send_requests_second_load[j]);
          SC_CHECK_MPI (mpiret);
          ++request_send_count;
        }
      }
    }
  }
#ifdef P4EST_DEBUG
  for (j = 0; j < num_procs; ++j) {
    SC_CHECK_ABORT (requests_first[j] == MPI_REQUEST_NULL, "Request A");
  }
#endif /* P4EST_DEBUG */
#endif /* P4EST_MPI */

  /* simulate send and receive with myself across tree boundaries */
  peer = sc_array_index (peers, rank);
  sc_array_sort (&peer->send_first, p4est_quadrant_compare_piggy);
  offset = peer->send_first.elem_count;
  peer->recv_first_count = peer->send_first_count = offset;
  obytes = offset * sizeof (p4est_quadrant_t);
  qarray = &peer->recv_both;
  sc_array_resize (qarray, offset);
  memcpy (qarray->array, peer->send_first.array, obytes);
  p4est_balance_response (p4est, rank, peer);
  qcount = peer->send_second.elem_count;
  peer->recv_second_count = peer->send_second_count = qcount;
  qbytes = qcount * sizeof (p4est_quadrant_t);
  sc_array_resize (qarray, offset + qcount);
  memcpy (qarray->array + obytes, peer->send_second.array, qbytes);

#ifdef P4EST_MPI
  /* receive second round appending to the same receive buffer */
  while (request_second_count > 0) {
    mpiret = MPI_Waitsome (num_procs, requests_second,
                           &outcount, wait_indices, recv_statuses);
    SC_CHECK_MPI (mpiret);
    P4EST_ASSERT (outcount != MPI_UNDEFINED);
    P4EST_ASSERT (outcount > 0);
    for (i = 0; i < outcount; ++i) {
      /* retrieve sender's rank */
      j = wait_indices[i];
      jstatus = &recv_statuses[i];
      wait_indices[i] = -1;
      P4EST_ASSERT (0 <= j && j < num_procs);
      P4EST_ASSERT (requests_second[j] == MPI_REQUEST_NULL);
      P4EST_ASSERT (jstatus->MPI_SOURCE == j);

      /* check if we are in receiving count or load */
      peer = sc_array_index (peers, j);
      P4EST_ASSERT (!peer->have_second_load);
      if (!peer->have_second_count) {
        /* verify message size */
        P4EST_ASSERT (jstatus->MPI_TAG == P4EST_COMM_BALANCE_SECOND_COUNT);
        mpiret = MPI_Get_count (jstatus, MPI_INT, &rcount);
        SC_CHECK_MPI (mpiret);
        SC_CHECK_ABORTF (rcount == 1, "Receive count mismatch B %d", rcount);

        /* process the count information received */
        peer->have_second_count = 1;
        qcount = peer->recv_second_count;
        if (qcount > 0) {
          /* received nonzero count, post receive for load */
          P4EST_LDEBUGF ("Balance B recv %d quadrants from %d\n", qcount, j);
          offset = peer->recv_both.elem_count;
          P4EST_ASSERT (offset == peer->recv_first_count);
          obytes = offset * sizeof (p4est_quadrant_t);
          sc_array_resize (&peer->recv_both, offset + qcount);
          total_recv_count += qcount;
          qbytes = qcount * sizeof (p4est_quadrant_t);
          P4EST_ASSERT (requests_second[j] == MPI_REQUEST_NULL);
          mpiret = MPI_Irecv (peer->recv_both.array + obytes, qbytes,
                              MPI_BYTE, j, P4EST_COMM_BALANCE_SECOND_LOAD,
                              p4est->mpicomm, &requests_second[j]);
          SC_CHECK_MPI (mpiret);
          ++recv_load[1];
        }
        else {
          /* will not receive load, close this request */
          P4EST_ASSERT (qcount == 0);
          P4EST_ASSERT (requests_second[j] == MPI_REQUEST_NULL);
          --request_second_count;
          ++recv_zero[1];
        }
      }
      else {
        /* verify received size */
        P4EST_ASSERT (jstatus->MPI_TAG == P4EST_COMM_BALANCE_SECOND_LOAD);
        P4EST_ASSERT (peer->recv_second_count > 0);
        mpiret = MPI_Get_count (jstatus, MPI_BYTE, &rcount);
        SC_CHECK_MPI (mpiret);
        SC_CHECK_ABORTF (rcount ==
                         peer->recv_second_count *
                         sizeof (p4est_quadrant_t),
                         "Receive load mismatch B %d %dx%lld", rcount,
                         peer->recv_second_count,
                         (long long) sizeof (p4est_quadrant_t));

        /* received load, close this request */
        peer->have_second_load = 1;
        P4EST_ASSERT (requests_second[j] == MPI_REQUEST_NULL);
        --request_second_count;

#ifdef P4EST_DEBUG
        checksum = p4est_quadrant_checksum (&peer->recv_both, &checkarray,
                                            peer->recv_first_count);
        P4EST_LDEBUGF ("Balance B recv checksum %x from %d\n", checksum, j);
#endif /* P4EST_DEBUG */
      }
    }
  }
#ifdef P4EST_DEBUG
  for (j = 0; j < num_procs; ++j) {
    SC_CHECK_ABORT (requests_second[j] == MPI_REQUEST_NULL, "Request B");
  }
#endif /* P4EST_DEBUG */

  /* print buffer statistics */
  P4EST_VERBOSEF ("first send Z %d L %d recv Z %d L %d\n",
                  send_zero[0], send_load[0], recv_zero[0], recv_load[0]);
  P4EST_VERBOSEF ("second send Z %d L %d recv Z %d L %d\n",
                  send_zero[1], send_load[1], recv_zero[1], recv_load[1]);
  P4EST_INFOF ("total send %d recv %d\n", total_send_count, total_recv_count);
  for (j = 0; j < num_procs; ++j) {
    peer = sc_array_index (peers, j);
    if (peer->send_first.elem_count > 0 || peer->recv_first_count > 0 ||
        peer->send_second.elem_count > 0 || peer->recv_second_count > 0) {
      P4EST_VERBOSEF ("peer %d first S %lu R %d second S %lu R %d\n",
                      j, (unsigned long) peer->send_first.elem_count,
                      peer->recv_first_count,
                      (unsigned long) peer->send_second.elem_count,
                      peer->recv_second_count);
    }
  }
  if (number_peer_windows == 1) {
    P4EST_ASSERT (send_zero[0] + send_load[0] == over_peer_count - rip);
    P4EST_ASSERT (send_zero[0] + recv_zero[1] + recv_load[1] ==
                  over_peer_count - rip);
  }
#endif /* P4EST_MPI */

  /* merge received quadrants */
  for (j = 0; j < num_procs; ++j) {
    /* access peer information */
    peer = sc_array_index (peers, j);
    qarray = &peer->recv_both;
    qcount = qarray->elem_count;
    P4EST_ASSERT (qcount == peer->recv_first_count + peer->recv_second_count);
    P4EST_ASSERT (peer->send_first_count == peer->send_first.elem_count);
    P4EST_ASSERT (peer->send_second_count == peer->send_second.elem_count);
    if (qcount == 0) {
      continue;
    }

    /* merge received quadrants into correct tree */
    for (k = 0; k < qcount; ++k) {
      s = sc_array_index (qarray, k);
      P4EST_ASSERT (p4est_quadrant_is_extended (s));
      qtree = s->p.piggy.which_tree;
      if (qtree < first_tree || qtree > last_tree) {
        /* this is a corner quadrant from the second pass of balance */
        P4EST_ASSERT (k >= peer->recv_first_count);
        P4EST_ASSERT (0 <= qtree && qtree < conn->num_trees);
        P4EST_ASSERT ((s->x < 0 && s->y < 0) || (s->x < 0 && s->y >= rh) ||
                      (s->x >= rh && s->y < 0) || (s->x >= rh && s->y >= rh));
        continue;
      }
      tree = sc_array_index (p4est->trees, qtree);
      treecount = tree->quadrants.elem_count;
      sc_array_resize (&tree->quadrants, treecount + 1);
      q = sc_array_index (&tree->quadrants, treecount);
      *q = *s;
      ++tree->quadrants_per_level[q->level];
      tree->maxlevel = (int8_t) SC_MAX (tree->maxlevel, q->level);
      ++p4est->local_num_quadrants;
      p4est_quadrant_init_data (p4est, qtree, q, init_fn);
    }
  }

  /* rebalance and clamp result back to original tree boundaries */
  p4est->local_num_quadrants = 0;
  for (j = first_tree; j <= last_tree; ++j) {
    /* check if we are the only processor in an isolated tree */
    tree = sc_array_index (p4est->trees, j);
    tquadrants = &tree->quadrants;
    if ((tree_flags[j] & fully_owned_flag) &&
        !(tree_flags[j] & any_face_flag)) {
      p4est->local_num_quadrants += tquadrants->elem_count;
      continue;
    }

    /* we have most probably received quadrants, run sort and balance */
    sc_array_sort (tquadrants, p4est_quadrant_compare);
    p4est_balance_subtree (p4est, tree, j, init_fn);
    P4EST_VERBOSEF ("Balance tree %d B %lu\n",
                    j, (unsigned long) tquadrants->elem_count);
    treecount = tquadrants->elem_count;

    /* figure out the new elements outside the original tree */
    for (first_index = 0; first_index < treecount; ++first_index) {
      q = sc_array_index (tquadrants, first_index);
      P4EST_ASSERT (p4est_quadrant_is_extended (q));
      if (p4est_quadrant_is_inside (q)) {
        break;
      }
    }
    if (j == first_tree) {
      for (; first_index < treecount; ++first_index) {
        q = sc_array_index (tquadrants, first_index);
        if (p4est_quadrant_compare (q, &mylow) >= 0 ||
            (q->x == mylow.x && q->y == mylow.y)) {
          break;
        }
      }
    }
    for (last_index = treecount - 1; last_index >= 0; --last_index) {
      q = sc_array_index (tquadrants, last_index);
      P4EST_ASSERT (p4est_quadrant_is_extended (q));
      if (p4est_quadrant_is_inside (q)) {
        break;
      }
    }
    if (j == next_tree) {
      for (; last_index >= 0; --last_index) {
        q = sc_array_index (tquadrants, last_index);
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
        q = sc_array_index (tquadrants, k);
        s = sc_array_index (tquadrants, first_index + k);
        if (k < first_index) {
          p4est_quadrant_free_data (p4est, q);
        }
        *q = *s;
        ++k;
      }
      while (k < first_index) {
        q = sc_array_index (tquadrants, k);
        p4est_quadrant_free_data (p4est, q);
        ++k;
      }
    }

    /* remove last part of tree */
    qcount = last_index - first_index + 1;
    for (k = last_index + 1; k < treecount; ++k) {
      q = sc_array_index (tquadrants, k);
      p4est_quadrant_free_data (p4est, q);
    }

    P4EST_ASSERT (qcount <= treecount);
    sc_array_resize (tquadrants, qcount);
    for (l = 0; l <= P4EST_MAXLEVEL; ++l) {
      tree->quadrants_per_level[l] = 0;
    }
    tree->maxlevel = 0;
    for (k = 0; k < qcount; ++k) {
      q = sc_array_index (tquadrants, k);
      P4EST_ASSERT (p4est_quadrant_is_valid (q));
      ++tree->quadrants_per_level[q->level];
      tree->maxlevel = (int8_t) SC_MAX (tree->maxlevel, q->level);
    }
    p4est->local_num_quadrants += qcount;

    P4EST_VERBOSEF ("Balance tree %d C %lu\n",
                    j, (unsigned long) tquadrants->elem_count);
    tquadrants = NULL;          /* safeguard */
  }

#ifdef P4EST_MPI
/* compute global sum of send and receive counts */
#ifdef P4EST_DEBUG
  gtotal[0] = gtotal[1] = 0;
  ltotal[0] = total_send_count;
  ltotal[1] = total_recv_count;
  if (p4est->mpicomm != MPI_COMM_NULL) {
    mpiret = MPI_Reduce (ltotal, gtotal, 2, MPI_INT,
                         MPI_SUM, 0, p4est->mpicomm);
    SC_CHECK_MPI (mpiret);
    P4EST_GLOBAL_STATISTICSF ("Global number of shipped quadrants %d\n",
                              gtotal[0]);
    P4EST_ASSERT (gtotal[0] == gtotal[1]);
  }
  else {
    P4EST_ASSERT (ltotal[0] == 0 && ltotal[1] == 0);
  }
#endif /* P4EST_DEBUG */

  /* wait for all send operations */
  if (request_send_count > 0) {
    mpiret = MPI_Waitall (4 * num_procs,
                          send_requests_first_count, MPI_STATUSES_IGNORE);
    SC_CHECK_MPI (mpiret);
  }
#endif /* P4EST_MPI */

  /* loop over all local trees to finalize balance */
  all_outcount = 0;
  for (j = first_tree; j <= last_tree; ++j) {
    tree = sc_array_index (p4est->trees, j);
    all_outcount += tree->quadrants.elem_count;

    /* final log message for this tree */
    P4EST_INFOF ("Done balance tree %d now %lu\n",
                 j, (unsigned long) tree->quadrants.elem_count);
  }

  /* cleanup temporary storage */
  P4EST_FREE (tree_flags);
  for (i = 0; i < num_procs; ++i) {
    peer = sc_array_index (peers, i);
    sc_array_reset (&peer->send_first);
    sc_array_reset (&peer->send_second);
    sc_array_reset (&peer->recv_both);
  }
  sc_array_destroy (peers);
  sc_array_reset (&corner_info);
#ifdef P4EST_MPI
  P4EST_FREE (peer_boundaries);
  P4EST_FREE (requests_first);  /* includes allocation for requests_second */
  P4EST_FREE (recv_statuses);
  P4EST_FREE (wait_indices);
#ifdef P4EST_DEBUG
  sc_array_reset (&checkarray);
#endif /* P4EST_DEBUG */
#endif /* P4EST_MPI */

  /* compute global number of quadrants */
  p4est_comm_count_quadrants (p4est);

  /* some sanity checks */
  P4EST_ASSERT (all_outcount == p4est->local_num_quadrants);
  P4EST_ASSERT (all_outcount >= all_incount);
  if (p4est->user_data_pool != NULL) {
    P4EST_ASSERT (data_pool_size + (all_outcount - all_incount) ==
                  p4est->user_data_pool->elem_count);
  }
  P4EST_ASSERT (p4est_is_valid (p4est));
  P4EST_GLOBAL_PRODUCTIONF ("Done p4est_balance with %lld total quadrants\n",
                            (long long) p4est->global_num_quadrants);
}

void
p4est_partition (p4est_t * p4est, p4est_weight_t weight_fn)
{
  p4est_gloidx_t      global_shipped = 0;
  const p4est_gloidx_t global_num_quadrants = p4est->global_num_quadrants;
#ifdef P4EST_MPI
  int                 mpiret;
  int                 low_source, high_source;
  const int           num_procs = p4est->mpisize;
  const int           rank = p4est->mpirank;
  const p4est_topidx_t first_tree = p4est->first_local_tree;
  const p4est_topidx_t last_tree = p4est->last_local_tree;
  const p4est_locidx_t local_num_quadrants = p4est->local_num_quadrants;
  int                 i, p;
  int                 send_lowest, send_highest;
  int                 num_sends, rcount, base_index;
  int32_t             j, k, qcount;
  int32_t            *num_quadrants_in_proc;
  int64_t             prev_quadrant, next_quadrant;
  int64_t             weight, weight_sum;
  int64_t             cut, my_lowcut, my_highcut;
  int64_t             send_index, recv_low, recv_high;
  int64_t            *local_weights;    /* cumulative weights by quadrant */
  int64_t            *global_weight_sums;
  int64_t            *send_array;
  p4est_quadrant_t   *q;
  p4est_tree_t       *tree;
  MPI_Request        *send_requests, recv_requests[2];
  MPI_Status          recv_statuses[2];
#endif /* P4EST_MPI */

  P4EST_ASSERT (p4est_is_valid (p4est));
  P4EST_GLOBAL_PRODUCTIONF
    ("Into p4est_partition with %lld total quadrants\n",
     (long long) p4est->global_num_quadrants);

#ifdef P4EST_MPI
  /* this function does nothing in a serial setup */
  if (p4est->mpicomm == MPI_COMM_NULL || p4est->mpisize == 1) {
    P4EST_GLOBAL_PRODUCTION ("Done p4est_partition no shipping\n");
    return;
  }

  /* allocate new quadrant distribution counts */
  num_quadrants_in_proc = P4EST_ALLOC (int32_t, num_procs);

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
    global_weight_sums = P4EST_ALLOC (int64_t, num_procs + 1);
    P4EST_VERBOSEF ("local quadrant count %d\n", p4est->local_num_quadrants);

    /* linearly sum weights across all trees */
    k = 0;
    local_weights[0] = 0;
    for (j = first_tree; j <= last_tree; ++j) {
      tree = sc_array_index (p4est->trees, j);
      for (i = 0; i < tree->quadrants.elem_count; ++i, ++k) {
        q = sc_array_index (&tree->quadrants, i);
        weight = weight_fn (p4est, j, q);
        P4EST_ASSERT (weight >= 0);
        local_weights[k + 1] = local_weights[k] + weight;
      }
    }
    P4EST_ASSERT (k == local_num_quadrants);
    weight_sum = local_weights[local_num_quadrants];
    P4EST_VERBOSEF ("local weight sum %lld\n", (long long) weight_sum);

    /* distribute local weight sums */
    global_weight_sums[0] = 0;
    mpiret = MPI_Allgather (&weight_sum, 1, MPI_LONG_LONG,
                            &global_weight_sums[1], 1, MPI_LONG_LONG,
                            p4est->mpicomm);
    SC_CHECK_MPI (mpiret);

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

    if (rank == 0) {
      for (i = 0; i <= num_procs; ++i) {
        P4EST_GLOBAL_VERBOSEF ("Global weight sum [%d] %lld\n",
                               i, (long long) global_weight_sums[i]);
      }
    }

    /* if all quadrants have zero weight we do nothing */
    if (weight_sum == 0) {
      P4EST_FREE (local_weights);
      P4EST_FREE (global_weight_sums);
      P4EST_FREE (num_quadrants_in_proc);
      P4EST_GLOBAL_PRODUCTION ("Done p4est_partition no shipping\n");
      return;
    }

    /* determine processor ids to send to */
    send_lowest = num_procs;
    send_highest = 0;
    for (i = 1; i <= num_procs; ++i) {
      cut = (weight_sum * i) / num_procs;
      if (global_weight_sums[rank] < cut &&
          cut <= global_weight_sums[rank + 1]) {
        send_lowest = SC_MIN (send_lowest, i);
        send_highest = SC_MAX (send_highest, i);
      }
    }
    /*
     * send low cut to send_lowest..send_highest
     * and high cut to send_lowest-1..send_highest-1
     */
    P4EST_LDEBUGF ("send bounds low %d high %d\n", send_lowest, send_highest);

    num_sends = 2 * (send_highest - send_lowest + 1);
    if (num_sends <= 0) {
      num_sends = 0;
      send_requests = NULL;
      send_array = NULL;
    }
    else {
      send_requests = P4EST_ALLOC (MPI_Request, num_sends);
      send_array = P4EST_ALLOC (int64_t, num_sends);
      k = 0;
      for (i = send_lowest; i <= send_highest; ++i) {
        base_index = 2 * (i - send_lowest);
        if (i < num_procs) {
          /* do binary search in the weight array */
          k = p4est_int64_lower_bound ((weight_sum * i) / num_procs,
                                       local_weights,
                                       local_num_quadrants + 1, k);
          P4EST_ASSERT (k > 0 && k <= local_num_quadrants);
          send_index = send_array[base_index + 1] = k +
            ((rank > 0) ? (p4est->global_last_quad_index[rank - 1] + 1) : 0);

          /* send low bound */
          mpiret = MPI_Isend (&send_array[base_index + 1], 1, MPI_LONG_LONG,
                              i, P4EST_COMM_PARTITION_WEIGHTED_LOW,
                              p4est->mpicomm, &send_requests[base_index + 1]);
          SC_CHECK_MPI (mpiret);
        }
        else {
          k = 0;
          send_index = global_num_quadrants;
          send_requests[base_index + 1] = MPI_REQUEST_NULL;
          send_array[base_index + 1] = -1;
        }
        P4EST_LDEBUGF ("send pos %d index %lld high %d low %d\n",
                       k, (long long) send_index, i - 1, i);

        /* send high bound */
        send_array[base_index] = send_index;
        mpiret = MPI_Isend (&send_array[base_index], 1, MPI_LONG_LONG,
                            i - 1, P4EST_COMM_PARTITION_WEIGHTED_HIGH,
                            p4est->mpicomm, &send_requests[base_index]);
        SC_CHECK_MPI (mpiret);
      }
    }

    /* determine processor ids to receive from and post irecv */
    i = 0;
    my_lowcut = (weight_sum * rank) / num_procs;
    if (my_lowcut == 0) {
      recv_low = 0;
      recv_requests[0] = MPI_REQUEST_NULL;
      low_source = -1;
    }
    else {
      for (; i < num_procs; ++i) {
        if (global_weight_sums[i] < my_lowcut &&
            my_lowcut <= global_weight_sums[i + 1]) {
          P4EST_LDEBUGF ("receive low cut from %d\n", i);
          mpiret = MPI_Irecv (&recv_low, 1, MPI_LONG_LONG, i,
                              P4EST_COMM_PARTITION_WEIGHTED_LOW,
                              p4est->mpicomm, &recv_requests[0]);
          SC_CHECK_MPI (mpiret);
          break;
        }
      }
      P4EST_ASSERT (i < num_procs);
      low_source = i;
    }
    my_highcut = (weight_sum * (rank + 1)) / num_procs;
    if (my_highcut == 0) {
      recv_high = 0;
      recv_requests[1] = MPI_REQUEST_NULL;
      high_source = -1;
    }
    else {
      for (; i < num_procs; ++i) {
        if (global_weight_sums[i] < my_highcut &&
            my_highcut <= global_weight_sums[i + 1]) {
          P4EST_LDEBUGF ("receive high cut from %d\n", i);
          mpiret = MPI_Irecv (&recv_high, 1, MPI_LONG_LONG, i,
                              P4EST_COMM_PARTITION_WEIGHTED_HIGH,
                              p4est->mpicomm, &recv_requests[1]);
          SC_CHECK_MPI (mpiret);
          break;
        }
      }
      P4EST_ASSERT (i < num_procs);
      high_source = i;
    }
    P4EST_LDEBUGF ("my cut low %lld high %lld\n",
                   (long long) my_lowcut, (long long) my_highcut);

    /* free temporary memory */
    P4EST_FREE (local_weights);
    P4EST_FREE (global_weight_sums);

    /* wait for sends and receives to complete */
    if (num_sends > 0) {
      mpiret = MPI_Waitall (num_sends, send_requests, MPI_STATUSES_IGNORE);
      SC_CHECK_MPI (mpiret);
      P4EST_FREE (send_requests);
      P4EST_FREE (send_array);
    }
    mpiret = MPI_Waitall (2, recv_requests, recv_statuses);
    SC_CHECK_MPI (mpiret);
    if (my_lowcut != 0) {
      SC_CHECK_ABORT (recv_statuses[0].MPI_SOURCE == low_source,
                      "Wait low source");
      SC_CHECK_ABORT (recv_statuses[0].MPI_TAG ==
                      P4EST_COMM_PARTITION_WEIGHTED_LOW, "Wait low tag");
      mpiret = MPI_Get_count (&recv_statuses[0], MPI_LONG_LONG, &rcount);
      SC_CHECK_MPI (mpiret);
      SC_CHECK_ABORTF (rcount == 1, "Wait low count %d", rcount);
    }
    if (my_highcut != 0) {
      SC_CHECK_ABORT (recv_statuses[1].MPI_SOURCE == high_source,
                      "Wait high source");
      SC_CHECK_ABORT (recv_statuses[1].MPI_TAG ==
                      P4EST_COMM_PARTITION_WEIGHTED_HIGH, "Wait high tag");
      mpiret = MPI_Get_count (&recv_statuses[1], MPI_LONG_LONG, &rcount);
      SC_CHECK_MPI (mpiret);
      SC_CHECK_ABORTF (rcount == 1, "Wait high count %d", rcount);
    }

    /* communicate the quadrant ranges */
    qcount = recv_high - recv_low;
    P4EST_ASSERT (qcount >= 0);
    P4EST_LDEBUGF ("weighted partition count %d\n", qcount);
    mpiret = MPI_Allgather (&qcount, 1, MPI_INT,
                            num_quadrants_in_proc, 1, MPI_INT,
                            p4est->mpicomm);
    SC_CHECK_MPI (mpiret);

#if(0)
    /* run through the count array and repair zero ranges */
    for (i = 0; i < num_procs; ++i) {
      if (num_quadrants_in_proc[i] == 0) {
        for (p = i - 1; p >= 0; --p) {
          P4EST_ASSERT (num_quadrants_in_proc[p] > 0);
          if (num_quadrants_in_proc[p] > 1) {
            --num_quadrants_in_proc[p];
            ++num_quadrants_in_proc[i];
            break;
          }
        }
        if (p < 0) {
          for (p = i + 1; p < num_procs; ++p) {
            P4EST_ASSERT (num_quadrants_in_proc[p] >= 0);
            if (num_quadrants_in_proc[p] > 1) {
              --num_quadrants_in_proc[p];
              ++num_quadrants_in_proc[i];
              break;
            }
          }
          P4EST_ASSERT (p < num_procs);
        }
      }
    }
#endif
  }

  /* run the partition algorithm with proper quadrant counts */
  global_shipped = p4est_partition_given (p4est, num_quadrants_in_proc);
  P4EST_FREE (num_quadrants_in_proc);

  /* check validity of the p4est */
  P4EST_ASSERT (p4est_is_valid (p4est));
#endif /* P4EST_MPI */

  P4EST_GLOBAL_PRODUCTIONF
    ("Done p4est_partition shipped %lld quadrants %.3g%%\n",
     (long long) global_shipped,
     global_shipped * 100. / global_num_quadrants);
}

unsigned
p4est_checksum (p4est_t * p4est)
{
  uLong               treecrc, crc;
  int32_t             j;
  p4est_tree_t       *tree;
  sc_array_t          checkarray;
#ifdef P4EST_MPI
  int                 mpiret;
  uint32_t            send[2];
  uint32_t           *gather;
#endif

  P4EST_ASSERT (p4est_is_valid (p4est));

  sc_array_init (&checkarray, 4);
  crc = adler32 (0, Z_NULL, 0);
  for (j = p4est->first_local_tree; j <= p4est->last_local_tree; ++j) {
    tree = sc_array_index (p4est->trees, j);
    treecrc = p4est_quadrant_checksum (&tree->quadrants, &checkarray, 0);
    crc = adler32_combine (crc, treecrc, checkarray.elem_count * 4);
  }
  sc_array_reset (&checkarray);

#ifdef P4EST_MPI
  if (p4est->mpirank == 0) {
    gather = P4EST_ALLOC (uint32_t, 2 * p4est->mpisize);
  }
  else {
    gather = NULL;
  }
  send[0] = (uint32_t) crc;
  send[1] = p4est->local_num_quadrants * 12;
  mpiret = MPI_Gather (send, 2, MPI_UNSIGNED, gather, 2, MPI_UNSIGNED,
                       0, p4est->mpicomm);
  SC_CHECK_MPI (mpiret);

  crc = 0;
  if (p4est->mpirank == 0) {
    crc = gather[0];
    for (j = 1; j < p4est->mpisize; ++j) {
      crc = adler32_combine (crc, gather[2 * j + 0], gather[2 * j + 1]);
    }
    P4EST_FREE (gather);
  }
#endif

  return (unsigned) crc;
}

/* EOF p4est.c */
