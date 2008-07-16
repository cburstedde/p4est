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
#include <p8est_algorithms.h>
#include <p8est_bits.h>
#include <p8est_communication.h>
#include <p8est_mesh.h>
#else
#include <p4est_algorithms.h>
#include <p4est_bits.h>
#include <p4est_communication.h>
#include <p4est_mesh.h>
#endif /* !P4_TO_P8 */

typedef struct
{
  int8_t              have_first_count, have_first_load;
  int8_t              have_second_count, have_second_load;
  int                 recv_first_count, recv_second_count;
  int                 send_first_count, send_second_count;
  sc_array_t          send_first, send_second, recv_both;
}
p4est_balance_peer_t;

#ifndef P4_TO_P8

p4est_locidx_t      p4est_initial_quadrants_per_processor = 15;
bool                p4est_refine_recursive = true;
bool                p4est_coarsen_recursive = true;
static int          p4est_uninitialized_key;
void               *P4EST_DATA_UNINITIALIZED = &p4est_uninitialized_key;

#endif /* P4_TO_P8 */

static const size_t number_toread_quadrants = 32;

#ifdef P4EST_MPI
static const int    number_peer_windows = 5;
#endif /* P4EST_MPI */

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
  p4est_quadrant_t   *global_first_position;

  P4EST_GLOBAL_PRODUCTION ("Into " P4EST_STRING "_new\n");
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
  for (level = 0; level < P4EST_QMAXLEVEL; ++level) {
    if (tree_num_quadrants >=
        (num_procs * (p4est_gloidx_t) p4est_initial_quadrants_per_processor)
        / num_trees) {
      break;
    }
    tree_num_quadrants *= P4EST_CHILDREN;
    P4EST_ASSERT (tree_num_quadrants > 0);
  }
  P4EST_ASSERT (level < P4EST_QMAXLEVEL
                && tree_num_quadrants <= (p4est_gloidx_t) P4EST_LOCIDX_MAX);

  /* compute global number of quadrants */
  global_num_quadrants = tree_num_quadrants * num_trees;
  P4EST_GLOBAL_PRODUCTIONF ("New " P4EST_STRING
                            " with %lld trees on %d processors\n",
                            (long long) num_trees, num_procs);
  P4EST_GLOBAL_INFOF ("Initial level %d potential global quadrants"
                      " %lld per tree %lld\n",
                      level, (long long) global_num_quadrants,
                      (long long) tree_num_quadrants);

  /* compute index of first tree for this processor */
  first_quadrant = (global_num_quadrants * rank) / num_procs;
  first_tree = first_quadrant / tree_num_quadrants;
  first_tree_quadrant = first_quadrant - first_tree * tree_num_quadrants;
  last_quadrant = (global_num_quadrants * (rank + 1)) / num_procs - 1;
  P4EST_VERBOSEF
    ("first tree %lld first quadrant %lld global quadrant %lld\n",
     (long long) first_tree, (long long) first_tree_quadrant,
     (long long) first_quadrant);
  P4EST_ASSERT (first_tree_quadrant < tree_num_quadrants);

  /* compute index of last tree for this processor */
  if (first_quadrant <= last_quadrant) {
    last_tree = last_quadrant / tree_num_quadrants;
    last_tree_quadrant = last_quadrant - last_tree * tree_num_quadrants;
    P4EST_VERBOSEF
      ("last tree %lld last quadrant %lld global quadrant %lld\n",
       (long long) last_tree, (long long) last_tree_quadrant,
       (long long) last_quadrant);

    /* check ranges of various integers to be 32bit compatible */
    P4EST_ASSERT (first_tree <= last_tree && last_tree < num_trees);
    P4EST_ASSERT (0 <= first_tree_quadrant && 0 <= last_tree_quadrant);
    P4EST_ASSERT (last_tree_quadrant < tree_num_quadrants);
    if (first_tree == last_tree) {
      P4EST_ASSERT (first_tree_quadrant <= last_tree_quadrant);
    }
  }
  else {
    P4EST_VERBOSE ("Empty processor");
    P4EST_ASSERT (0 <= first_tree && 0 <= first_tree_quadrant);
    first_tree = -1;
    last_tree = -2;
    last_tree_quadrant = -1;
  }

  /* allocate trees and quadrants */
  p4est->trees = sc_array_new (sizeof (p4est_tree_t));
  sc_array_resize (p4est->trees, num_trees);
  for (jl = 0; jl < num_trees; ++jl) {
    tree = sc_array_index (p4est->trees, jl);
    sc_array_init (&tree->quadrants, sizeof (p4est_quadrant_t));
    for (i = 0; i <= P4EST_QMAXLEVEL; ++i) {
      tree->quadrants_per_level[i] = 0;
    }
    for (; i <= P4EST_MAXLEVEL; ++i) {
      tree->quadrants_per_level[i] = -1;
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
#ifdef P4_TO_P8
    P4EST_LDEBUGF ("tree %lld first morton 0x%llx 0x%llx 0x%llx\n",
                   (long long) jl, (long long) a.x,
                   (long long) a.y, (long long) a.z);
#else
    P4EST_LDEBUGF ("tree %lld first morton 0x%llx 0x%llx\n",
                   (long long) jl, (long long) a.x, (long long) a.y);
#endif
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
#ifdef P4_TO_P8
      P4EST_LDEBUGF ("tree %lld last morton 0x%llx 0x%llx 0x%llx\n",
                     (long long) jl, (long long) b.x,
                     (long long) b.y, (long long) b.z);
#else
      P4EST_LDEBUGF ("tree %lld last morton 0x%llx 0x%llx\n",
                     (long long) jl, (long long) b.x, (long long) b.y);
#endif
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
  global_first_position =
    P4EST_ALLOC_ZERO (p4est_quadrant_t, (size_t) num_procs + 1);
  for (i = 0; i <= num_procs; ++i) {
    first_quadrant = (global_num_quadrants * i) / num_procs;
    first_tree = first_quadrant / tree_num_quadrants;
    first_tree_quadrant = first_quadrant - first_tree * tree_num_quadrants;
    p4est_quadrant_set_morton (&c, level, first_tree_quadrant);
    global_first_position[i].x = c.x;
    global_first_position[i].y = c.y;
#ifdef P4_TO_P8
    global_first_position[i].z = c.z;
#endif
    global_first_position[i].level = P4EST_QMAXLEVEL;
    global_first_position[i].p.which_tree = first_tree;
  }
  p4est->global_first_position = global_first_position;

  /* print more statistics */
  P4EST_INFOF ("total local quadrants %lld\n",
               (long long) p4est->local_num_quadrants);

  P4EST_ASSERT (p4est_is_valid (p4est));
  P4EST_GLOBAL_PRODUCTIONF ("Done " P4EST_STRING
                            "_new with %lld total quadrants\n",
                            (long long) p4est->global_num_quadrants);
  return p4est;
}

void
p4est_destroy (p4est_t * p4est)
{
#ifdef P4EST_DEBUG
  size_t              qz;
#endif
  p4est_topidx_t      jt;
  p4est_tree_t       *tree;

  for (jt = 0; jt < p4est->connectivity->num_trees; ++jt) {
    tree = p4est_array_index_topidx (p4est->trees, jt);

#ifdef P4EST_DEBUG
    for (qz = 0; qz < tree->quadrants.elem_count; ++qz) {
      p4est_quadrant_t   *quad = sc_array_index (&tree->quadrants, qz);
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
p4est_copy (p4est_t * input, bool copy_data)
{
  const p4est_topidx_t num_trees = input->connectivity->num_trees;
  const p4est_topidx_t first_tree = input->first_local_tree;
  const p4est_topidx_t last_tree = input->last_local_tree;
  size_t              icount;
  size_t              zz;
  p4est_topidx_t      jt;
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
  for (jt = 0; jt < num_trees; ++jt) {
    itree = p4est_array_index_topidx (input->trees, jt);
    ptree = p4est_array_index_topidx (p4est->trees, jt);
    memcpy (ptree, itree, sizeof (p4est_tree_t));
    sc_array_init (&ptree->quadrants, sizeof (p4est_quadrant_t));
  }
  for (jt = first_tree; jt <= last_tree; ++jt) {
    itree = p4est_array_index_topidx (input->trees, jt);
    iquadrants = &itree->quadrants;
    icount = iquadrants->elem_count;
    ptree = p4est_array_index_topidx (p4est->trees, jt);
    pquadrants = &ptree->quadrants;
    sc_array_resize (pquadrants, icount);
    memcpy (pquadrants->array, iquadrants->array,
            icount * sizeof (p4est_quadrant_t));
    if (p4est->data_size > 0) {
      for (zz = 0; zz < icount; ++zz) {
        iq = sc_array_index (iquadrants, zz);
        pq = sc_array_index (pquadrants, zz);
        pq->p.user_data = sc_mempool_alloc (p4est->user_data_pool);
        memcpy (pq->p.user_data, iq->p.user_data, p4est->data_size);
      }
    }
    else {
      for (zz = 0; zz < icount; ++zz) {
        pq = sc_array_index (pquadrants, zz);
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
  p4est->global_first_position = P4EST_ALLOC (p4est_quadrant_t,
                                              p4est->mpisize + 1);
  memcpy (p4est->global_first_position, input->global_first_position,
          (p4est->mpisize + 1) * sizeof (p4est_quadrant_t));

  /* check for valid p4est and return */
  P4EST_ASSERT (p4est_is_valid (p4est));

  return p4est;
}

void
p4est_refine (p4est_t * p4est, p4est_refine_t refine_fn, p4est_init_t init_fn)
{
  size_t              quadrant_pool_size, data_pool_size;
  bool                dorefine;
  int                 skey, *key = &skey;
  int                 i, maxlevel;
  p4est_topidx_t      nt;
  size_t              incount, current, restpos, movecount;
  sc_list_t          *list;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *q, *qalloc, *qpop;
  p4est_quadrant_t   *c0, *c1, *c2, *c3;
#ifdef P4_TO_P8
  p4est_quadrant_t   *c4, *c5, *c6, *c7;
#endif
  sc_array_t         *tquadrants;

  P4EST_GLOBAL_PRODUCTIONF ("Into " P4EST_STRING
                            "_refine with %lld total quadrants\n",
                            (long long) p4est->global_num_quadrants);
  P4EST_ASSERT (p4est_is_valid (p4est));

  /*
     q points to a quadrant that is an array member
     qalloc is a quadrant that has been allocated through quadrant_pool
     qpop is a quadrant that has been allocated through quadrant_pool
     never mix these two types of quadrant pointers
   */
  list = sc_list_new (NULL);
  p4est->local_num_quadrants = 0;

  /* loop over all local trees */
  for (nt = p4est->first_local_tree; nt <= p4est->last_local_tree; ++nt) {
    tree = p4est_array_index_topidx (p4est->trees, nt);
    tquadrants = &tree->quadrants;
    quadrant_pool_size = p4est->quadrant_pool->elem_count;
    data_pool_size = 0;
    if (p4est->user_data_pool != NULL) {
      data_pool_size = p4est->user_data_pool->elem_count;
    }

    /* initial log message for this tree */
    P4EST_INFOF ("Into refine tree %lld with %llu\n",
                 (long long) nt, (unsigned long long) tquadrants->elem_count);

    /* reset the quadrant counters */
    maxlevel = 0;
    for (i = 0; i <= P4EST_QMAXLEVEL; ++i) {
      tree->quadrants_per_level[i] = 0;
    }

    /* run through the array to find first quadrant to be refined */
    q = NULL;
    dorefine = false;
    incount = tquadrants->elem_count;
    for (current = 0; current < incount; ++current) {
      q = sc_array_index (tquadrants, current);
      dorefine = ((q->level < P4EST_QMAXLEVEL) && refine_fn (p4est, nt, q));
      if (dorefine) {
        break;
      }
      maxlevel = SC_MAX (maxlevel, (int) q->level);
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
          ((p4est_refine_recursive || qpop->p.user_data != key) &&
           qpop->level < P4EST_QMAXLEVEL && refine_fn (p4est, nt, qpop))) {
        dorefine = false;
        sc_array_resize (tquadrants,
                         tquadrants->elem_count + P4EST_CHILDREN - 1);

        /* compute children and prepend them to the list */
        if (qpop->p.user_data != key) {
          p4est_quadrant_free_data (p4est, qpop);
        }
        c0 = qpop;
        c1 = sc_mempool_alloc (p4est->quadrant_pool);
        c2 = sc_mempool_alloc (p4est->quadrant_pool);
        c3 = sc_mempool_alloc (p4est->quadrant_pool);
#ifdef P4_TO_P8
        c4 = sc_mempool_alloc (p4est->quadrant_pool);
        c5 = sc_mempool_alloc (p4est->quadrant_pool);
        c6 = sc_mempool_alloc (p4est->quadrant_pool);
        c7 = sc_mempool_alloc (p4est->quadrant_pool);

        p8est_quadrant_children (qpop, c0, c1, c2, c3, c4, c5, c6, c7);
#else
        p4est_quadrant_children (qpop, c0, c1, c2, c3);
#endif

        c0->p.user_data = key;
        c1->p.user_data = key;
        c2->p.user_data = key;
        c3->p.user_data = key;
#ifdef P4_TO_P8
        c4->p.user_data = key;
        c5->p.user_data = key;
        c6->p.user_data = key;
        c7->p.user_data = key;

        sc_list_prepend (list, c7);
        sc_list_prepend (list, c6);
        sc_list_prepend (list, c5);
        sc_list_prepend (list, c4);
#endif
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
          p4est_quadrant_init_data (p4est, nt, qpop, init_fn);
        }
        q = sc_array_index (tquadrants, current);
        *q = *qpop;
        maxlevel = SC_MAX (maxlevel, (int) qpop->level);
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
    P4EST_INFOF ("Done refine tree %lld now %llu\n",
                 (long long) nt, (unsigned long long) tquadrants->elem_count);
  }

  sc_list_destroy (list);

  /* compute global number of quadrants */
  p4est_comm_count_quadrants (p4est);

  P4EST_ASSERT (p4est_is_valid (p4est));
  P4EST_GLOBAL_PRODUCTIONF ("Done " P4EST_STRING
                            "_refine with %lld total quadrants\n",
                            (long long) p4est->global_num_quadrants);
}

void
p4est_coarsen (p4est_t * p4est, p4est_coarsen_t coarsen_fn,
               p4est_init_t init_fn)
{
  int                 i, maxlevel;
  bool                couldbegood;
  size_t              zz;
  size_t              data_pool_size;
  size_t              incount, removed;
  size_t              cidz, first, last, rest, before;
  p4est_locidx_t      num_quadrants;
  p4est_topidx_t      jt;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *c[P4EST_CHILDREN];
  p4est_quadrant_t   *cfirst, *clast;
  sc_array_t         *tquadrants;

  P4EST_GLOBAL_PRODUCTIONF ("Into " P4EST_STRING
                            "_coarsen with %lld total quadrants\n",
                            (long long) p4est->global_num_quadrants);
  P4EST_ASSERT (p4est_is_valid (p4est));

  /* loop over all local trees */
  for (jt = p4est->first_local_tree; jt <= p4est->last_local_tree; ++jt) {
    tree = p4est_array_index_topidx (p4est->trees, jt);
    tquadrants = &tree->quadrants;
    data_pool_size = 0;
    if (p4est->user_data_pool != NULL) {
      data_pool_size = p4est->user_data_pool->elem_count;
    }
    removed = 0;

    /* initial log message for this tree */
    P4EST_INFOF ("Into coarsen tree %lld with %llu\n",
                 (long long) jt, (unsigned long long) tquadrants->elem_count);

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
    while (rest + P4EST_CHILDREN - 1 - before < incount) {
      couldbegood = true;
      for (zz = 0; zz < P4EST_CHILDREN; ++zz) {
        if (zz < before) {
          c[zz] = sc_array_index (tquadrants, first + zz);
          if (zz != (size_t) p4est_quadrant_child_id (c[zz])) {
            couldbegood = false;
            break;
          }
        }
        else {
          c[zz] = sc_array_index (tquadrants, rest + zz - before);
        }
      }
      if (couldbegood && p4est_quadrant_is_familypv (c) &&
#ifdef P4_TO_P8
          coarsen_fn (p4est, jt, c[0], c[1], c[2], c[3],
                      c[4], c[5], c[6], c[7]) &&
#else
          coarsen_fn (p4est, jt, c[0], c[1], c[2], c[3]) &&
#endif
          true) {
        /* coarsen now */
        for (zz = 0; zz < P4EST_CHILDREN; ++zz) {
          p4est_quadrant_free_data (p4est, c[zz]);
        }
        tree->quadrants_per_level[c[0]->level] -= P4EST_CHILDREN;
        cfirst = c[0];
        p4est_quadrant_parent (c[0], cfirst);
        p4est_quadrant_init_data (p4est, jt, cfirst, init_fn);
        tree->quadrants_per_level[cfirst->level] += 1;
        p4est->local_num_quadrants -= P4EST_CHILDREN - 1;
        removed += P4EST_CHILDREN - 1;

        rest += P4EST_CHILDREN - before;
        if (p4est_coarsen_recursive) {
          last = first;
          cidz = (size_t) p4est_quadrant_child_id (cfirst);
          if (cidz > first)
            first = 0;
          else
            first -= cidz;
        }
        else {
          /* don't coarsen again, move the counters and the hole */
          P4EST_ASSERT (first == last && before == 1);
          if (rest < incount) {
            ++first;
            cfirst = sc_array_index (tquadrants, first);
            clast = sc_array_index (tquadrants, rest);
            *cfirst = *clast;
            last = first;
            ++rest;
          }
        }
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
    for (i = 0; i <= P4EST_QMAXLEVEL; ++i) {
      P4EST_ASSERT (tree->quadrants_per_level[i] >= 0);
      num_quadrants += tree->quadrants_per_level[i];    /* same type */
      if (tree->quadrants_per_level[i] > 0) {
        maxlevel = i;
      }
    }
    tree->maxlevel = (int8_t) maxlevel;

    /* do some sanity checks */
    P4EST_ASSERT (num_quadrants == (p4est_locidx_t) tquadrants->elem_count);
    P4EST_ASSERT (tquadrants->elem_count == incount - removed);
    if (p4est->user_data_pool != NULL) {
      P4EST_ASSERT (data_pool_size - removed ==
                    p4est->user_data_pool->elem_count);
    }
    P4EST_ASSERT (p4est_tree_is_sorted (tree));
    P4EST_ASSERT (p4est_tree_is_complete (tree));

    /* final log message for this tree */
    P4EST_INFOF ("Done coarsen tree %lld now %llu\n",
                 (long long) jt, (unsigned long long) tquadrants->elem_count);
  }

  /* compute global number of quadrants */
  p4est_comm_count_quadrants (p4est);

  P4EST_ASSERT (p4est_is_valid (p4est));
  P4EST_GLOBAL_PRODUCTIONF ("Done " P4EST_STRING
                            "_coarsen with %lld total quadrants\n",
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
                        p4est_topidx_t qtree, bool inter_tree,
                        const p4est_quadrant_t * q,
                        const p4est_quadrant_t * insul,
                        int *first_peer, int *last_peer)
{
  const int           rank = p4est->mpirank;
  bool                found;
  int                 back, pos;
  int                 owner, first_owner, last_owner;
  p4est_quadrant_t    ld, *s;
  p4est_balance_peer_t *peer;

  P4EST_QUADRANT_INIT (&ld);

  /* querying insul is equivalent to querying first descendent */
  first_owner = p4est_comm_find_owner (p4est, qtree, insul, rank);
  /* querying last descendent */
  p4est_quadrant_last_descendent (insul, &ld, P4EST_QMAXLEVEL);
  last_owner = p4est_comm_find_owner (p4est, qtree, &ld, rank);

  /* send to all processors possibly intersecting insulation */
  for (owner = first_owner; owner <= last_owner; ++owner) {
    if (owner == rank && !inter_tree) {
      continue;
    }
    peer = sc_array_index_int (peers, owner);
    /* avoid duplicates in the send array */
    found = false;
    for (back = 0; back < P4EST_INSUL - 1; ++back) {
      pos = (int) peer->send_first.elem_count - back - 1;
      if (pos < 0) {
        break;
      }
      s = sc_array_index_int (&peer->send_first, pos);
      if (p4est_quadrant_is_equal (s, q) && s->p.piggy2.which_tree == qtree) {
        found = true;
        break;
      }
    }
    if (found) {
      continue;
    }

    /* copy quadrant into shipping list */
    s = sc_array_push (&peer->send_first);
    *s = *q;
    s->p.piggy2.which_tree = qtree;     /* piggy back tree id */

    /* update lowest and highest peer */
    *first_peer = SC_MIN (owner, *first_peer);
    *last_peer = SC_MAX (owner, *last_peer);
  }
}

static void
p4est_balance_response (p4est_t * p4est, int peer_id,
                        p4est_balance_peer_t * peer)
{
#ifdef P4EST_DEBUG
  const p4est_topidx_t first_tree = p4est->first_local_tree;
  const p4est_topidx_t last_tree = p4est->last_local_tree;
#endif /* P4EST_DEBUG */
  size_t              zz, qcount;
  p4est_topidx_t      prev, qtree;
  p4est_topidx_t      num_receive_trees, nt;
  p4est_topidx_t     *pi;
  sc_array_t         *qarray, tree_array;
  p4est_quadrant_t   *q;

  qarray = &peer->recv_both;
  qcount = qarray->elem_count;
  P4EST_ASSERT ((size_t) peer->recv_first_count == qcount);

  /* build list of received tree ids */
  prev = -1;
  num_receive_trees = 0;
  sc_array_init (&tree_array, sizeof (p4est_topidx_t));
  for (zz = 0; zz < qcount; ++zz) {
    q = sc_array_index (qarray, zz);
    qtree = q->p.piggy2.which_tree;
    P4EST_ASSERT (first_tree <= qtree && qtree <= last_tree);
    P4EST_ASSERT (qtree >= prev);
    if (qtree > prev) {
      sc_array_resize (&tree_array, (size_t) num_receive_trees + 1);
      pi = p4est_array_index_topidx (&tree_array, num_receive_trees);
      *pi = qtree;
      ++num_receive_trees;
      prev = qtree;
    }
  }
  P4EST_LDEBUGF ("first load from %d into %lld trees\n", peer_id,
                 (long long) num_receive_trees);
  P4EST_ASSERT (num_receive_trees == (p4est_topidx_t) tree_array.elem_count);

  /* loop over the trees to receive into and update overlap quadrants */
  for (nt = 0; nt < num_receive_trees; ++nt) {
    pi = p4est_array_index_topidx (&tree_array, nt);
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
  int                 i, j, k, l, m, which;
  int                 face;
  int                 first_peer, last_peer;
  int                 over_peer_count;
  bool                quad_contact[2 * P4EST_DIM];
  bool                any_face, tree_contact[2 * P4EST_DIM];
  bool                tree_fully_owned, full_tree[2];
  int8_t             *tree_flags;
  size_t              zz, treecount, ctree;
  size_t              qcount, qbytes, offset, obytes;
  size_t              data_pool_size;
  size_t              all_incount, all_outcount;
  p4est_qcoord_t      qh;
  const p4est_qcoord_t rh = P4EST_ROOT_LEN;
  p4est_topidx_t      qtree, nt;
  p4est_topidx_t      first_tree, last_tree;
  p4est_locidx_t      skipped;
  p4est_balance_peer_t *peer;
  p4est_tree_t       *tree;
  p4est_quadrant_t    mylow, nextlow;
  p4est_quadrant_t    tosend, insulq, tempq;
  p4est_quadrant_t   *q, *s;
  p4est_connectivity_t *conn = p4est->connectivity;
  sc_array_t         *peers, *qarray, *tquadrants;
#ifndef P4_TO_P8
  int                 transform, corner;
  p4est_corner_transform_t *ct;
  sc_array_t          ctransforms, *cta;
#else
  int                 ftransform[9], edge, corner;
  bool                face_axis[3], contact_face_only, contact_edge_only;
  size_t              etree;
  p8est_edge_info_t   ei;
  p8est_edge_transform_t *et;
  p8est_corner_info_t ci;
  p8est_corner_transform_t *ct;
  sc_array_t         *eta, *cta;
#endif
#ifdef P4EST_MPI
#ifdef P4EST_DEBUG
  unsigned            checksum;
  sc_array_t          checkarray;
  int                 ltotal[2], gtotal[2];
#endif /* P4EST_DEBUG */
  const int           twopeerw = 2 * number_peer_windows;
  int                 mpiret, rcount;
  int                 first_bound, last_bound;
  int                 request_first_count, request_second_count, outcount;
  int                 request_send_count, total_send_count, total_recv_count;
  int                 lastw, nwin;
  int                 send_zero[2], send_load[2];
  int                 recv_zero[2], recv_load[2];
  int                *wait_indices;
  int                 prev, start, end;
  int                 rip, eff_peer_count;
  int                 length, shortest_window, shortest_length;
  int                 peer_windows[twopeerw];
  int                *peer_boundaries;
  MPI_Request        *requests_first, *requests_second;
  MPI_Request        *send_requests_first_count, *send_requests_first_load;
  MPI_Request        *send_requests_second_count, *send_requests_second_load;
  MPI_Status         *recv_statuses, *jstatus;
#endif /* P4EST_MPI */

  P4EST_GLOBAL_PRODUCTIONF ("Into " P4EST_STRING
                            "_balance with %lld total quadrants\n",
                            (long long) p4est->global_num_quadrants);
  P4EST_ASSERT (p4est_is_valid (p4est));

  /* prepare sanity checks */
  data_pool_size = 0;
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
  for (nt = 0; nt < conn->num_trees; ++nt) {
    tree_flags[nt] = 0x00;
  }

#ifdef P4EST_MPI
  /* will contain first and last peer (inclusive) for each processor */
  peer_boundaries = P4EST_ALLOC (int, twopeerw * num_procs);
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
  sc_array_resize (peers, (size_t) num_procs);
  for (i = 0; i < num_procs; ++i) {
    peer = sc_array_index_int (peers, i);
    sc_array_init (&peer->send_first, sizeof (p4est_quadrant_t));
    sc_array_init (&peer->send_second, sizeof (p4est_quadrant_t));
    sc_array_init (&peer->recv_both, sizeof (p4est_quadrant_t));
    peer->send_first_count = peer->send_second_count = 0;
    peer->recv_first_count = peer->recv_second_count = 0;
    peer->have_first_count = peer->have_first_load = 0;
    peer->have_second_count = peer->have_second_load = 0;
  }
#ifndef P4_TO_P8
  cta = &ctransforms;
  sc_array_init (cta, sizeof (p4est_corner_transform_t));
#else
  eta = &ei.edge_transforms;
  sc_array_init (eta, sizeof (p8est_edge_transform_t));
  cta = &ci.corner_transforms;
  sc_array_init (cta, sizeof (p8est_corner_transform_t));
#endif /* !P4_TO_P8 */

  /* compute first quadrant on finest level */
  mylow.x = p4est->global_first_position[rank].x;
  mylow.y = p4est->global_first_position[rank].y;
#ifdef P4_TO_P8
  mylow.z = p4est->global_first_position[rank].z;
#endif
  mylow.level = P4EST_QMAXLEVEL;

  /* and the first finest quadrant of the next processor */
  nextlow.x = p4est->global_first_position[rank + 1].x;
  nextlow.y = p4est->global_first_position[rank + 1].y;
#ifdef P4_TO_P8
  nextlow.z = p4est->global_first_position[rank + 1].z;
#endif
  nextlow.level = P4EST_QMAXLEVEL;

  /* loop over all local trees to assemble first send list */
  first_tree = p4est->first_local_tree;
  last_tree = p4est->last_local_tree;
  first_peer = num_procs;
  last_peer = -1;
  all_incount = 0;
  skipped = 0;
  for (nt = first_tree; nt <= last_tree; ++nt) {
    p4est_comm_tree_info (p4est, nt, full_tree, tree_contact, NULL, NULL);
    tree_fully_owned = full_tree[0] && full_tree[1];
    any_face = false;
    for (face = 0; face < 2 * P4EST_DIM; ++face) {
      any_face = any_face || tree_contact[face];
    }
    if (any_face) {
      tree_flags[nt] |= any_face_flag;
    }
    tree = p4est_array_index_topidx (p4est->trees, nt);
    tquadrants = &tree->quadrants;
    all_incount += tquadrants->elem_count;

    /* initial log message for this tree */
    P4EST_INFOF ("Into balance tree %lld with %llu\n",
                 (long long) nt, (unsigned long long) tquadrants->elem_count);

    /* local balance first pass */
    p4est_balance_subtree (p4est, nt, init_fn);
    treecount = tquadrants->elem_count;
    P4EST_VERBOSEF ("Balance tree %lld A %llu\n",
                    (long long) nt, (unsigned long long) treecount);

    /* check if this tree is not shared with other processors */
    if (tree_fully_owned) {
      /* all quadrants in this tree are owned by me */
      tree_flags[nt] |= fully_owned_flag;
      if (!any_face) {
        /* this tree is isolated, no balance between trees */
        continue;
      }
    }

    /* identify boundary quadrants and prepare them to be sent */
    for (zz = 0; zz < treecount; ++zz) {
      /* this quadrant may be on the boundary with a range of processors */
      q = sc_array_index (tquadrants, zz);
      qh = P4EST_QUADRANT_LEN (q->level);
      if (p4est_comm_neighborhood_owned (p4est, nt,
                                         full_tree, tree_contact, q)) {
        /* this quadrant's 3x3 neighborhood is onwed by this processor */
        ++skipped;
        continue;
      }

#ifdef P4_TO_P8
      for (m = 0; m < 3; ++m) {
#if 0
      }
#endif
#else
      m = 0;
#endif
      for (k = 0; k < 3; ++k) {
        for (l = 0; l < 3; ++l) {
          which = m * 9 + k * 3 + l;    /* 2D: 0..8, 3D: 0..26 */
          /* exclude myself from the queries */
          if (which == P4EST_INSUL / 2) {
            continue;
          }
          /* may modify insulq below, never modify q itself! */
          insulq = *q;
          insulq.x += (l - 1) * qh;
          insulq.y += (k - 1) * qh;
#ifdef P4_TO_P8
          insulq.z += (m - 1) * qh;
#endif

          /* check boundary status of insulation quadrant */
#ifndef P4_TO_P8
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
              p4est_find_corner_transform (conn, nt, corner, cta);
              for (ctree = 0; ctree < cta->elem_count; ++ctree) {
                ct = sc_array_index (cta, ctree);
                tosend = *q;
                p4est_quadrant_transform_corner (&tosend, (int) ct->ncorner,
                                                 false);
                p4est_quadrant_transform_corner (&insulq, (int) ct->ncorner,
                                                 true);
                p4est_balance_schedule (p4est, peers, ct->ntree, true,
                                        &tosend, &insulq, &first_peer,
                                        &last_peer);
                /* insulq is now invalid don't use it below */
              }
            }
            else {
              /* this quadrant goes across a face */
              for (face = 0; face < 4; ++face) {
                if (quad_contact[face] && tree_contact[face]) {
                  qtree = conn->tree_to_tree[4 * nt + face];
                  break;
                }
              }
              if (face == 4) {
                /* this quadrant ran across a face with no neighbor */
                continue;
              }
              /* transform both q and insulq into the neighbor's coordinates */
              transform = p4est_find_face_transform (conn, nt, face);
              tempq = *q;
              p4est_quadrant_translate_face (&tempq, face);
              p4est_quadrant_transform_face (&tempq, &tosend, transform);
              p4est_quadrant_translate_face (&insulq, face);
              p4est_quadrant_transform_face (&insulq, &tempq, transform);
              p4est_balance_schedule (p4est, peers, qtree, true,
                                      &tosend, &tempq,
                                      &first_peer, &last_peer);
            }
          }
#else /* P4_TO_P8 */
          quad_contact[0] = (insulq.x < 0);
          quad_contact[1] = (insulq.x >= rh);
          quad_contact[2] = (insulq.y < 0);
          quad_contact[3] = (insulq.y >= rh);
          quad_contact[4] = (insulq.z < 0);
          quad_contact[5] = (insulq.z >= rh);
          face_axis[0] = quad_contact[0] || quad_contact[1];
          face_axis[1] = quad_contact[2] || quad_contact[3];
          face_axis[2] = quad_contact[4] || quad_contact[5];
          contact_face_only = contact_edge_only = false;
          face = edge = -1;
          if (face_axis[0] || face_axis[1] || face_axis[2]) {
            /* this quadrant is relevant for inter-tree balancing */
            if (!face_axis[1] && !face_axis[2]) {
              contact_face_only = true;
              face = 0 + quad_contact[1];
            }
            else if (!face_axis[0] && !face_axis[2]) {
              contact_face_only = true;
              face = 2 + quad_contact[3];
            }
            else if (!face_axis[0] && !face_axis[1]) {
              contact_face_only = true;
              face = 4 + quad_contact[5];
            }
            else if (!face_axis[0]) {
              contact_edge_only = true;
              edge = 0 + 2 * quad_contact[5] + quad_contact[3];
            }
            else if (!face_axis[1]) {
              contact_edge_only = true;
              edge = 4 + 2 * quad_contact[5] + quad_contact[1];
            }
            else if (!face_axis[2]) {
              contact_edge_only = true;
              edge = 8 + 2 * quad_contact[3] + quad_contact[1];
            }
            if (contact_face_only) {
              /* square contact across a face */
              P4EST_ASSERT (!contact_edge_only && face >= 0 && face < 6);
              P4EST_ASSERT (quad_contact[face]);
              qtree = p8est_find_face_transform (conn, nt, face, ftransform);
              if (qtree >= 0) {
                P4EST_ASSERT (tree_contact[face]);
                p8est_quadrant_transform_face (q, &tosend, ftransform);
                p8est_quadrant_transform_face (&insulq, &tempq, ftransform);
                p4est_balance_schedule (p4est, peers, qtree, true,
                                        &tosend, &tempq,
                                        &first_peer, &last_peer);
              }
              else {
                /* goes across a face with no neighbor */
                P4EST_ASSERT (!tree_contact[face]);
              }
            }
            else if (contact_edge_only) {
              /* this quadrant crosses an edge */
              P4EST_ASSERT (!contact_face_only && edge >= 0 && edge < 12);
              p8est_find_edge_transform (conn, nt, edge, &ei);
              for (etree = 0; etree < eta->elem_count; ++etree) {
                et = sc_array_index (eta, etree);
                p8est_quadrant_transform_edge (q, &tosend, &ei, et, false);
                p8est_quadrant_transform_edge (&insulq, &tempq, &ei, et,
                                               true);
                p4est_balance_schedule (p4est, peers, et->ntree, true,
                                        &tosend, &tempq,
                                        &first_peer, &last_peer);
              }
            }
            else {
              /* this quadrant crosses a corner */
              P4EST_ASSERT (face_axis[0] && face_axis[1] && face_axis[2]);
              corner =
                4 * quad_contact[5] + 2 * quad_contact[3] + quad_contact[1];
              P4EST_ASSERT (p8est_quadrant_touches_corner (q, corner));
              P4EST_ASSERT (p8est_quadrant_touches_corner (&insulq, corner));
              p8est_find_corner_transform (conn, nt, corner, &ci);
              for (ctree = 0; ctree < cta->elem_count; ++ctree) {
                ct = sc_array_index (cta, ctree);
                tosend = *q;
                p8est_quadrant_transform_corner (&tosend, (int) ct->ncorner,
                                                 false);
                tempq = insulq;
                p8est_quadrant_transform_corner (&tempq, (int) ct->ncorner,
                                                 true);
                p4est_balance_schedule (p4est, peers, ct->ntree, true,
                                        &tosend, &tempq, &first_peer,
                                        &last_peer);
              }
            }
          }
#endif /* !P4_TO_P8 */
          else {
            /* no inter-tree contact */
            p4est_balance_schedule (p4est, peers, nt, false,
                                    q, &insulq, &first_peer, &last_peer);
          }
        }
      }
#ifdef P4_TO_P8
#if 0
      {
#endif
      }
#endif
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
    peer = sc_array_index_int (peers, j);
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
    peer = sc_array_index_int (peers, j);
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
      P4EST_LDEBUGF ("Balance A send %llu quadrants to %d\n",
                     (unsigned long long) qcount, j);
      ++send_load[0];
    }
    else {
      ++send_zero[0];
    }
    peer->send_first_count = (int) qcount;
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
      mpiret = MPI_Isend (peer->send_first.array, (int) qbytes, MPI_BYTE,
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
    peer = sc_array_index_int (peers, j);
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
    peer = sc_array_index_int (peers, j);
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
      peer = sc_array_index_int (peers, j);
      P4EST_ASSERT (!peer->have_first_load);
      if (!peer->have_first_count) {
        /* verify message size */
        P4EST_ASSERT (jstatus->MPI_TAG == P4EST_COMM_BALANCE_FIRST_COUNT);
        mpiret = MPI_Get_count (jstatus, MPI_INT, &rcount);
        SC_CHECK_MPI (mpiret);
        SC_CHECK_ABORTF (rcount == 1, "Receive count mismatch A %d", rcount);

        /* process the count information received */
        peer->have_first_count = 1;
        qcount = (size_t) peer->recv_first_count;
        if (qcount > 0) {
          /* received nonzero count, post receive for load */
          P4EST_LDEBUGF ("Balance A recv %llu quadrants from %d\n",
                         (unsigned long long) qcount, j);
          P4EST_ASSERT (peer->recv_both.elem_count == 0);
          sc_array_resize (&peer->recv_both, qcount);
          total_recv_count += qcount;
          qbytes = qcount * sizeof (p4est_quadrant_t);
          P4EST_ASSERT (requests_first[j] == MPI_REQUEST_NULL);
          mpiret = MPI_Irecv (peer->recv_both.array, (int) qbytes, MPI_BYTE,
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
                         (int) sizeof (p4est_quadrant_t),
                         "Receive load mismatch A %d %dx%llu", rcount,
                         peer->recv_first_count,
                         (unsigned long long) sizeof (p4est_quadrant_t));

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
          P4EST_LDEBUGF ("Balance B send %llu quadrants to %d\n",
                         (unsigned long long) qcount, j);
          ++send_load[1];
        }
        else {
          ++send_zero[1];
        }
        peer->send_second_count = (int) qcount;
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
          mpiret = MPI_Isend (peer->send_second.array, (int) qbytes, MPI_BYTE,
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
  peer = sc_array_index_int (peers, rank);
  sc_array_sort (&peer->send_first, p4est_quadrant_compare_piggy);
  offset = peer->send_first.elem_count;
  peer->recv_first_count = peer->send_first_count = (int) offset;
  obytes = offset * sizeof (p4est_quadrant_t);
  qarray = &peer->recv_both;
  sc_array_resize (qarray, offset);
  memcpy (qarray->array, peer->send_first.array, obytes);
  p4est_balance_response (p4est, rank, peer);
  qcount = peer->send_second.elem_count;
  peer->recv_second_count = peer->send_second_count = (int) qcount;
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
      peer = sc_array_index_int (peers, j);
      P4EST_ASSERT (!peer->have_second_load);
      if (!peer->have_second_count) {
        /* verify message size */
        P4EST_ASSERT (jstatus->MPI_TAG == P4EST_COMM_BALANCE_SECOND_COUNT);
        mpiret = MPI_Get_count (jstatus, MPI_INT, &rcount);
        SC_CHECK_MPI (mpiret);
        SC_CHECK_ABORTF (rcount == 1, "Receive count mismatch B %d", rcount);

        /* process the count information received */
        peer->have_second_count = 1;
        qcount = (size_t) peer->recv_second_count;
        if (qcount > 0) {
          /* received nonzero count, post receive for load */
          P4EST_LDEBUGF ("Balance B recv %llu quadrants from %d\n",
                         (unsigned long long) qcount, j);
          offset = peer->recv_both.elem_count;
          P4EST_ASSERT (offset == (size_t) peer->recv_first_count);
          obytes = offset * sizeof (p4est_quadrant_t);
          sc_array_resize (&peer->recv_both, offset + qcount);
          total_recv_count += qcount;
          qbytes = qcount * sizeof (p4est_quadrant_t);
          P4EST_ASSERT (requests_second[j] == MPI_REQUEST_NULL);
          mpiret = MPI_Irecv (peer->recv_both.array + obytes, (int) qbytes,
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
                         (int) sizeof (p4est_quadrant_t),
                         "Receive load mismatch B %d %dx%llu", rcount,
                         peer->recv_second_count,
                         (unsigned long long) sizeof (p4est_quadrant_t));

        /* received load, close this request */
        peer->have_second_load = 1;
        P4EST_ASSERT (requests_second[j] == MPI_REQUEST_NULL);
        --request_second_count;

#ifdef P4EST_DEBUG
        checksum = p4est_quadrant_checksum (&peer->recv_both, &checkarray,
                                            (size_t) peer->recv_first_count);
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
    peer = sc_array_index_int (peers, j);
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
    peer = sc_array_index_int (peers, j);
    qarray = &peer->recv_both;
    qcount = qarray->elem_count;
    P4EST_ASSERT (qcount ==
                  (size_t) (peer->recv_first_count +
                            peer->recv_second_count));
    P4EST_ASSERT (peer->send_first_count ==
                  (int) peer->send_first.elem_count);
    P4EST_ASSERT (peer->send_second_count ==
                  (int) peer->send_second.elem_count);
    if (qcount == 0) {
      continue;
    }

    /* merge received quadrants into correct tree */
    for (zz = 0; zz < qcount; ++zz) {
      s = sc_array_index (qarray, zz);
      P4EST_ASSERT (p4est_quadrant_is_extended (s));
      qtree = s->p.piggy2.which_tree;
      if (qtree < first_tree || qtree > last_tree) {
        /* this is a corner/edge quadrant from the second pass of balance */
        P4EST_ASSERT (zz >= (size_t) peer->recv_first_count);
        P4EST_ASSERT (0 <= qtree && qtree < conn->num_trees);
#ifndef P4_TO_P8
        P4EST_ASSERT ((s->x < 0 && s->y < 0) || (s->x < 0 && s->y >= rh) ||
                      (s->x >= rh && s->y < 0) || (s->x >= rh && s->y >= rh));
#else
        face_axis[0] = (s->x < 0 || s->x >= rh);
        face_axis[1] = (s->y < 0 || s->y >= rh);
        face_axis[2] = (s->z < 0 || s->z >= rh);
        P4EST_ASSERT ((face_axis[0] && face_axis[1]) ||
                      (face_axis[0] && face_axis[2]) ||
                      (face_axis[1] && face_axis[2]));
#endif
        continue;
      }
      tree = sc_array_index (p4est->trees, qtree);
      q = sc_array_push (&tree->quadrants);
      *q = *s;
      ++tree->quadrants_per_level[q->level];
      tree->maxlevel = (int8_t) SC_MAX (tree->maxlevel, q->level);
      ++p4est->local_num_quadrants;
      p4est_quadrant_init_data (p4est, qtree, q, init_fn);
    }
  }

  /* rebalance and clamp result back to original tree boundaries */
  p4est->local_num_quadrants = 0;
  for (nt = first_tree; nt <= last_tree; ++nt) {
    /* check if we are the only processor in an isolated tree */
    tree = p4est_array_index_topidx (p4est->trees, nt);
    tquadrants = &tree->quadrants;
    treecount = tquadrants->elem_count;
    if (!(tree_flags[nt] & fully_owned_flag) ||
        (tree_flags[nt] & any_face_flag)) {
      /* we have most probably received quadrants, run sort and balance */
      sc_array_sort (tquadrants, p4est_quadrant_compare);
      p4est_balance_subtree (p4est, nt, init_fn);
      P4EST_VERBOSEF ("Balance tree %lld B %llu to %llu\n",
                      (long long) nt,
                      (unsigned long long) treecount,
                      (unsigned long long) tquadrants->elem_count);
    }
    p4est->local_num_quadrants += tquadrants->elem_count;
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
  for (nt = first_tree; nt <= last_tree; ++nt) {
    tree = p4est_array_index_topidx (p4est->trees, nt);
    all_outcount += tree->quadrants.elem_count;

    /* final log message for this tree */
    P4EST_INFOF ("Done balance tree %lld now %llu\n",
                 (long long) nt,
                 (unsigned long long) tree->quadrants.elem_count);
  }

  /* cleanup temporary storage */
  P4EST_FREE (tree_flags);
  for (i = 0; i < num_procs; ++i) {
    peer = sc_array_index_int (peers, i);
    sc_array_reset (&peer->send_first);
    sc_array_reset (&peer->send_second);
    sc_array_reset (&peer->recv_both);
  }
  sc_array_destroy (peers);

#ifdef P4_TO_P8
  sc_array_reset (eta);
#endif
  sc_array_reset (cta);

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
  P4EST_ASSERT ((p4est_locidx_t) all_outcount == p4est->local_num_quadrants);
  P4EST_ASSERT (all_outcount >= all_incount);
  if (p4est->user_data_pool != NULL) {
    P4EST_ASSERT (data_pool_size + all_outcount - all_incount ==
                  p4est->user_data_pool->elem_count);
  }
  P4EST_ASSERT (p4est_is_valid (p4est));
  P4EST_ASSERT (p4est_is_balanced (p4est));
  P4EST_VERBOSEF ("Balance skipped %lld\n", (long long) skipped);
  P4EST_GLOBAL_PRODUCTIONF ("Done " P4EST_STRING
                            "_balance with %lld total quadrants\n",
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
  size_t              lz;
  ssize_t             lowers;
  p4est_topidx_t      nt;
  p4est_locidx_t      kl, qlocal;
  p4est_locidx_t     *num_quadrants_in_proc;
  p4est_gloidx_t      prev_quadrant, next_quadrant;
  p4est_gloidx_t      send_index, recv_low, recv_high, qcount;
  p4est_gloidx_t     *send_array;
  int64_t             weight, weight_sum;
  int64_t             cut, my_lowcut, my_highcut;
  int64_t            *local_weights;    /* cumulative weights by quadrant */
  int64_t            *global_weight_sums;
  p4est_quadrant_t   *q;
  p4est_tree_t       *tree;
  MPI_Request        *send_requests, recv_requests[2];
  MPI_Status          recv_statuses[2];
#endif /* P4EST_MPI */

  P4EST_ASSERT (p4est_is_valid (p4est));
  P4EST_GLOBAL_PRODUCTIONF
    ("Into " P4EST_STRING
     "_partition with %lld total quadrants\n",
     (long long) p4est->global_num_quadrants);

#ifdef P4EST_MPI
  /* this function does nothing in a serial setup */
  if (p4est->mpicomm == MPI_COMM_NULL || p4est->mpisize == 1) {
    P4EST_GLOBAL_PRODUCTION ("Done " P4EST_STRING "_partition no shipping\n");
    return;
  }

  /* allocate new quadrant distribution counts */
  num_quadrants_in_proc = P4EST_ALLOC (p4est_locidx_t, num_procs);

  if (weight_fn == NULL) {
    /* Divide up the quadants equally */
    for (p = 0, next_quadrant = 0; p < num_procs; ++p) {
      prev_quadrant = next_quadrant;
      next_quadrant = (global_num_quadrants * (p + 1)) / num_procs;
      qcount = next_quadrant - prev_quadrant;
      P4EST_ASSERT (0 <= qcount
                    && qcount <= (p4est_gloidx_t) P4EST_LOCIDX_MAX);
      num_quadrants_in_proc[p] = (p4est_locidx_t) (qcount);
    }
  }
  else {
    /* do a weighted partition */
    local_weights = P4EST_ALLOC (int64_t, local_num_quadrants + 1);
    global_weight_sums = P4EST_ALLOC (int64_t, num_procs + 1);
    P4EST_VERBOSEF ("local quadrant count %lld\n",
                    (long long) local_num_quadrants);

    /* linearly sum weights across all trees */
    kl = 0;
    local_weights[0] = 0;
    for (nt = first_tree; nt <= last_tree; ++nt) {
      tree = p4est_array_index_topidx (p4est->trees, nt);
      for (lz = 0; lz < tree->quadrants.elem_count; ++lz, ++kl) {
        q = sc_array_index (&tree->quadrants, lz);
        weight = (int64_t) weight_fn (p4est, nt, q);
        P4EST_ASSERT (weight >= 0);
        local_weights[kl + 1] = local_weights[kl] + weight;
      }
    }
    P4EST_ASSERT (kl == local_num_quadrants);
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
      for (kl = 0; kl <= local_num_quadrants; ++kl) {
        local_weights[kl] += weight_sum;
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
      P4EST_GLOBAL_PRODUCTION ("Done " P4EST_STRING
                               "_partition no shipping\n");
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
      send_array = P4EST_ALLOC (p4est_gloidx_t, num_sends);
      lowers = 0;
      for (i = send_lowest; i <= send_highest; ++i) {
        base_index = 2 * (i - send_lowest);
        if (i < num_procs) {
          /* do binary search in the weight array */
          lowers = p4est_int64_lower_bound ((weight_sum * i) / num_procs,
                                            local_weights,
                                            (size_t) local_num_quadrants + 1,
                                            (size_t) lowers);
          P4EST_ASSERT (lowers > 0
                        && (p4est_locidx_t) lowers <= local_num_quadrants);
          send_index = send_array[base_index + 1] =
            (p4est_gloidx_t) lowers +
            ((rank > 0) ? (p4est->global_last_quad_index[rank - 1] + 1) : 0);

          /* send low bound */
          mpiret =
            MPI_Isend (&send_array[base_index + 1], 1, P4EST_MPI_GLOIDX, i,
                       P4EST_COMM_PARTITION_WEIGHTED_LOW, p4est->mpicomm,
                       &send_requests[base_index + 1]);
          SC_CHECK_MPI (mpiret);
        }
        else {
          lowers = 0;
          send_index = global_num_quadrants;
          send_requests[base_index + 1] = MPI_REQUEST_NULL;
          send_array[base_index + 1] = -1;
        }
        P4EST_LDEBUGF ("send pos %lld index %lld high %d low %d\n",
                       (long long) lowers, (long long) send_index, i - 1, i);

        /* send high bound */
        send_array[base_index] = send_index;
        mpiret = MPI_Isend (&send_array[base_index], 1, P4EST_MPI_GLOIDX,
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
          mpiret = MPI_Irecv (&recv_low, 1, P4EST_MPI_GLOIDX, i,
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
          mpiret = MPI_Irecv (&recv_high, 1, P4EST_MPI_GLOIDX, i,
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
      mpiret = MPI_Get_count (&recv_statuses[0], P4EST_MPI_GLOIDX, &rcount);
      SC_CHECK_MPI (mpiret);
      SC_CHECK_ABORTF (rcount == 1, "Wait low count %d", rcount);
    }
    if (my_highcut != 0) {
      SC_CHECK_ABORT (recv_statuses[1].MPI_SOURCE == high_source,
                      "Wait high source");
      SC_CHECK_ABORT (recv_statuses[1].MPI_TAG ==
                      P4EST_COMM_PARTITION_WEIGHTED_HIGH, "Wait high tag");
      mpiret = MPI_Get_count (&recv_statuses[1], P4EST_MPI_GLOIDX, &rcount);
      SC_CHECK_MPI (mpiret);
      SC_CHECK_ABORTF (rcount == 1, "Wait high count %d", rcount);
    }

    /* communicate the quadrant ranges */
    qcount = recv_high - recv_low;
    P4EST_LDEBUGF ("weighted partition count %lld\n", (long long) qcount);
    P4EST_ASSERT (qcount >= 0 && qcount <= (p4est_gloidx_t) P4EST_LOCIDX_MAX);
    qlocal = (p4est_locidx_t) qcount;
    mpiret = MPI_Allgather (&qlocal, 1, P4EST_MPI_LOCIDX,
                            num_quadrants_in_proc, 1, P4EST_MPI_LOCIDX,
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
    ("Done " P4EST_STRING "_partition shipped %lld quadrants %.3g%%\n",
     (long long) global_shipped,
     global_shipped * 100. / global_num_quadrants);
}

unsigned
p4est_checksum (p4est_t * p4est)
{
  uLong               treecrc, crc;
  size_t              scount, ssum;
  p4est_topidx_t      nt;
  p4est_tree_t       *tree;
  sc_array_t          checkarray;
#ifdef P4EST_MPI
  int                 p;
  int                 mpiret;
  uint64_t            send[2];
  uint64_t           *gather;
#endif

  P4EST_ASSERT (p4est_is_valid (p4est));

  sc_array_init (&checkarray, 4);
  crc = adler32 (0, Z_NULL, 0);
  ssum = 0;
  for (nt = p4est->first_local_tree; nt <= p4est->last_local_tree; ++nt) {
    tree = p4est_array_index_topidx (p4est->trees, nt);
    treecrc = p4est_quadrant_checksum (&tree->quadrants, &checkarray, 0);
    scount = 4 * checkarray.elem_count;
    ssum += scount;
    crc = adler32_combine (crc, treecrc, (z_off_t) scount);
  }
  sc_array_reset (&checkarray);
  P4EST_ASSERT ((p4est_locidx_t) ssum ==
                p4est->local_num_quadrants * 4 * (P4EST_DIM + 1));

#ifdef P4EST_MPI
  send[0] = (uint64_t) crc;
  send[1] = (uint64_t) ssum;
  gather = NULL;
  if (p4est->mpirank == 0) {
    gather = P4EST_ALLOC (uint64_t, 2 * p4est->mpisize);
    if (p4est->mpicomm == MPI_COMM_NULL)
      memcpy (gather, send, 2 * sizeof (uint64_t));
  }
  if (p4est->mpicomm != MPI_COMM_NULL) {
    mpiret = MPI_Gather (send, 2, MPI_UNSIGNED_LONG_LONG,
                         gather, 2, MPI_UNSIGNED_LONG_LONG, 0,
                         p4est->mpicomm);
    SC_CHECK_MPI (mpiret);
  }

  crc = 0;
  if (p4est->mpirank == 0) {
    crc = gather[0];
    for (p = 1; p < p4est->mpisize; ++p) {
      crc = adler32_combine (crc, (uLong) gather[2 * p + 0],
                             (z_off_t) gather[2 * p + 1]);
    }
    P4EST_FREE (gather);
  }
#endif

  return (unsigned) crc;
}

/* EOF p4est.c */
