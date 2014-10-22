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

/********************************************************************
 *                          IMPORTANT NOTE                          *
 *                                                                  *
 * The p4est_points functionality depends on sc/src/sc_sort.        *
 * That parallel bitonic sort is still buggy (see sc/bugs).         *
 * If you want to use this code you have to fix the sort first.     *
 ********************************************************************/

#ifdef P4_TO_P8
#include <p8est_algorithms.h>
#include <p8est_bits.h>
#include <p8est_communication.h>
#include <p8est_extended.h>
#include <p8est_points.h>
#else
#include <p4est_algorithms.h>
#include <p4est_bits.h>
#include <p4est_communication.h>
#include <p4est_extended.h>
#include <p4est_points.h>
#endif /* !P4_TO_P8 */
#include <sc_allgather.h>
#include <sc_sort.h>

typedef struct
{
  p4est_quadrant_t   *points;
  p4est_locidx_t      num_points, max_points, current;
  int                 maxlevel;
}
p4est_points_state_t;

static void
p4est_points_init (p4est_t * p4est, p4est_topidx_t which_tree,
                   p4est_quadrant_t * quadrant)
{
  p4est_points_state_t *s = (p4est_points_state_t *) p4est->user_pointer;
  p4est_locidx_t     *qdata = (p4est_locidx_t *) quadrant->p.user_data;
  p4est_quadrant_t   *p;

  qdata[0] = s->current;
  while (s->current < s->num_points) {
    p = s->points + s->current;
    P4EST_ASSERT (p->p.which_tree >= which_tree);
    if (p->p.which_tree > which_tree) {
      break;
    }
    P4EST_ASSERT (p4est_quadrant_is_node (p, 1));
    P4EST_ASSERT (p4est_quadrant_compare (quadrant, p) < 0);
    if (!p4est_quadrant_contains_node (quadrant, p)) {
      break;
    }
    ++s->current;
  }
  qdata[1] = s->current;
}

static int
p4est_points_refine (p4est_t * p4est, p4est_topidx_t which_tree,
                     p4est_quadrant_t * quadrant)
{
  p4est_points_state_t *s = (p4est_points_state_t *) p4est->user_pointer;
  p4est_locidx_t     *qdata = (p4est_locidx_t *) quadrant->p.user_data;

  P4EST_ASSERT (s->max_points >= 0);
  if (qdata[1] > qdata[0] + s->max_points) {
    s->current = qdata[0];
    return 1;
  }
  else {
    return 0;
  }
}

p4est_t            *
p4est_new_points (sc_MPI_Comm mpicomm, p4est_connectivity_t * connectivity,
                  int maxlevel, p4est_quadrant_t * points,
                  p4est_locidx_t num_points, p4est_locidx_t max_points,
                  size_t data_size, p4est_init_t init_fn, void *user_pointer)
{
  int                 mpiret;
  int                 num_procs, rank;
  int                 i, isizet;
  size_t              lcount;
  size_t             *nmemb;
#ifdef P4EST_ENABLE_DEBUG
  size_t              zz;
#endif
  p4est_topidx_t      jt, num_trees;
  p4est_topidx_t      first_tree, last_tree, next_tree;
  p4est_quadrant_t   *first_quad, *next_quad, *quad;
  p4est_quadrant_t    a, b, c, f, l, n;
  p4est_tree_t       *tree;
  p4est_t            *p4est;
  p4est_points_state_t ppstate;

  P4EST_GLOBAL_PRODUCTIONF ("Into " P4EST_STRING
                            "_new_points with max level %d max points %lld\n",
                            maxlevel, (long long) max_points);
  p4est_log_indent_push ();
  P4EST_ASSERT (p4est_connectivity_is_valid (connectivity));
  P4EST_ASSERT (max_points >= -1);

  /* retrieve MPI information */
  mpiret = sc_MPI_Comm_size (mpicomm, &num_procs);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (mpicomm, &rank);
  SC_CHECK_MPI (mpiret);

  /* This implementation runs in O(P/p * maxlevel)
   * with P the total number of points, p the number of processors.
   * Two optimizations are possible:
   * 1. At startup remove points that lead to duplicate quadrants.
   * 2. Use complete_region between successive points instead of
   *    the call to refine. This should give O(N/p) * maxlevel
   *    with N the total number of quadrants.
   */

  /* parallel sort the incoming points */
  lcount = (size_t) num_points;
  nmemb = P4EST_ALLOC_ZERO (size_t, num_procs);
  isizet = (int) sizeof (size_t);
  mpiret = sc_MPI_Allgather (&lcount, isizet, sc_MPI_BYTE,
                             nmemb, isizet, sc_MPI_BYTE, mpicomm);
  SC_CHECK_MPI (mpiret);
  sc_psort (mpicomm, points, nmemb, sizeof (p4est_quadrant_t),
            p4est_quadrant_compare_piggy);
  P4EST_FREE (nmemb);
#ifdef P4EST_ENABLE_DEBUG
  first_quad = points;
  for (zz = 1; zz < lcount; ++zz) {
    next_quad = points + zz;
    P4EST_ASSERT (p4est_quadrant_compare_piggy (first_quad, next_quad) <= 0);
    first_quad = next_quad;
  }
#endif

  /* create the p4est */
  p4est = P4EST_ALLOC_ZERO (p4est_t, 1);
  ppstate.points = points;
  ppstate.num_points = num_points;
  ppstate.max_points = max_points;
  ppstate.current = 0;
  ppstate.maxlevel = maxlevel;

  /* assign some data members */
  p4est->mpicomm = mpicomm;
  p4est->mpisize = num_procs;
  p4est->mpirank = rank;
  p4est->data_size = 2 * sizeof (p4est_locidx_t);       /* temporary */
  p4est->user_pointer = &ppstate;
  p4est->connectivity = connectivity;
  num_trees = connectivity->num_trees;

  /* allocate memory pools */
  p4est->user_data_pool = sc_mempool_new (p4est->data_size);
  p4est->quadrant_pool = sc_mempool_new (sizeof (p4est_quadrant_t));

  P4EST_GLOBAL_PRODUCTIONF ("New " P4EST_STRING
                            " with %lld trees on %d processors\n",
                            (long long) num_trees, num_procs);

  /* allocate trees */
  p4est->trees = sc_array_new (sizeof (p4est_tree_t));
  sc_array_resize (p4est->trees, num_trees);
  for (jt = 0; jt < num_trees; ++jt) {
    tree = p4est_tree_array_index (p4est->trees, jt);
    sc_array_init (&tree->quadrants, sizeof (p4est_quadrant_t));
    P4EST_QUADRANT_INIT (&tree->first_desc);
    P4EST_QUADRANT_INIT (&tree->last_desc);
    tree->quadrants_offset = 0;
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

  /* create point based partition */
  P4EST_QUADRANT_INIT (&f);
  p4est->global_first_position =
    P4EST_ALLOC_ZERO (p4est_quadrant_t, num_procs + 1);
  if (num_points == 0) {
    P4EST_VERBOSE ("Empty processor");
    first_tree = p4est->first_local_tree = -1;
    first_quad = NULL;
  }
  else {
    /* we are probably not empty */
    if (rank == 0) {
      first_tree = p4est->first_local_tree = 0;
      p4est_quadrant_set_morton (&f, maxlevel, 0);
    }
    else {
      first_tree = p4est->first_local_tree = points->p.which_tree;
      p4est_node_to_quadrant (points, maxlevel, &f);
    }
    first_quad = &f;
  }
  last_tree = p4est->last_local_tree = -2;
  p4est_comm_global_partition (p4est, first_quad);
  first_quad = p4est->global_first_position + rank;
  next_quad = p4est->global_first_position + (rank + 1);
  next_tree = next_quad->p.which_tree;
  if (first_tree >= 0 &&
      p4est_quadrant_is_equal (first_quad, next_quad) &&
      first_quad->p.which_tree == next_quad->p.which_tree) {
    /* if all our points are consumed by the next processor we are empty */
    first_tree = p4est->first_local_tree = -1;
  }
  if (first_tree >= 0) {
    /* we are definitely not empty */
    if (next_quad->x == 0 && next_quad->y == 0
#ifdef P4_TO_P8
        && next_quad->z == 0
#endif
      ) {
      last_tree = p4est->last_local_tree = next_tree - 1;
    }
    else {
      last_tree = p4est->last_local_tree = next_tree;
    }
    P4EST_ASSERT (first_tree <= last_tree);
  }

  /* fill the local trees */
  P4EST_QUADRANT_INIT (&a);
  P4EST_QUADRANT_INIT (&b);
  P4EST_QUADRANT_INIT (&c);
  P4EST_QUADRANT_INIT (&l);
  n = *next_quad;
  n.level = (int8_t) maxlevel;
  for (jt = first_tree; jt <= last_tree; ++jt) {
    int                 onlyone = 0;
    int                 includeb = 0;

    tree = p4est_tree_array_index (p4est->trees, jt);

    /* determine first local quadrant of this tree */
    if (jt == first_tree) {
      a = *first_quad;
      a.level = (int8_t) maxlevel;
      first_quad = next_quad = NULL;    /* free to use further down */
      P4EST_ASSERT (p4est_quadrant_is_valid (&a));
    }
    else {
      p4est_quadrant_set_morton (&a, maxlevel, 0);
      P4EST_ASSERT (jt < next_tree || p4est_quadrant_compare (&a, &n) < 0);
    }

    /* enlarge first local quadrant if possible */
    if (jt < next_tree) {
      while (p4est_quadrant_child_id (&a) == 0 && a.level > 0) {
        p4est_quadrant_parent (&a, &a);
      }
      P4EST_ASSERT (jt == first_tree || a.level == 0);
    }
    else {
      for (c = a; p4est_quadrant_child_id (&c) == 0; a = c) {
        p4est_quadrant_parent (&c, &c);
        p4est_quadrant_last_descendant (&c, &l, maxlevel);
        if (p4est_quadrant_compare (&l, &n) >= 0) {
          break;
        }
      }
      P4EST_ASSERT (a.level > 0);
      P4EST_ASSERT ((p4est_quadrant_last_descendant (&a, &l, maxlevel),
                     p4est_quadrant_compare (&l, &n) < 0));
    }
    p4est_quadrant_first_descendant (&a, &tree->first_desc, P4EST_QMAXLEVEL);

    /* determine largest possible last quadrant of this tree */
    if (jt < next_tree) {
      p4est_quadrant_last_descendant (&a, &l, maxlevel);
      p4est_quadrant_set_morton (&b, 0, 0);
      p4est_quadrant_last_descendant (&b, &b, maxlevel);
      if (p4est_quadrant_is_equal (&l, &b)) {
        onlyone = 1;
      }
      else {
        includeb = 1;
        for (c = b; p4est_quadrant_child_id (&c) == P4EST_CHILDREN - 1; b = c) {
          p4est_quadrant_parent (&c, &c);
          p4est_quadrant_first_descendant (&c, &f, maxlevel);
          if (p4est_quadrant_compare (&l, &f) >= 0) {
            break;
          }
        }
      }
    }
    else {
      b = n;
    }

    /* create a complete tree */
    if (onlyone) {
      quad = p4est_quadrant_array_push (&tree->quadrants);
      *quad = a;
      p4est_quadrant_init_data (p4est, jt, quad, p4est_points_init);
      tree->maxlevel = a.level;
      ++tree->quadrants_per_level[a.level];
    }
    else {
      p4est_complete_region (p4est, &a, 1, &b, includeb,
                             tree, jt, p4est_points_init);
      quad = p4est_quadrant_array_index (&tree->quadrants,
                                         tree->quadrants.elem_count - 1);
    }
    tree->quadrants_offset = p4est->local_num_quadrants;
    p4est->local_num_quadrants += tree->quadrants.elem_count;
    p4est_quadrant_last_descendant (quad, &tree->last_desc, P4EST_QMAXLEVEL);

    /* verification */
#ifdef P4EST_ENABLE_DEBUG
    first_quad = p4est_quadrant_array_index (&tree->quadrants, 0);
    for (zz = 1; zz < tree->quadrants.elem_count; ++zz) {
      next_quad = p4est_quadrant_array_index (&tree->quadrants, zz);
      P4EST_ASSERT (((p4est_locidx_t *) first_quad->p.user_data)[1] ==
                    ((p4est_locidx_t *) next_quad->p.user_data)[0]);
      first_quad = next_quad;
    }
#endif
  }
  if (last_tree >= 0) {
    for (; jt < num_trees; ++jt) {
      tree = p4est_tree_array_index (p4est->trees, jt);
      tree->quadrants_offset = p4est->local_num_quadrants;
    }
  }

  /* compute some member variables */
  p4est->global_first_quadrant = P4EST_ALLOC (p4est_gloidx_t, num_procs + 1);
  p4est_comm_count_quadrants (p4est);

  /* print more statistics */
  P4EST_VERBOSEF ("total local quadrants %lld\n",
                  (long long) p4est->local_num_quadrants);

  P4EST_ASSERT (p4est_is_valid (p4est));
  p4est_log_indent_pop ();
  P4EST_GLOBAL_PRODUCTIONF ("Done " P4EST_STRING
                            "_new_points with %lld total quadrants\n",
                            (long long) p4est->global_num_quadrants);

  /* refine to have one point per quadrant */
  if (max_points >= 0) {
    p4est_refine_ext (p4est, 1, maxlevel, p4est_points_refine,
                      p4est_points_init, NULL);
#ifdef P4EST_ENABLE_DEBUG
    for (jt = first_tree; jt <= last_tree; ++jt) {
      tree = p4est_tree_array_index (p4est->trees, jt);
      first_quad = p4est_quadrant_array_index (&tree->quadrants, 0);
      for (zz = 1; zz < tree->quadrants.elem_count; ++zz) {
        next_quad = p4est_quadrant_array_index (&tree->quadrants, zz);
        P4EST_ASSERT (((p4est_locidx_t *) first_quad->p.user_data)[1] ==
                      ((p4est_locidx_t *) next_quad->p.user_data)[0]);
        first_quad = next_quad;
      }
    }
#endif
  }

  /* initialize user pointer and data size */
  p4est_reset_data (p4est, data_size, init_fn, user_pointer);

  return p4est;
}
