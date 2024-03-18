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

#ifndef P4_TO_P8
#include <p4est_bits.h>
#include <p4est_communication.h>
#include <p4est_search.h>
#else
#include <p8est_bits.h>
#include <p8est_communication.h>
#include <p8est_search.h>
#endif
#include <sc_search.h>

#define P4EST_COMM_IS_EMPTY_GFQ_GFP(gfq, gfp, num_procs, p) \
                (((gfq) != NULL) ? \
                p4est_comm_is_empty_gfq ((gfq), (num_procs), (p)) :\
                p4est_comm_is_empty_gfp ((gfp), (num_procs), (p)))

/** A callback function that describes the search window.
 *  The idea is to define the type of an array entry as type 1, if
 *  my_begin <= array[i], my_end > array[i] and as type 2, if the entry
 *  is not of type 1 and my_end <= array[i]. The remaining cases are
 *  indicated by type 0. This function is passed into sc_array_split by
 *  \ref p4est_find_partition. Note that this function, as well as
 *  p4est_find_partition, is dimension-independent; still we duplicate
 *  it in 3D in the usual way.
 */
static              size_t
type_fn_global_quad_index (sc_array_t * array, size_t index, void *data_array)
{
  p4est_gloidx_t     *my_begin_end;

  my_begin_end = (p4est_gloidx_t *) data_array;
  if (*(p4est_gloidx_t *) sc_array_index (array, index) < my_begin_end[0])
    return 0;
  if ((my_begin_end[0] <= *(p4est_gloidx_t *) sc_array_index (array, index))
      && (my_begin_end[1] >
          *(p4est_gloidx_t *) sc_array_index (array, index)))
    return 1;
  if (my_begin_end[1] <= *(p4est_gloidx_t *) sc_array_index (array, index))
    return 2;
  return 0;
}

void
p4est_find_partition (const int num_procs, p4est_gloidx_t * search_in,
                      p4est_gloidx_t my_begin, p4est_gloidx_t my_end,
                      p4est_gloidx_t * begin, p4est_gloidx_t * end)
{
  sc_array_t          offsets[1];
  sc_array_t          view[1];
  p4est_gloidx_t      my_begin_end[2];

  P4EST_ASSERT (my_begin <= my_end);

  sc_array_init (offsets, sizeof (size_t));
  sc_array_init_data (view, search_in, sizeof (p4est_gloidx_t),
                      (size_t) num_procs);

  my_begin_end[0] = my_begin;
  my_begin_end[1] = my_end;

  sc_array_split (view, offsets, 3, type_fn_global_quad_index, my_begin_end);

  *begin = *(p4est_gloidx_t *) sc_array_index (offsets, 1);
  *end = *(p4est_gloidx_t *) sc_array_index (offsets, 2);

  sc_array_reset (offsets);
  sc_array_reset (view);
}

ssize_t
p4est_find_lower_bound (sc_array_t * array,
                        const p4est_quadrant_t * q, size_t guess)
{
  int                 comp;
  size_t              count;
  size_t              quad_low, quad_high;
  p4est_quadrant_t   *cur;

  count = array->elem_count;
  if (count == 0)
    return -1;

  quad_low = 0;
  quad_high = count - 1;

  for (;;) {
    P4EST_ASSERT (quad_low <= quad_high);
    P4EST_ASSERT (quad_low < count && quad_high < count);
    P4EST_ASSERT (quad_low <= guess && guess <= quad_high);

    /* compare two quadrants */
    cur = p4est_quadrant_array_index (array, guess);
    comp = p4est_quadrant_compare (q, cur);

    /* check if guess is higher or equal q and there's room below it */
    if (comp <= 0 && (guess > 0 && p4est_quadrant_compare (q, cur - 1) <= 0)) {
      quad_high = guess - 1;
      guess = (quad_low + quad_high + 1) / 2;
      continue;
    }

    /* check if guess is lower than q */
    if (comp > 0) {
      quad_low = guess + 1;
      if (quad_low > quad_high)
        return -1;

      guess = (quad_low + quad_high) / 2;
      continue;
    }

    /* otherwise guess is the correct quadrant */
    break;
  }

  return (ssize_t) guess;
}

ssize_t
p4est_find_higher_bound (sc_array_t * array,
                         const p4est_quadrant_t * q, size_t guess)
{
  int                 comp;
  size_t              count;
  size_t              quad_low, quad_high;
  p4est_quadrant_t   *cur;

  count = array->elem_count;
  if (count == 0)
    return -1;

  quad_low = 0;
  quad_high = count - 1;

  for (;;) {
    P4EST_ASSERT (quad_low <= quad_high);
    P4EST_ASSERT (quad_low < count && quad_high < count);
    P4EST_ASSERT (quad_low <= guess && guess <= quad_high);

    /* compare two quadrants */
    cur = p4est_quadrant_array_index (array, guess);
    comp = p4est_quadrant_compare (cur, q);

    /* check if guess is lower or equal q and there's room above it */
    if (comp <= 0 &&
        (guess < count - 1 && p4est_quadrant_compare (cur + 1, q) <= 0)) {
      quad_low = guess + 1;
      guess = (quad_low + quad_high) / 2;
      continue;
    }

    /* check if guess is higher than q */
    if (comp > 0) {
      if (guess == 0)
        return -1;

      quad_high = guess - 1;
      if (quad_high < quad_low)
        return -1;

      guess = (quad_low + quad_high + 1) / 2;
      continue;
    }

    /* otherwise guess is the correct quadrant */
    break;
  }

  return (ssize_t) guess;
}

p4est_quadrant_t   *
p4est_find_quadrant_cumulative (p4est_t * p4est, p4est_locidx_t cumulative_id,
                                p4est_topidx_t * which_tree,
                                p4est_locidx_t * quadrant_id)
{
  int                 the_quadrant_id;
  p4est_topidx_t      low_tree, high_tree, guess_tree;
  p4est_tree_t       *tree;

  P4EST_ASSERT (0 <= cumulative_id &&
                cumulative_id < p4est->local_num_quadrants);

  low_tree = p4est->first_local_tree;
  high_tree = p4est->last_local_tree;
  if (which_tree != NULL && *which_tree != -1) {
    guess_tree = *which_tree;
    P4EST_ASSERT (0 <= guess_tree &&
                  guess_tree < p4est->connectivity->num_trees);
  }
  else {
    guess_tree = (low_tree + high_tree) / 2;
  }
  for (;;) {
    P4EST_ASSERT (p4est->first_local_tree <= low_tree);
    P4EST_ASSERT (high_tree <= p4est->last_local_tree);
    P4EST_ASSERT (low_tree <= guess_tree && guess_tree <= high_tree);

    tree = p4est_tree_array_index (p4est->trees, guess_tree);
    if (cumulative_id < tree->quadrants_offset) {
      high_tree = guess_tree - 1;
    }
    else if (cumulative_id >= tree->quadrants_offset +
             (p4est_locidx_t) tree->quadrants.elem_count) {
      low_tree = guess_tree + 1;
    }
    else {
      the_quadrant_id = cumulative_id - tree->quadrants_offset;
      P4EST_ASSERT (0 <= the_quadrant_id);

      if (which_tree != NULL) {
        *which_tree = guess_tree;
      }
      if (quadrant_id != NULL) {
        *quadrant_id = the_quadrant_id;
      }
      return p4est_quadrant_array_index (&tree->quadrants,
                                         (size_t) the_quadrant_id);
    }
    guess_tree = (low_tree + high_tree) / 2;
  }
}

static              size_t
p4est_array_split_ancestor_id (sc_array_t * array, size_t zindex, void *data)
{
  int                *levelp = (int *) data;
  p4est_quadrant_t   *q = p4est_quadrant_array_index (array, zindex);

  return ((size_t) p4est_quadrant_ancestor_id (q, *levelp));
}

void
p4est_split_array (sc_array_t * array, int level, size_t indices[])
{
  size_t              count = array->elem_count;
  sc_array_t          view;
#ifdef P4EST_ENABLE_DEBUG
  p4est_quadrant_t   *test1, test2;
  p4est_quadrant_t   *cur;
#endif

  P4EST_ASSERT (0 <= level && level < P4EST_QMAXLEVEL);
  /** If empty, return all zeroes */
  if (count == 0) {
    indices[0] = indices[1] = indices[2] = indices[3] = indices[4] =
#ifdef P4_TO_P8
      indices[5] = indices[6] = indices[7] = indices[8] =
#endif
      0;
    return;
  }

  P4EST_ASSERT (sc_array_is_sorted (array, p4est_quadrant_compare));
#ifdef P4EST_ENABLE_DEBUG
  cur = p4est_quadrant_array_index (array, 0);
  P4EST_ASSERT ((int) cur->level > level);
  test1 = p4est_quadrant_array_index (array, count - 1);
  P4EST_ASSERT ((int) test1->level > level);
  p4est_nearest_common_ancestor (cur, test1, &test2);
  P4EST_ASSERT ((int) test2.level >= level);
#endif

  sc_array_init_data (&view, indices, sizeof (size_t), P4EST_CHILDREN + 1);
  level++;
  sc_array_split (array, &view, P4EST_CHILDREN, p4est_array_split_ancestor_id,
                  &level);
}

/** If we suppose a range of quadrants touches a corner of a tree, then it must
 * also touch the faces (and edges) that touch that corner.
 */
#ifndef P4_TO_P8
/* *INDENT-OFF* */
static int32_t p4est_corner_boundaries[4] =
{             /*                           |corners | faces */
  0x00000015, /* 0000 0000 0000 0000 0000 0000| 0001| 0101  */
  0x00000026, /* 0000 0000 0000 0000 0000 0000| 0010| 0110  */
  0x00000049, /* 0000 0000 0000 0000 0000 0000| 0100| 1001  */
  0x0000008a  /* 0000 0000 0000 0000 0000 0000| 1000| 1010  */
};
/* *INDENT-ON* */
static int32_t      p4est_all_boundaries = 0x000000ff;
#else
/* *INDENT-OFF* */
static int32_t p4est_corner_boundaries[8] =
{             /*        |corners   |edges          |faces   */
  0x00044455, /* 0000 00|00 0000 01|00 0100 0100 01|01 0101 */
  0x00088856, /* 0000 00|00 0000 10|00 1000 1000 01|01 0110 */
  0x00110499, /* 0000 00|00 0001 00|01 0000 0100 10|01 1001 */
  0x0022089a, /* 0000 00|00 0010 00|10 0000 1000 10|01 1010 */
  0x00405125, /* 0000 00|00 0100 00|00 0101 0001 00|10 0101 */
  0x0080a126, /* 0000 00|00 1000 00|00 1010 0001 00|10 0110 */
  0x01011229, /* 0000 00|01 0000 00|01 0001 0010 00|10 1001 */
  0x0202222a  /* 0000 00|10 0000 00|10 0010 0010 00|10 1010 */
};
/* *INDENT-ON* */
static int32_t      p4est_all_boundaries = 0x03ffffff;
#endif

static              int32_t
p4est_limit_boundaries (p4est_quadrant_t * q, int dir, int limit,
                        int last_level, int level, int32_t touch,
                        int32_t mask)
{
  int                 cid;
  int32_t             next;

  P4EST_ASSERT (q->level == P4EST_QMAXLEVEL);
  P4EST_ASSERT (level <= P4EST_QMAXLEVEL);
  P4EST_ASSERT (level <= last_level);
  if ((mask & ~touch) == 0) {
    return touch;
  }
  cid = p4est_quadrant_ancestor_id (q, level);
  next = p4est_corner_boundaries[cid] & mask;
  cid += dir;
  while (cid != limit) {
    touch |= (p4est_corner_boundaries[cid] & mask);
    cid += dir;
  }
  if (level == last_level) {
    return (touch | next);
  }
  return p4est_limit_boundaries (q, dir, limit, last_level, level + 1, touch,
                                 next);
}

static              int32_t
p4est_range_boundaries (p4est_quadrant_t * lq, p4est_quadrant_t * uq,
                        int alevel, int level, int32_t mask)
{
  int                 i, lcid, ucid, cid;
  int32_t             lnext, unext, touch;
  p4est_qcoord_t      x, y, a;
#ifdef P4_TO_P8
  p4est_qcoord_t      z;
#endif
  const p4est_qcoord_t shift = P4EST_QUADRANT_LEN (P4EST_QMAXLEVEL);
  int                 count;
  int                 last_level;

  P4EST_ASSERT (level <= alevel + 1);

  if (mask == 0) {
    return 0;
  }
  if (level == alevel + 1) {
    lcid = p4est_quadrant_ancestor_id (lq, level);
    ucid = p4est_quadrant_ancestor_id (uq, level);
    P4EST_ASSERT (lcid < ucid);
    lnext = (p4est_corner_boundaries[lcid] & mask);
    unext = (p4est_corner_boundaries[ucid] & mask);
    touch = 0;
    for (i = lcid + 1; i < ucid; i++) {
      touch |= (p4est_corner_boundaries[i] & mask);
    }

    cid = p4est_quadrant_child_id (lq);
    x = lq->x + ((cid & 1) ? shift : 0);
    y = lq->y + (((cid >> 1) & 1) ? shift : 0);
#ifdef P4_TO_P8
    z = lq->z + ((cid >> 2) ? shift : 0);
#endif
    a = ~(x | y
#ifdef P4_TO_P8
          | z
#endif
      );
    count = 0;
    while ((a & ((p4est_qcoord_t) 1)) && count <= P4EST_MAXLEVEL) {
      a >>= 1;
      count++;
    }
    last_level = (P4EST_MAXLEVEL - count) + 1;
    if (last_level <= level) {
      touch |= lnext;
    }
    else {
      P4EST_ASSERT (last_level <= P4EST_QMAXLEVEL);
      touch |= p4est_limit_boundaries (lq, 1, P4EST_CHILDREN, last_level,
                                       level + 1, touch, lnext);
    }

    cid = p4est_quadrant_child_id (uq);
    x = uq->x + ((cid & 1) ? shift : 0);
    y = uq->y + (((cid >> 1) & 1) ? shift : 0);
#ifdef P4_TO_P8
    z = uq->z + ((cid >> 2) ? shift : 0);
#endif
    a = ~(x | y
#ifdef P4_TO_P8
          | z
#endif
      );
    count = 0;
    while ((a & ((p4est_qcoord_t) 1)) && count <= P4EST_MAXLEVEL) {
      a >>= 1;
      count++;
    }
    last_level = (P4EST_MAXLEVEL - count) + 1;
    if (last_level <= level) {
      touch |= unext;
    }
    else {
      P4EST_ASSERT (last_level <= P4EST_QMAXLEVEL);
      touch |= p4est_limit_boundaries (uq, -1, -1, last_level, level + 1,
                                       touch, unext);
    }

    return touch;
  }
  lcid = p4est_quadrant_ancestor_id (lq, level);
  P4EST_ASSERT (p4est_quadrant_ancestor_id (uq, level) == lcid);
  return p4est_range_boundaries (lq, uq, alevel, level + 1,
                                 (p4est_corner_boundaries[lcid] & mask));
}

int32_t
p4est_find_range_boundaries (p4est_quadrant_t * lq, p4est_quadrant_t * uq,
                             int level, int faces[],
#ifdef P4_TO_P8
                             int edges[],
#endif
                             int corners[])
{
  int                 i;
  p4est_quadrant_t    a;
  int                 alevel;
  int32_t             touch;
  int32_t             mask = 0x00000001;
  p4est_qcoord_t      x, y, all;
#ifdef P4_TO_P8
  p4est_qcoord_t      z;
#endif
  const p4est_qcoord_t shift = P4EST_QUADRANT_LEN (P4EST_QMAXLEVEL);
  int                 count;
  int                 last_level;
  int                 cid;

  P4EST_ASSERT (level >= 0 && level <= P4EST_QMAXLEVEL);
  if ((lq == NULL && uq == NULL) || level == P4EST_QMAXLEVEL) {
    touch = p4est_all_boundaries;
    goto find_range_boundaries_exit;
  }

  if (lq == NULL) {
    P4EST_ASSERT (uq->level == P4EST_QMAXLEVEL);

    cid = p4est_quadrant_child_id (uq);
    x = uq->x + ((cid & 1) ? shift : 0);
    y = uq->y + (((cid >> 1) & 1) ? shift : 0);
#ifdef P4_TO_P8
    z = uq->z + ((cid >> 2) ? shift : 0);
#endif
    all = ~(x | y
#ifdef P4_TO_P8
            | z
#endif
      );
    count = 0;
    while ((all & ((p4est_qcoord_t) 1)) && count <= P4EST_MAXLEVEL) {
      all >>= 1;
      count++;
    }
    last_level = (P4EST_MAXLEVEL - count) + 1;
    last_level = (last_level <= level) ? level + 1 : last_level;

    P4EST_ASSERT (last_level <= P4EST_QMAXLEVEL);

    touch = p4est_limit_boundaries (uq, -1, -1, last_level, level + 1, 0,
                                    p4est_all_boundaries);
  }
  else if (uq == NULL) {
    P4EST_ASSERT (lq->level == P4EST_QMAXLEVEL);

    cid = p4est_quadrant_child_id (lq);
    x = lq->x + ((cid & 1) ? shift : 0);
    y = lq->y + (((cid >> 1) & 1) ? shift : 0);
#ifdef P4_TO_P8
    z = lq->z + ((cid >> 2) ? shift : 0);
#endif
    all = ~(x | y
#ifdef P4_TO_P8
            | z
#endif
      );
    count = 0;
    while ((all & ((p4est_qcoord_t) 1)) && count <= P4EST_MAXLEVEL) {
      all >>= 1;
      count++;
    }
    last_level = (P4EST_MAXLEVEL - count) + 1;
    last_level = (last_level <= level) ? level + 1 : last_level;

    P4EST_ASSERT (last_level <= P4EST_QMAXLEVEL);

    touch = p4est_limit_boundaries (lq, 1, P4EST_CHILDREN, last_level,
                                    level + 1, 0, p4est_all_boundaries);
  }
  else {
    P4EST_ASSERT (uq->level == P4EST_QMAXLEVEL);
    P4EST_ASSERT (lq->level == P4EST_QMAXLEVEL);
    p4est_nearest_common_ancestor (lq, uq, &a);
    alevel = (int) a.level;
    P4EST_ASSERT (alevel >= level);
    touch = p4est_range_boundaries (lq, uq, alevel, level + 1,
                                    p4est_all_boundaries);
  }

find_range_boundaries_exit:
  if (faces != NULL) {
    for (i = 0; i < P4EST_FACES; i++) {
      faces[i] = (touch & mask) != 0;
      mask <<= 1;
    }
  }
  else {
    mask <<= P4EST_FACES;
  }
#ifdef P4_TO_P8
  if (edges != NULL) {
    for (i = 0; i < P8EST_EDGES; i++) {
      edges[i] = (touch & mask) != 0;
      mask <<= 1;
    }
  }
  else {
    mask <<= P8EST_EDGES;
  }
#endif
  if (corners != NULL) {
    for (i = 0; i < P4EST_CHILDREN; i++) {
      corners[i] = (touch & mask) != 0;
      mask <<= 1;
    }
  }

  return touch;
}

/** This recursion context saves on the number of parameters passed. */
typedef struct p4est_local_recursion
{
  p4est_t            *p4est;            /**< Forest being traversed. */
  p4est_topidx_t      which_tree;       /**< Current tree number. */
  int                 call_post;        /**< Boolean to call quadrant twice. */
  int                 skip;             /**< Boolean to avoid skipping levels in parallel, if desired*/
  p4est_search_reorder_t children_fn;   /**< Reorder children if not NULL. */
  p4est_search_local_t quadrant_fn;     /**< The quadrant callback for backwards compatibility, if any. */
  p4est_search_local_t pre_quadrant_fn; /**< The pre recursion quadrant callback, if any. */
  p4est_search_local_t post_quadrant_fn;/**< The post recursion quadrant callback, if any. */
  p4est_search_local_t point_fn;        /**< The point callback, if any. */
  sc_array_t         *points;           /**< Array of points to search. */
}
p4est_local_recursion_t;

/* The recursion may overwrite the \a quadrant input argument contents. */
static void
p4est_local_recursion (const p4est_local_recursion_t * rec,
                       p4est_quadrant_t * quadrant,
                       sc_array_t * quadrants, sc_array_t * actives)
{
  int                 i;
  int                 is_leaf, is_match;
  int                 level;
  size_t              qcount, act_count;
  size_t              zz, *pz, *qz;
  size_t              split[P4EST_CHILDREN + 1];
  p4est_locidx_t      local_num;
  p4est_quadrant_t   *q, *lq, child;
  sc_array_t          child_quadrants, child_actives, *chact;

  /*
   * Invariants of the recursion:
   * 1. quadrant is larger or equal in size than those in the array.
   * 2. quadrant is equal to or an ancestor of those in the array.
   */
  P4EST_ASSERT (rec != NULL);
  P4EST_ASSERT (rec->children_fn == NULL);
  P4EST_ASSERT (quadrant != NULL && quadrants != NULL);
  P4EST_ASSERT (quadrants->elem_size == sizeof (p4est_quadrant_t));
  P4EST_ASSERT (rec->skip);     /* search_local tries to skip always */
  P4EST_ASSERT (rec->pre_quadrant_fn == NULL);  /* unused here */
  P4EST_ASSERT (rec->post_quadrant_fn == NULL); /* unused here */

  qcount = quadrants->elem_count;

  /* As an optimization we pass a NULL actives array to every root. */
  if (rec->points != NULL && actives == NULL) {
    act_count = rec->points->elem_count;
  }
  else {
    P4EST_ASSERT ((rec->points == NULL) == (actives == NULL));
    P4EST_ASSERT (rec->points == NULL ||
                  actives->elem_count <= rec->points->elem_count);
    act_count = actives == NULL ? 0 : actives->elem_count;
  }

  /* return if there are no quadrants or active points */
  if (qcount == 0 || (rec->points != NULL && act_count == 0))
    return;

  /* determine leaf situation */
  q = p4est_quadrant_array_index (quadrants, 0);
  if (qcount > 1) {
    P4EST_ASSERT (p4est_quadrant_is_ancestor (quadrant, q));
    is_leaf = 0;
    local_num = -1;
    lq = p4est_quadrant_array_index (quadrants, quadrants->elem_count - 1);
    P4EST_ASSERT (!p4est_quadrant_is_equal (q, lq) &&
                  p4est_quadrant_is_ancestor (quadrant, lq));

    /* skip unnecessary intermediate levels if possible */
    level = (int) quadrant->level;
    if (p4est_quadrant_ancestor_id (q, level + 1) ==
        p4est_quadrant_ancestor_id (lq, level + 1)) {
      p4est_nearest_common_ancestor (q, lq, quadrant);
      P4EST_ASSERT (level < (int) quadrant->level);
      P4EST_ASSERT (p4est_quadrant_is_ancestor (quadrant, q));
      P4EST_ASSERT (p4est_quadrant_is_ancestor (quadrant, lq));
    }
  }
  else {
    p4est_locidx_t      offset;
    p4est_tree_t       *tree;

    P4EST_ASSERT (p4est_quadrant_is_equal (quadrant, q) ||
                  p4est_quadrant_is_ancestor (quadrant, q));
    is_leaf = 1;

    /* determine offset of quadrant in local forest */
    tree = p4est_tree_array_index (rec->p4est->trees, rec->which_tree);
    offset = (p4est_locidx_t) ((quadrants->array - tree->quadrants.array)
                               / sizeof (p4est_quadrant_t));
    P4EST_ASSERT (offset >= 0 &&
                  (size_t) offset < tree->quadrants.elem_count);
    local_num = tree->quadrants_offset + offset;

    /* skip unnecessary intermediate levels if possible */
    quadrant = q;
  }

  /* execute pre-quadrant callback if present, which may stop the recursion */
  if (rec->quadrant_fn != NULL &&
      !rec->quadrant_fn (rec->p4est, rec->which_tree,
                         quadrant, local_num, NULL)) {
    return;
  }

  /* check out points */
  if (rec->points == NULL) {
    /* we have called the callback already.  For leaves we are done */
    if (is_leaf) {
      return;
    }
    chact = NULL;
  }
  else {
    /* query callback for all points and return if none remain */
    chact = &child_actives;
    sc_array_init (chact, sizeof (size_t));
    for (zz = 0; zz < act_count; ++zz) {
      pz = actives == NULL ? &zz : (size_t *) sc_array_index (actives, zz);
      is_match = rec->point_fn (rec->p4est, rec->which_tree,
                                quadrant, local_num,
                                sc_array_index (rec->points, *pz));
      if (!is_leaf && is_match) {
        qz = (size_t *) sc_array_push (chact);
        *qz = *pz;
      }
    }

    /* call post-quadrant callback, which may also terminate the recursion */
    if (rec->call_post && rec->quadrant_fn != NULL &&
        !rec->quadrant_fn (rec->p4est, rec->which_tree,
                           quadrant, local_num, NULL)) {
      /* clears memory and will trigger the return below */
      sc_array_reset (chact);
    }

    if (chact->elem_count == 0) {
      /* with zero members there is no need to call sc_array_reset */
      return;
    }
  }

  /* leaf situation has returned above */
  P4EST_ASSERT (!is_leaf);
  P4EST_ASSERT (quadrant->level < P4EST_QMAXLEVEL);

  /* split quadrant array and run recursion */
  p4est_split_array (quadrants, (int) quadrant->level, split);
  for (i = 0; i < P4EST_CHILDREN; ++i) {
    if (split[i] < split[i + 1]) {
      sc_array_init_view (&child_quadrants, quadrants,
                          split[i], split[i + 1] - split[i]);
      p4est_quadrant_child (quadrant, &child, i);
      p4est_local_recursion (rec, &child, &child_quadrants, chact);
      sc_array_reset (&child_quadrants);
    }
  }
  if (chact != NULL) {
    sc_array_reset (chact);
  }
}

void
p4est_search_local (p4est_t * p4est,
                    int call_post, p4est_search_local_t quadrant_fn,
                    p4est_search_local_t point_fn, sc_array_t * points)
{
  p4est_topidx_t      jt;
  p4est_tree_t       *tree;
  p4est_quadrant_t    root;
  p4est_local_recursion_t srec, *rec = &srec;
  sc_array_t         *tquadrants;

  /* correct call convention? */
  P4EST_ASSERT (p4est != NULL);
  P4EST_ASSERT (points == NULL || point_fn != NULL);

  /* we do nothing if there is nothing we can do */
  if (quadrant_fn == NULL && points == NULL) {
    return;
  }

  /* set recursion context */
  rec->p4est = p4est;
  rec->which_tree = -1;
  rec->call_post = call_post;
  rec->children_fn = NULL;
  rec->quadrant_fn = quadrant_fn;
  rec->pre_quadrant_fn = NULL;
  rec->post_quadrant_fn = NULL;
  rec->point_fn = point_fn;
  rec->points = points;
  rec->skip = 1;
  for (jt = p4est->first_local_tree; jt <= p4est->last_local_tree; ++jt) {
    rec->which_tree = jt;

    /* grab complete tree quadrant array */
    tree = p4est_tree_array_index (p4est->trees, jt);
    tquadrants = &tree->quadrants;

    /* the recursion shrinks the search quadrant whenever possible */
    p4est_quadrant_set_morton (&root, 0, 0);
    p4est_local_recursion (rec, &root, tquadrants, NULL);
  }
}

/* The recursion may overwrite the \a quadrant input argument contents. */
static void
p4est_reorder_recursion (const p4est_local_recursion_t * rec,
                         p4est_quadrant_t * quadrant,
                         sc_array_t * quadrants, sc_array_t * actives)
{
  int                 i;
  int                 is_leaf;
  int                 conchildren;
  int                 level;
  int                 is_same;
  size_t              qcount, act_count;
  size_t              zz, *pz, *qz;
  size_t              split[P4EST_CHILDREN + 1];
  p4est_locidx_t      local_num;
  p4est_quadrant_t   *q;
  sc_array_t          child_quadrants, child_actives, *chact;
  sc_array_t          child_indices;

  /*
   * Invariants of the recursion:
   * 1. quadrant is larger or equal in size than those in the array.
   * 2. quadrant is equal to or an ancestor of those in the array.
   */
  P4EST_ASSERT (rec != NULL);
  P4EST_ASSERT (rec->quadrant_fn == NULL);      /* unused here */
  P4EST_ASSERT (quadrant != NULL && quadrants != NULL);
  P4EST_ASSERT (quadrants->elem_size == sizeof (p4est_quadrant_t));
  qcount = quadrants->elem_count;

  /* As an optimization we pass a NULL actives array to every root. */
  if (rec->points != NULL && actives == NULL) {
    P4EST_ASSERT (quadrant->level == 0);
    act_count = rec->points->elem_count;
  }
  else {
    P4EST_ASSERT ((rec->points == NULL) == (actives == NULL));
    P4EST_ASSERT (rec->points == NULL ||
                  actives->elem_count <= rec->points->elem_count);
    act_count = actives == NULL ? 0 : actives->elem_count;
  }

  /* return if there are no quadrants or active points */
  if (qcount == 0 || (rec->points != NULL && act_count == 0)) {
    return;
  }

  /* determine leaf situation: can the search quadrant be shrunk? */
  local_num = -1;
  is_leaf = is_same = 0;
  q = p4est_quadrant_array_index (quadrants, 0);
  if (qcount > 1) {
    /* There are multiple quadrants in the search window.
       This is definitely not a leaf situation. */

    if (rec->skip) {
      p4est_quadrant_t   *lq;

      /* access last leaf in the search window in addition to the first */
      P4EST_ASSERT (p4est_quadrant_is_ancestor (quadrant, q));
      lq = p4est_quadrant_array_index (quadrants, quadrants->elem_count - 1);
      P4EST_ASSERT (!p4est_quadrant_is_equal (q, lq) &&
                    p4est_quadrant_is_ancestor (quadrant, lq));

      /* skip unnecessary intermediate levels if possible */
      level = (int) quadrant->level;
      if (p4est_quadrant_ancestor_id (q, level + 1) ==
          p4est_quadrant_ancestor_id (lq, level + 1)) {
        p4est_nearest_common_ancestor (q, lq, quadrant);
        P4EST_ASSERT (level < (int) quadrant->level);
        P4EST_ASSERT (p4est_quadrant_is_ancestor (quadrant, q));
        P4EST_ASSERT (p4est_quadrant_is_ancestor (quadrant, lq));
      }
    }
  }
  else {
    P4EST_ASSERT (p4est_quadrant_is_equal (quadrant, q) ||
                  p4est_quadrant_is_ancestor (quadrant, q));

    /* q is the only leaf in the search quadrant.
       It can be reached directly if we skip all intermediate levels. */
    is_leaf = 1;
    is_same = p4est_quadrant_is_equal (quadrant, q);
    if (rec->skip && !is_same) {
      /* we are asked to optimize by skipping intermediate levels */
      is_same = 1;
    }
    if (is_same) {
      p4est_locidx_t      offset;
      p4est_tree_t       *tree;

      /* determine offset of leaf quadrant in local forest */
      tree = p4est_tree_array_index (rec->p4est->trees, rec->which_tree);
      offset = (p4est_locidx_t) ((quadrants->array - tree->quadrants.array)
                                 / sizeof (p4est_quadrant_t));
      P4EST_ASSERT (offset >= 0 &&
                    (size_t) offset < tree->quadrants.elem_count);
      local_num = tree->quadrants_offset + offset;

      /* make sure we pass the leaf itself to subsequent callbacks */
      quadrant = q;
    }
  }
  P4EST_ASSERT (!is_same || is_leaf);
  P4EST_ASSERT (!is_same == !p4est_quadrant_is_equal (quadrant, q));
  P4EST_ASSERT (is_same || quadrant->level < P4EST_QMAXLEVEL);

  /* execute pre-quadrant callback if present, which may stop the recursion */
  if (rec->pre_quadrant_fn != NULL && !rec->pre_quadrant_fn
      (rec->p4est, rec->which_tree, quadrant, local_num, NULL)) {
    return;
  }

  /* call point callback on remaining points, skip to post if none remain */
  chact = NULL;
  conchildren = 1;
  if (rec->points != NULL) {
    /* query callback for all points and skip to post if none remain */
    chact = &child_actives;
    sc_array_init (chact, sizeof (size_t));
    for (zz = 0; zz < act_count; ++zz) {
      pz = actives == NULL ? &zz : (size_t *) sc_array_index (actives, zz);
      if (rec->point_fn (rec->p4est, rec->which_tree, quadrant, local_num,
                         sc_array_index (rec->points, *pz)) && !is_same) {
        qz = (size_t *) sc_array_push (chact);
        *qz = *pz;
      }
    }
    if (chact->elem_count == 0) {
      /* without members (forced for any leaf) no need to call sc_array_reset */
      chact = NULL;
      conchildren = 0;
    }
  }

  /* execute search recursion, either when points remain or there never were any */
  if (conchildren && !is_leaf) {
    /* identify search children, figure out which are relevant */
    P4EST_ASSERT (qcount > 1);
    p4est_split_array (quadrants, (int) quadrant->level, split);
    sc_array_init (&child_indices, sizeof (p4est_topidx_t));
    for (i = 0; i < P4EST_CHILDREN; ++i) {
      if (split[i] < split[i + 1]) {
        *(p4est_topidx_t *) sc_array_push (&child_indices) = i;
      }
    }
    P4EST_ASSERT (child_indices.elem_count > 0);

    /* reorder/reduce search children, skip to post if callback returns false */
    if (rec->children_fn != NULL) {
      sc_array_t          childrenview;
      p4est_quadrant_t    children[P4EST_CHILDREN];
      p4est_quadrant_childrenv (quadrant, children);
      for (i = 0; i < P4EST_CHILDREN; ++i) {
        children[i].p.piggy1.which_tree = rec->which_tree;
      }
      sc_array_init_data (&childrenview, children,
                          sizeof (p4est_quadrant_t), P4EST_CHILDREN);
      conchildren =
        rec->children_fn (rec->p4est, &childrenview, &child_indices);
    }
    P4EST_ASSERT (child_indices.elem_count <= P4EST_CHILDREN);

    /* go into recursion in potentially reordered child order */
    if (conchildren) {
      p4est_quadrant_t    child;

      for (zz = 0; zz < child_indices.elem_count; ++zz) {
        i = (int) *(p4est_topidx_t *) sc_array_index (&child_indices, zz);
        P4EST_ASSERT (0 <= i && i < P4EST_CHILDREN);
        P4EST_ASSERT (split[i] < split[i + 1]);

        sc_array_init_view (&child_quadrants, quadrants,
                            split[i], split[i + 1] - split[i]);
        p4est_quadrant_child (quadrant, &child, i);
        p4est_reorder_recursion (rec, &child, &child_quadrants, chact);
        sc_array_reset (&child_quadrants);
      }
    }
    sc_array_reset (&child_indices);
  }

  /* follow down single quadrant, either when points remain or there never were any */
  if (conchildren && is_leaf && !is_same) {
    int                 level_diff;
    p4est_quadrant_t    child;

    P4EST_ASSERT (!rec->skip);
    P4EST_ASSERT (qcount == 1);
    level_diff = q->level - quadrant->level;
    P4EST_ASSERT (level_diff > 0);

    /* construct the search window one level finer than current quadrant */
    if (level_diff > 1) {
      p4est_quadrant_ancestor (q, quadrant->level + 1, &child);
      P4EST_ASSERT (quadrant->level < child.level && child.level < q->level);
    }
    else {
      P4EST_ASSERT (level_diff == 1);
      child = *q;               /* child is not part of the actual forest storage */
    }
    P4EST_ASSERT (p4est_quadrant_is_parent (quadrant, &child));
    p4est_reorder_recursion (rec, &child, quadrants, chact);
  }

  /* free temporary point storage for children */
  if (chact != NULL) {
    sc_array_reset (chact);
  }

  /* always call post callback on search quadrant, ignoring its return value. */
  if (rec->post_quadrant_fn != NULL) {
    rec->post_quadrant_fn (rec->p4est, rec->which_tree,
                           quadrant, local_num, NULL);
  }
}

void
p4est_search_reorder (p4est_t * p4est, int skip_levels,
                      p4est_search_reorder_t reorder_fn,
                      p4est_search_query_t pre_quadrant_fn,
                      p4est_search_query_t post_quadrant_fn,
                      p4est_search_query_t point_fn, sc_array_t * points)
{
  sc_array_t         *tquadrants;
  sc_array_t         *root_indices;
  p4est_topidx_t      num_local_trees, tt, rit;
  p4est_tree_t       *tree;
  p4est_quadrant_t    root;
  p4est_local_recursion_t srec, *rec = &srec;

  /* check call convention */
  P4EST_ASSERT (p4est != NULL);
  P4EST_ASSERT (points == NULL || point_fn != NULL);

  /* we do nothing if there is nothing we can do */
  if (pre_quadrant_fn == NULL && post_quadrant_fn == NULL && point_fn == NULL) {
    return;
  }

  /* do nothing if there are no local quadrants */
  if (p4est->first_local_tree > p4est->last_local_tree) {
    return;
  }
  num_local_trees = p4est->last_local_tree - p4est->first_local_tree + 1;
  P4EST_ASSERT (num_local_trees > 0);

  /* reorder trees if so desired */
  root_indices = NULL;
  if (reorder_fn != NULL) {
    p4est_quadrant_t   *proot;
    int                 contrees;

    tquadrants =
      sc_array_new_count (sizeof (p4est_quadrant_t), num_local_trees);
    root_indices =
      sc_array_new_count (sizeof (p4est_topidx_t), num_local_trees);
    for (tt = p4est->first_local_tree; tt <= p4est->last_local_tree; ++tt) {
      rit = tt - p4est->first_local_tree;
      proot = (p4est_quadrant_t *) sc_array_index (tquadrants, rit);
      p4est_quadrant_set_morton (proot, 0, 0);
      proot->p.piggy1.which_tree = tt;
      *(p4est_topidx_t *) sc_array_index (root_indices, rit) = rit;
    }

    /* invoke reorder callback */
    contrees = reorder_fn (p4est, tquadrants, root_indices);
    sc_array_destroy (tquadrants);
    if (!contrees) {
      /* returning false breaks recursion */
      sc_array_destroy (root_indices);
      return;
    }
  }

  /* set recursion context */
  rec->p4est = p4est;
  rec->which_tree = -1;
  rec->call_post = 1;
  rec->children_fn = reorder_fn;
  rec->quadrant_fn = NULL;
  rec->pre_quadrant_fn = pre_quadrant_fn;
  rec->post_quadrant_fn = post_quadrant_fn;
  rec->point_fn = point_fn;
  rec->points = points;
  rec->skip = skip_levels;

  /* go through tree recursions in proper order */
  for (tt = p4est->first_local_tree; tt <= p4est->last_local_tree; ++tt) {
    rec->which_tree = root_indices != NULL ? *(p4est_topidx_t *)
      sc_array_index (root_indices, tt - p4est->first_local_tree) +
      p4est->first_local_tree : tt;
    P4EST_ASSERT (p4est->first_local_tree <= rec->which_tree &&
                  rec->which_tree <= p4est->last_local_tree);

    /* grab complete tree quadrant array */
    tree = p4est_tree_array_index (p4est->trees, rec->which_tree);
    tquadrants = &tree->quadrants;

    /* the recursion shrinks the search quadrant whenever possible */
    p4est_quadrant_set_morton (&root, 0, 0);
    p4est_reorder_recursion (rec, &root, tquadrants, NULL);
  }

  /* cleanup sorted tree indices */
  P4EST_ASSERT ((reorder_fn == NULL) == (root_indices == NULL));
  if (reorder_fn != NULL) {
    sc_array_destroy (root_indices);
  }
}

void
p4est_search (p4est_t * p4est, p4est_search_query_t quadrant_fn,
              p4est_search_query_t point_fn, sc_array_t * points)
{
  p4est_search_local (p4est, 0, quadrant_fn, point_fn, points);
}

static              size_t
p4est_traverse_array_index (sc_array_t * array, p4est_topidx_t tt)
{
  P4EST_ASSERT (array != NULL);
  P4EST_ASSERT (array->elem_size == sizeof (size_t));
  P4EST_ASSERT (tt >= 0);

  return *(size_t *) sc_array_index (array, (size_t) tt);
}

static              size_t
p4est_traverse_type_tree (sc_array_t * array, size_t pindex, void *data)
{
  p4est_quadrant_t   *pos;

  P4EST_ASSERT (data == NULL);
  P4EST_ASSERT (array != NULL);
  pos = p4est_quadrant_array_index (array, pindex);
  P4EST_ASSERT (pos->p.which_tree >= 0);

  return (size_t) pos->p.which_tree;
}

/** Check whether a processor begins entering a given quadrant.
 * \param [in] p4est    Used for its partition markers.
 * \param [in] quadrant This quadrant's q->p.which_tree field must match p's.
 * \param [in] p        Processor number that must start in quadrant's tree.
 * \return              True if and only if \b p starts with \p quadrant.
 */
static int
p4est_traverse_is_clean_start (const p4est_quadrant_t *gfp, int num_procs,
                               p4est_topidx_t num_trees,
                               p4est_quadrant_t *quadrant, int p)
{
  const p4est_quadrant_t *marker;

  P4EST_ASSERT (gfp != NULL);
  P4EST_ASSERT (quadrant != NULL);
  P4EST_ASSERT (0 <= p && p <= num_procs);

  marker = gfp + p;
  P4EST_ASSERT (marker->level == P4EST_QMAXLEVEL);
  P4EST_ASSERT (0 <= marker->p.which_tree &&
                marker->p.which_tree <= num_trees);
  P4EST_ASSERT (marker->p.which_tree == quadrant->p.which_tree);

  return marker->x == quadrant->x && marker->y == quadrant->y
#ifdef P4_TO_P8
    && marker->z == quadrant->z
#endif
    ;
}

static              size_t
p4est_traverse_type_childid (sc_array_t * array, size_t pindex, void *data)
{
  p4est_quadrant_t   *quadrant = (p4est_quadrant_t *) data;
  p4est_quadrant_t   *pos;

  P4EST_ASSERT (array != NULL);
  P4EST_ASSERT (data != NULL);
  P4EST_ASSERT (p4est_quadrant_is_valid (quadrant));

  pos = p4est_quadrant_array_index (array, pindex);
  P4EST_ASSERT (pos->p.which_tree >= 0);
  P4EST_ASSERT (p4est_quadrant_is_ancestor (quadrant, pos));

  return (size_t) p4est_quadrant_ancestor_id (pos, quadrant->level + 1);
}

#ifdef P4EST_ENABLE_DEBUG

static int          p4est_traverse_is_valid_quadrant
  (const p4est_quadrant_t *gfp, int num_procs, p4est_topidx_t num_trees,
   p4est_topidx_t which_tree, const p4est_quadrant_t *quadrant,
   int pfirst, int plast)
{
  p4est_quadrant_t    desc;

  P4EST_ASSERT (gfp != NULL);
  P4EST_ASSERT (0 <= which_tree && which_tree < num_trees);
  P4EST_ASSERT (p4est_quadrant_is_valid (quadrant));
  P4EST_ASSERT (0 <= pfirst && pfirst <= plast && plast < num_procs);

  /* check that pfirst is really the first owner in this quadrant */
  p4est_quadrant_first_descendant (quadrant, &desc, P4EST_QMAXLEVEL);
  if (!p4est_comm_is_owner_gfp (gfp, num_procs, num_trees,
                                which_tree, &desc, pfirst)) {
    return 0;
  }

  /* check that plast is really the last owner in this quadrant */
  p4est_quadrant_last_descendant (quadrant, &desc, P4EST_QMAXLEVEL);
  if (!p4est_comm_is_owner_gfp (gfp, num_procs, num_trees,
                                which_tree, &desc, plast)) {
    return 0;
  }

#if 0
  /* this is redundant after the checks above */
  /* we're optimizing this away to avoid passing gfq in addition */
  P4EST_ASSERT (!p4est_comm_is_empty_gfq (gfq, num_procs, pfirst) &&
                !p4est_comm_is_empty_gfq (gfq, num_procs, plast));
#endif

  return 1;
}

static int          p4est_traverse_is_valid_tree
  (const p4est_quadrant_t *gfp, int num_procs, p4est_topidx_t num_trees,
   p4est_topidx_t which_tree, int pfirst, int plast)
{
  p4est_quadrant_t    root;

  P4EST_ASSERT (gfp != NULL);
  P4EST_ASSERT (0 <= which_tree && which_tree < num_trees);
  P4EST_ASSERT (0 <= pfirst && pfirst <= plast && plast < num_procs);

  p4est_quadrant_set_morton (&root, 0, 0);

  return p4est_traverse_is_valid_quadrant
    (gfp, num_procs, num_trees, which_tree, &root, pfirst, plast);
}

#endif /* P4EST_ENABLE_DEBUG */

/** This recursion context saves on the number of parameters passed. */
typedef struct p4est_partition_recursion
{
  p4est_t            *user_p4est;       /**< Forest passed to callbacks. */
  const p4est_gloidx_t *gfq;            /**< Global first quadrant array. */
  const p4est_quadrant_t *gfp;          /**< Global first positiion array. */
  int                 num_procs;        /**< Number of ranks in partition. */
  p4est_topidx_t      num_trees;        /**< Number of trees in partition. */
  p4est_topidx_t      which_tree;       /**< Current tree number. */
  int                 call_post;        /**< Boolean to call quadrant twice. */
  p4est_search_partition_t quadrant_fn; /**< Per-quadrant callback. */
  p4est_search_partition_t point_fn;    /**< Per-point callback. */
  sc_array_t         *points;           /**< Array of points to search. */
  sc_array_t         *position_array;   /**< Array view of partition's
                                             global_first_position */
}
p4est_partition_recursion_t;

static void
p4est_partition_recursion (const p4est_partition_recursion_t * rec,
                           p4est_quadrant_t * quadrant, int pfirst, int plast,
                           sc_array_t * actives)
{
  int                 i;
  int                 is_match;
  int                 cpfirst, cplast, cpnext;
  size_t              zz, *pz, *qz;
  size_t              act_count;
  p4est_quadrant_t    child;
  sc_array_t          pview, offsets;
  sc_array_t          child_actives, *chact;

  P4EST_ASSERT (rec != NULL);
  P4EST_ASSERT (0 <= rec->which_tree && rec->which_tree < rec->num_trees);
  P4EST_ASSERT (p4est_traverse_is_valid_quadrant
                (rec->gfp, rec->num_procs, rec->num_trees,
                 rec->which_tree, quadrant, pfirst, plast));
  P4EST_ASSERT (quadrant != NULL && p4est_quadrant_is_valid (quadrant));
  P4EST_ASSERT (quadrant->p.which_tree == rec->which_tree);

  /* As an optimization we pass a NULL actives array to every root. */
  if (rec->points != NULL && actives == NULL) {
    act_count = rec->points->elem_count;
  }
  else {
    P4EST_ASSERT ((rec->points == NULL) == (actives == NULL));
    P4EST_ASSERT (rec->points == NULL ||
                  actives->elem_count <= rec->points->elem_count);
    act_count = actives == NULL ? 0 : actives->elem_count;
  }

  /* return if there are no active points */
  if (rec->points != NULL && act_count == 0)
    return;

  /* execute pre-quadrant callback if present, which may stop the recursion */
  if (rec->quadrant_fn != NULL &&
      !rec->quadrant_fn (rec->user_p4est, rec->which_tree,
                         quadrant, pfirst, plast, NULL)) {
    return;
  }

  /* check out points */
  if (rec->points == NULL) {
    /* we have called the callback already.  Maybe we are done */
    if (pfirst == plast) {
      return;
    }
    chact = NULL;
  }
  else {
    /* query callback for all points and return if none remain */
    chact = &child_actives;
    sc_array_init (chact, sizeof (size_t));
    for (zz = 0; zz < act_count; ++zz) {
      pz = actives == NULL ? &zz : (size_t *) sc_array_index (actives, zz);
      is_match = rec->point_fn (rec->user_p4est, rec->which_tree,
                                quadrant, pfirst, plast,
                                sc_array_index (rec->points, *pz));
      if (!(pfirst == plast) && is_match) {
        qz = (size_t *) sc_array_push (chact);
        *qz = *pz;
      }
    }

    /* call post-quadrant callback, which may also terminate the recursion */
    if (rec->call_post && rec->quadrant_fn != NULL &&
        !rec->quadrant_fn (rec->user_p4est, rec->which_tree,
                           quadrant, pfirst, plast, NULL)) {
      /* clears memory and will trigger the return below */
      sc_array_reset (chact);
    }

    if (chact->elem_count == 0) {
      /* with zero members there is no need to call sc_array_reset */
      return;
    }
  }

  /* the one-processor branch has returned above */
  P4EST_ASSERT (!(pfirst == plast));
  P4EST_ASSERT (quadrant->level < P4EST_QMAXLEVEL);

  /* find the processors for all children of the quadrant */
  sc_array_init_view (&pview, rec->position_array,
                      pfirst + 1, plast - pfirst);
  sc_array_init_size (&offsets, sizeof (size_t), P4EST_CHILDREN + 1);
  sc_array_split (&pview, &offsets, P4EST_CHILDREN,
                  p4est_traverse_type_childid, quadrant);
  P4EST_ASSERT (offsets.elem_count == (size_t) (P4EST_CHILDREN + 1));
  P4EST_ASSERT (p4est_traverse_array_index
                (&offsets, P4EST_CHILDREN) == (size_t) (plast - pfirst));
  P4EST_ASSERT (p4est_traverse_array_index (&offsets, 0) == 0);

  /* go through the quadrant's children */
  child.p.which_tree = rec->which_tree;
  for (cpfirst = pfirst + 1, i = 0; i < P4EST_CHILDREN; cpfirst = cpnext, ++i) {
    p4est_quadrant_child (quadrant, &child, i);

    /* determine the exclusive upper bound of processors starting in child */
    cpnext = p4est_traverse_array_index (&offsets, i + 1) + pfirst + 1;
    P4EST_ASSERT (cpfirst <= cpnext && cpnext <= plast + 1);

    /* fix the last processor in child, which is known at this point */
    P4EST_ASSERT (cpnext > 0);
    cplast = cpnext - 1;

    /* now check multiple cases for the beginning processor */
    if (cpfirst < cpnext) {
      /* at least one processor starts in this child */

      if (p4est_traverse_is_clean_start
          (rec->gfp, rec->num_procs, rec->num_trees, &child, cpfirst)) {
        /* cpfirst starts at the tree's first descendant but may be empty */
        P4EST_ASSERT (i > 0);
        while (P4EST_COMM_IS_EMPTY_GFQ_GFP
               (rec->gfq, rec->gfp, rec->num_procs, cpfirst)) {
          ++cpfirst;
          P4EST_ASSERT (p4est_traverse_type_childid
                        (rec->position_array, cpfirst, quadrant) ==
                        (size_t) i);
        }
      }
      else {
        /* there must be exactly one processor before us in this child */
        --cpfirst;
        P4EST_ASSERT (cpfirst == pfirst ||
                      p4est_traverse_type_childid
                      (rec->position_array, cpfirst, quadrant) < (size_t) i);
      }
    }
    else {
      /* this whole child is owned by one processor */
      cpfirst = cplast;
    }

    /* we should have found tight bounds on processors for this child */
    P4EST_ASSERT (i > 0 || pfirst == cpfirst);
    P4EST_ASSERT (i < P4EST_CHILDREN - 1 || plast == cplast);
    P4EST_ASSERT (pfirst <= cpfirst && cpfirst <= cplast && cplast <= plast);

    /* we know these are non-negative; check before casting to unsigned */
    P4EST_ASSERT (cplast >= 0 && cpnext >= 0 && plast + 1 >= 0);

    /* These casts remove compiler warnings due to the assumption of the
     * compiler under -O3 that there cannot happen a signed overflow.
     */
    P4EST_ASSERT ((unsigned) cplast <= (unsigned) cpnext
                  && (unsigned) cpnext <= (unsigned) plast + 1);
    P4EST_ASSERT (cplast == pfirst ||
                  p4est_traverse_type_childid
                  (rec->position_array, cplast, quadrant) <= (size_t) i);
    P4EST_ASSERT (p4est_traverse_is_valid_quadrant
                  (rec->gfp, rec->num_procs, rec->num_trees,
                   rec->which_tree, &child, cpfirst, cplast));

    /* go deeper into the recursion */
    p4est_partition_recursion (rec, &child, cpfirst, cplast, chact);
  }

  /* this is it */
  if (chact != NULL) {
    sc_array_reset (chact);
  }
  sc_array_reset (&offsets);
  sc_array_reset (&pview);
}

static void         p4est_search_partition_internal
  (const p4est_gloidx_t *gfq, const p4est_quadrant_t *gfp,
   int nmemb, p4est_topidx_t num_trees, int call_post, p4est_t *user_p4est,
   p4est_search_partition_t quadrant_fn, p4est_search_partition_t point_fn,
   sc_array_t *points);

void
p4est_search_partition (p4est_t *p4est, int call_post,
                        p4est_search_partition_t quadrant_fn,
                        p4est_search_partition_t point_fn,
                        sc_array_t *points)
{
  P4EST_ASSERT (p4est != NULL);
  P4EST_ASSERT (p4est->connectivity != NULL);

  p4est_search_partition_internal
    (p4est->global_first_quadrant, p4est->global_first_position,
     p4est->mpisize, p4est->connectivity->num_trees,
     call_post, p4est, quadrant_fn, point_fn, points);
}

void
p4est_search_partition_gfx (const p4est_gloidx_t *gfq,
                            const p4est_quadrant_t *gfp,
                            int nmemb, p4est_topidx_t num_trees,
                            int call_post, void *user,
                            p4est_search_partition_t quadrant_fn,
                            p4est_search_partition_t point_fn,
                            sc_array_t *points)
{
  p4est_t             p, *user_p4est = &p;

  /* sanity checks on global first quadrant */
  P4EST_ASSERT (gfq != NULL);
  P4EST_ASSERT (gfq[0] == 0);
  P4EST_ASSERT (gfq[nmemb] >= (p4est_gloidx_t) num_trees);

  /* sanity checks on global first position */
  P4EST_ASSERT (gfp != NULL);
  P4EST_ASSERT (gfp[0].p.which_tree == 0);
  P4EST_ASSERT (gfp[nmemb].x == 0);
  P4EST_ASSERT (gfp[nmemb].y == 0);
#ifdef P4_TO_P8
  P4EST_ASSERT (gfp[nmemb].z == 0);
#endif
  P4EST_ASSERT (gfp[nmemb].p.which_tree == num_trees);

  /* conjure up call convention for partition search */
  memset (user_p4est, 0, sizeof (p4est_t));
  user_p4est->user_pointer = user;
  p4est_search_partition_internal
    (gfq, gfp, nmemb, num_trees, call_post, user_p4est,
     quadrant_fn, point_fn, points);
}

void
p4est_search_partition_gfp (const p4est_quadrant_t *gfp, int nmemb,
                            p4est_topidx_t num_trees, int call_post,
                            void *user, p4est_search_partition_t quadrant_fn,
                            p4est_search_partition_t point_fn,
                            sc_array_t *points)
{
  p4est_t             p, *user_p4est = &p;

  /* sanity checks on global first position */
  P4EST_ASSERT (gfp != NULL);
  P4EST_ASSERT (gfp[0].p.which_tree == 0);
  P4EST_ASSERT (gfp[nmemb].x == 0);
  P4EST_ASSERT (gfp[nmemb].y == 0);
#ifdef P4_TO_P8
  P4EST_ASSERT (gfp[nmemb].z == 0);
#endif
  P4EST_ASSERT (gfp[nmemb].p.which_tree == num_trees);

  /* conjure up call convention for partition search */
  memset (user_p4est, 0, sizeof (p4est_t));
  user_p4est->user_pointer = user;
  p4est_search_partition_internal
    (NULL, gfp, nmemb, num_trees, call_post, user_p4est,
     quadrant_fn, point_fn, points);
}

void                p4est_search_partition_internal
  (const p4est_gloidx_t *gfq, const p4est_quadrant_t *gfp,
   int nmemb, p4est_topidx_t num_trees, int call_post, p4est_t *user_p4est,
   p4est_search_partition_t quadrant_fn, p4est_search_partition_t point_fn,
   sc_array_t *points)
{
  const int           num_procs = nmemb;
  int                 pfirst, plast, pnext;
  sc_array_t          position_array;
  sc_array_t         *tree_offsets;
  p4est_topidx_t      tt;
  p4est_quadrant_t    root;
  p4est_partition_recursion_t srec, *rec = &srec;

  /* we do nothing if there is nothing to be done */
  P4EST_ASSERT (gfp != NULL);
  P4EST_ASSERT (points == NULL || point_fn != NULL);
  if (quadrant_fn == NULL && points == NULL) {
    return;
  }

  /* array to split is the p4est partition marker.  A const crime */
  /* it is important to include the highest tree number plus one */
  sc_array_init_data (&position_array, (p4est_quadrant_t *) gfp,
                      sizeof (p4est_quadrant_t), num_procs + 1);

  /* the enumerable type is the tree number -- we know the size already */
  tree_offsets = sc_array_new_size (sizeof (size_t), num_trees + 2);

  /* split processors into tree-wise sections, going one beyond */
  sc_array_split (&position_array, tree_offsets, num_trees + 1,
                  p4est_traverse_type_tree, NULL);
  P4EST_ASSERT (tree_offsets->elem_count == (size_t) (num_trees + 2));
  P4EST_ASSERT (p4est_traverse_array_index
                (tree_offsets, num_trees + 1) == (size_t) num_procs + 1);
  P4EST_ASSERT (p4est_traverse_array_index
                (tree_offsets, num_trees) <= (size_t) num_procs);
  P4EST_ASSERT (p4est_traverse_array_index (tree_offsets, 0) == 0);

  /* now loop through all trees, local or not */
  rec->user_p4est = user_p4est;
  rec->gfq = gfq;
  rec->gfp = gfp;
  rec->num_procs = num_procs;
  rec->num_trees = num_trees;
  rec->which_tree = -1;
  rec->call_post = call_post;
  rec->quadrant_fn = quadrant_fn;
  rec->point_fn = point_fn;
  rec->points = points;
  rec->position_array = &position_array;
  p4est_quadrant_set_morton (&root, 0, 0);
  for (pfirst = 0, tt = 0; tt < num_trees; pfirst = pnext, ++tt) {
    /* pfirst is the first processor indexed for this tree */
    rec->which_tree = root.p.which_tree = tt;

    /* determine the exclusive upper bound of processors starting in this tree */
    pnext = p4est_traverse_array_index (tree_offsets, tt + 1);
    P4EST_ASSERT (pfirst <= pnext && pnext <= num_procs);

    /* fix the last processor in the tree, which is known at this point */
    P4EST_ASSERT (pnext > 0);
    plast = pnext - 1;

    /* now check multiple cases for the beginning processor */
    if (pfirst < pnext) {
      /* at least one processor starts in this tree */

      if (p4est_traverse_is_clean_start
          (rec->gfp, rec->num_procs, rec->num_trees, &root, pfirst)) {
        /* pfirst starts at the tree's first descendant but may be empty */
        while (P4EST_COMM_IS_EMPTY_GFQ_GFP
               (rec->gfq, rec->gfp, rec->num_procs, pfirst)) {
          ++pfirst;
          P4EST_ASSERT (p4est_traverse_type_tree
                        (&position_array, pfirst, NULL) == (size_t) tt);
        }
      }
      else {
        /* there must be exactly one processor before us in this tree */
        --pfirst;
        P4EST_ASSERT (p4est_traverse_type_tree
                      (&position_array, pfirst, NULL) < (size_t) tt);
      }
    }
    else {
      /* this whole tree is owned by one processor */
      pfirst = plast;
    }

    /* we should have found tight bounds on processors for this tree */
    P4EST_ASSERT (pfirst <= plast && plast < num_procs);

    /* we know these are non-negative; check before casting to unsigned */
    P4EST_ASSERT (plast >= 0 && pnext >= 0 && num_procs >= 0);

    /* These casts remove compiler warnings due to the assumption of the
     * compiler under -O3 that there cannot happen a signed overflow.
     */
    P4EST_ASSERT ((unsigned) plast <= (unsigned) pnext
                  && (unsigned) pnext <= (unsigned) num_procs);
    P4EST_ASSERT (p4est_traverse_type_tree
                  (&position_array, plast, NULL) <= (size_t) tt);
    P4EST_ASSERT (p4est_traverse_is_valid_tree
                  (rec->gfp, rec->num_procs, rec->num_trees,
                   tt, pfirst, plast));

    /* go into recursion for this tree */
    p4est_partition_recursion (rec, &root, pfirst, plast, NULL);
  }

  /* cleanup */
  sc_array_destroy (tree_offsets);
  sc_array_reset (&position_array);
}

/** This recursion context saves on the number of parameters passed. */
typedef struct p4est_all_recursion
{
  p4est_t            *p4est;            /**< Forest being traversed. */
  const p4est_gloidx_t *gfq;            /**< Global first quadrant array. */
  const p4est_quadrant_t *gfp;          /**< Global first positiion array. */
  int                 num_procs;        /**< Number of ranks in partition. */
  p4est_topidx_t      num_trees;        /**< Number of trees in partition. */
  p4est_topidx_t      which_tree;       /**< Current tree number. */
  int                 call_post;        /**< Boolean to call quadrant twice. */
  p4est_search_all_t  quadrant_fn;      /**< Per-quadrant callback. */
  p4est_search_all_t  point_fn;         /**< Per-point callback. */
  sc_array_t         *points;           /**< Array of points to search. */
  sc_array_t         *position_array;   /**< Array view of p4est's
                                             global_first_position */
}
p4est_all_recursion_t;

static void
p4est_all_recursion (const p4est_all_recursion_t * rec,
                     p4est_quadrant_t * quadrant, int pfirst, int plast,
                     sc_array_t * quadrants, sc_array_t * actives)
{
  int                 i;
  int                 proceed;
  int                 is_leaf, is_match;
  int                 cpfirst, cplast, cpnext;
  size_t              qcount, act_count;
  size_t              zz, *pz, *qz;
  size_t              split[P4EST_CHILDREN + 1];
  p4est_locidx_t      local_num;
  p4est_quadrant_t   *q, child;
  sc_array_t          pview, offsets;
  sc_array_t          child_quadrants, *chpass, child_actives, *chact;

  P4EST_ASSERT (rec != NULL);
  P4EST_ASSERT (0 <= rec->which_tree &&
                rec->which_tree < rec->p4est->connectivity->num_trees);
  P4EST_ASSERT (p4est_traverse_is_valid_quadrant
                (rec->gfp, rec->num_procs, rec->num_trees,
                 rec->which_tree, quadrant, pfirst, plast));
  P4EST_ASSERT (quadrant != NULL && p4est_quadrant_is_valid (quadrant));
  P4EST_ASSERT (quadrant->p.which_tree == rec->which_tree);

  /* As an optimization we pass a NULL actives array to every root. */
  if (rec->points != NULL && actives == NULL) {
    act_count = rec->points->elem_count;
  }
  else {
    P4EST_ASSERT ((rec->points == NULL) == (actives == NULL));
    P4EST_ASSERT (rec->points == NULL ||
                  actives->elem_count <= rec->points->elem_count);
    act_count = actives == NULL ? 0 : actives->elem_count;
  }

  /* return if there are no active points */
  if (rec->points != NULL && act_count == 0)
    return;

  /* check out local quadrant portion */
  qcount = 0;
  is_leaf = 0;
  local_num = -1;
  if (quadrants != NULL && (qcount = quadrants->elem_count) > 0) {
    /* determine leaf situation */
    q = p4est_quadrant_array_index (quadrants, 0);
    if (!p4est_quadrant_is_equal (quadrant, q)) {
      P4EST_ASSERT (p4est_quadrant_is_ancestor (quadrant, q));

      /* since the parallel partition needs to be followed, we cannot
       * optimize for the case that the quadrants array may be contained
       * in a descendent of quadrant */
    }
    else {
      p4est_locidx_t      offset;
      p4est_tree_t       *tree;

      /* we have reached a leaf of the forest */
      P4EST_ASSERT (qcount == 1);
      P4EST_ASSERT (pfirst == plast);
      P4EST_ASSERT (pfirst == rec->p4est->mpirank);
      is_leaf = 1;

      /* determine offset of quadrant in local forest */
      tree = p4est_tree_array_index (rec->p4est->trees, rec->which_tree);
      offset = (p4est_locidx_t) ((quadrants->array - tree->quadrants.array)
                                 / sizeof (p4est_quadrant_t));
      P4EST_ASSERT (offset >= 0 &&
                    (size_t) offset < tree->quadrants.elem_count);
      local_num = tree->quadrants_offset + offset;

      /* use the actual quadrant stored in the forest for the callbacks */
      quadrant = q;
    }
  }
  P4EST_ASSERT ((!is_leaf) == (local_num == -1));

  /* execute pre-quadrant callback if present, which may stop the recursion */
  if (rec->quadrant_fn != NULL &&
      !rec->quadrant_fn (rec->p4est, rec->which_tree,
                         quadrant, pfirst, plast, local_num, NULL)) {
    return;
  }

  /*
   * We continue the recursion if and only if any of these conditions holds:
   * pfirst < plast
   * pfirst == plast == mpirank and not a leaf
   */
  proceed = pfirst < plast || (pfirst == rec->p4est->mpirank && !is_leaf);

  /* check out the points */
  if (rec->points == NULL) {
    /* we have called the callback already.  Maybe we are done */
    if (!proceed) {
      /* if this is a leaf we necessarily enter here */
      return;
    }
    chact = NULL;
  }
  else {
    /* query callback for all points and return if none remain */
    chact = &child_actives;
    sc_array_init (chact, sizeof (size_t));
    for (zz = 0; zz < act_count; ++zz) {
      pz = actives == NULL ? &zz : (size_t *) sc_array_index (actives, zz);
      is_match = rec->point_fn (rec->p4est, rec->which_tree,
                                quadrant, pfirst, plast, local_num,
                                sc_array_index (rec->points, *pz));
      if (proceed && is_match) {
        qz = (size_t *) sc_array_push (chact);
        *qz = *pz;
      }
    }

    /* call post-quadrant callback, which may also terminate the recursion */
    if (rec->call_post && rec->quadrant_fn != NULL &&
        !rec->quadrant_fn (rec->p4est, rec->which_tree,
                           quadrant, pfirst, plast, local_num, NULL)) {
      /* clears memory and will trigger the return below */
      sc_array_reset (chact);
    }

    if (chact->elem_count == 0) {
      /* with zero members there is no need to call sc_array_reset */
      return;
    }
  }

  /* the stopping condition (including we being a leaf) has returned above */
  P4EST_ASSERT (proceed);
  P4EST_ASSERT (quadrant->level < P4EST_QMAXLEVEL);

  /* find the processors for all children of the quadrant */
  sc_array_init_view (&pview, rec->position_array,
                      pfirst + 1, plast - pfirst);
  sc_array_init_size (&offsets, sizeof (size_t), P4EST_CHILDREN + 1);
  sc_array_split (&pview, &offsets, P4EST_CHILDREN,
                  p4est_traverse_type_childid, quadrant);
  P4EST_ASSERT (offsets.elem_count == (size_t) (P4EST_CHILDREN + 1));
  P4EST_ASSERT (p4est_traverse_array_index
                (&offsets, P4EST_CHILDREN) == (size_t) (plast - pfirst));
  P4EST_ASSERT (p4est_traverse_array_index (&offsets, 0) == 0);

  /* split quadrant array for local portion */
  if (quadrants != NULL) {
    p4est_split_array (quadrants, (int) quadrant->level, split);
  }

  /* go through the quadrant's children */
  child.p.which_tree = rec->which_tree;
  for (cpfirst = pfirst + 1, i = 0; i < P4EST_CHILDREN; cpfirst = cpnext, ++i) {
    p4est_quadrant_child (quadrant, &child, i);

    /* determine the exclusive upper bound of processors starting in child */
    cpnext = p4est_traverse_array_index (&offsets, i + 1) + pfirst + 1;
    P4EST_ASSERT (cpfirst <= cpnext && cpnext <= plast + 1);

    /* fix the last processor in child, which is known at this point */
    P4EST_ASSERT (cpnext > 0);
    cplast = cpnext - 1;

    /* now check multiple cases for the beginning processor */
    if (cpfirst < cpnext) {
      /* at least one processor starts in this child */

      if (p4est_traverse_is_clean_start
          (rec->gfp, rec->num_procs, rec->num_trees, &child, cpfirst)) {
        /* cpfirst starts at the tree's first descendant but may be empty */
        P4EST_ASSERT (i > 0);
        while (p4est_comm_is_empty (rec->p4est, cpfirst)) {
          ++cpfirst;
          P4EST_ASSERT (p4est_traverse_type_childid
                        (rec->position_array, cpfirst, quadrant) ==
                        (size_t) i);
        }
      }
      else {
        /* there must be exactly one processor before us in this child */
        --cpfirst;
        P4EST_ASSERT (cpfirst == pfirst ||
                      p4est_traverse_type_childid
                      (rec->position_array, cpfirst, quadrant) < (size_t) i);
      }
    }
    else {
      /* this whole child is owned by one processor */
      cpfirst = cplast;
    }

    /* we should have found tight bounds on processors for this child */
    P4EST_ASSERT (i > 0 || pfirst == cpfirst);
    P4EST_ASSERT (i < P4EST_CHILDREN - 1 || plast == cplast);
    P4EST_ASSERT (pfirst <= cpfirst && cpfirst <= cplast && cplast <= plast);

    /* we know these are non-negative; check before casting to unsigned */
    P4EST_ASSERT (cplast >= 0 && cpnext >= 0 && plast + 1 >= 0);

    /* These casts remove compiler warnings due to the assumption of the
     * compiler under -O3 that there cannot happen a signed overflow.
     */
    P4EST_ASSERT ((unsigned) cplast <= (unsigned) cpnext
                  && (unsigned) cpnext <= (unsigned) plast + 1);
    P4EST_ASSERT (cplast == pfirst ||
                  p4est_traverse_type_childid
                  (rec->position_array, cplast, quadrant) <= (size_t) i);
    P4EST_ASSERT (p4est_traverse_is_valid_quadrant
                  (rec->gfp, rec->num_procs, rec->num_trees,
                   rec->which_tree, &child, cpfirst, cplast));

    /* designate the subarray of local quadrants */
    chpass = NULL;
    if (quadrants != NULL && split[i] < split[i + 1]) {
      chpass = &child_quadrants;
      sc_array_init_view (chpass, quadrants,
                          split[i], split[i + 1] - split[i]);
    }

    /* go deeper into the recursion */
    p4est_all_recursion (rec, &child, cpfirst, cplast, chpass, chact);
    if (chpass != NULL) {
      sc_array_reset (&child_quadrants);
    }
  }

  /* this is it */
  if (chact != NULL) {
    sc_array_reset (chact);
  }
  sc_array_reset (&offsets);
  sc_array_reset (&pview);
}

void
p4est_search_all (p4est_t * p4est,
                  int call_post, p4est_search_all_t quadrant_fn,
                  p4est_search_all_t point_fn, sc_array_t * points)
{
  const int           num_procs = p4est->mpisize;
  const p4est_topidx_t num_trees = p4est->connectivity->num_trees;
  int                 pfirst, plast, pnext;
  sc_array_t          position_array;
  sc_array_t         *tree_offsets;
  sc_array_t         *tquadrants;
  p4est_topidx_t      tt;
  p4est_tree_t       *tree;
  p4est_quadrant_t    root;
  p4est_all_recursion_t srec, *rec = &srec;

  /* we do nothing if there is nothing to be done */
  P4EST_ASSERT (p4est != NULL);
  P4EST_ASSERT (points == NULL || point_fn != NULL);
  if (quadrant_fn == NULL && points == NULL) {
    return;
  }

  /* array to split is the p4est partition marker */
  /* it is important to include the highest tree number plus one */
  sc_array_init_data (&position_array, p4est->global_first_position,
                      sizeof (p4est_quadrant_t), num_procs + 1);

  /* the enumerable type is the tree number -- we know the size already */
  tree_offsets = sc_array_new_size (sizeof (size_t), num_trees + 2);

  /* split processors into tree-wise sections, going one beyond */
  sc_array_split (&position_array, tree_offsets, num_trees + 1,
                  p4est_traverse_type_tree, NULL);
  P4EST_ASSERT (tree_offsets->elem_count == (size_t) (num_trees + 2));
  P4EST_ASSERT (p4est_traverse_array_index
                (tree_offsets, num_trees + 1) == (size_t) num_procs + 1);
  P4EST_ASSERT (p4est_traverse_array_index
                (tree_offsets, num_trees) <= (size_t) num_procs);
  P4EST_ASSERT (p4est_traverse_array_index (tree_offsets, 0) == 0);

  /* now loop through all trees, local or not */
  rec->p4est = p4est;
  rec->gfq = p4est->global_first_quadrant;
  rec->gfp = p4est->global_first_position;
  rec->num_procs = num_procs;
  rec->num_trees = num_trees;
  rec->which_tree = -1;
  rec->call_post = call_post;
  rec->quadrant_fn = quadrant_fn;
  rec->point_fn = point_fn;
  rec->points = points;
  rec->position_array = &position_array;
  p4est_quadrant_set_morton (&root, 0, 0);
  for (pfirst = 0, tt = 0; tt < num_trees; pfirst = pnext, ++tt) {
    /* pfirst is the first processor indexed for this tree */
    rec->which_tree = root.p.which_tree = tt;

    /* determine the exclusive upper bound of processors starting in this tree */
    pnext = p4est_traverse_array_index (tree_offsets, tt + 1);
    P4EST_ASSERT (pfirst <= pnext && pnext <= num_procs);

    /* fix the last processor in the tree, which is known at this point */
    P4EST_ASSERT (pnext > 0);
    plast = pnext - 1;

    /* now check multiple cases for the beginning processor */
    if (pfirst < pnext) {
      /* at least one processor starts in this tree */

      if (p4est_traverse_is_clean_start
          (rec->gfp, rec->num_procs, rec->num_trees, &root, pfirst)) {
        /* pfirst starts at the tree's first descendant but may be empty */
        while (p4est_comm_is_empty (p4est, pfirst)) {
          ++pfirst;
          P4EST_ASSERT (p4est_traverse_type_tree
                        (&position_array, pfirst, NULL) == (size_t) tt);
        }
      }
      else {
        /* there must be exactly one processor before us in this tree */
        --pfirst;
        P4EST_ASSERT (p4est_traverse_type_tree
                      (&position_array, pfirst, NULL) < (size_t) tt);
      }
    }
    else {
      /* this whole tree is owned by one processor */
      pfirst = plast;
    }

    /* we should have found tight bounds on processors for this tree */
    P4EST_ASSERT (pfirst <= plast && plast < num_procs);

    /* we know these are non-negative; check before casting to unsigned */
    P4EST_ASSERT (plast >= 0 && pnext >= 0 && num_procs >= 0);

    /* These casts remove compiler warnings due to the assumption of the
     * compiler under -O3 that there cannot happen a signed overflow.
     */
    P4EST_ASSERT ((unsigned) plast <= (unsigned) pnext
                  && (unsigned) pnext <= (unsigned) num_procs);
    P4EST_ASSERT (p4est_traverse_type_tree
                  (&position_array, plast, NULL) <= (size_t) tt);
    P4EST_ASSERT (p4est_traverse_is_valid_tree
                  (rec->gfp, rec->num_procs, rec->num_trees,
                   tt, pfirst, plast));

    /* if this tree is at least partially local, get the local quadrants */
    if (p4est->first_local_tree <= tt && tt <= p4est->last_local_tree) {
      P4EST_ASSERT (pfirst <= p4est->mpirank && p4est->mpirank <= plast);

      /* grab complete tree quadrant array */
      tree = p4est_tree_array_index (p4est->trees, tt);
      tquadrants = &tree->quadrants;

      /* we must not shrink the root quadrant in p4est_search_all */
    }
    else {
      /* this processor is empty or entirely before or after this tree */
      tquadrants = NULL;
    }

    /* go into recursion for this tree */
    p4est_all_recursion (rec, &root, pfirst, plast, tquadrants, NULL);
  }

  /* cleanup */
  sc_array_destroy (tree_offsets);
  sc_array_reset (&position_array);
}
