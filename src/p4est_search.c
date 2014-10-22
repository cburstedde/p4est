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

#ifndef P4_TO_P8
#include <p4est_bits.h>
#include <p4est_search.h>
#else
#include <p8est_bits.h>
#include <p8est_search.h>
#endif

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

static void
p4est_search_recursion (p4est_t * p4est, p4est_topidx_t which_tree,
                        p4est_quadrant_t * quadrant,
                        p4est_search_query_t search_quadrant_fn,
                        p4est_search_query_t search_point_fn,
                        sc_array_t * quadrants,
                        sc_array_t * points, sc_array_t * actives)
{
  int                 i;
  int                 is_leaf, is_match;
  size_t              qcount = quadrants->elem_count;
  size_t              zz, *pz, *qz;
  size_t              split[P4EST_CHILDREN + 1];
  p4est_locidx_t      local_num;
  p4est_quadrant_t   *q, *lq, children[P4EST_CHILDREN];
  sc_array_t          child_quadrants, child_actives;

  /*
   * Invariants of the recursion:
   * 1. quadrant is larger or equal in size than those in the array.
   * 2. quadrant is equal to or an ancestor of those in the array.
   */

  P4EST_ASSERT (actives->elem_count <= points->elem_count);

  /* return if there are no quadrants or active points */
  if (qcount == 0 || actives->elem_count == 0)
    return;

  /* determine leaf situation */
  q = p4est_quadrant_array_index (quadrants, 0);
  if (qcount > 1) {
    P4EST_ASSERT (p4est_quadrant_is_ancestor (quadrant, q));
    is_leaf = 0;
    local_num = -1;
    lq = p4est_quadrant_array_index (quadrants, quadrants->elem_count - 1);
    if (p4est_quadrant_ancestor_id (q, quadrant->level + 1) ==
        p4est_quadrant_ancestor_id (lq, quadrant->level + 1)) {
      p4est_nearest_common_ancestor (q, lq, quadrant);
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
    tree = p4est_tree_array_index (p4est->trees, which_tree);
    offset = (p4est_locidx_t) ((quadrants->array - tree->quadrants.array)
                               / sizeof (p4est_quadrant_t));
    P4EST_ASSERT (offset >= 0 &&
                  (size_t) offset < tree->quadrants.elem_count);
    local_num = tree->quadrants_offset + offset;
    quadrant = q;
  }

  /* execute quadrant callback if present, which may stop the recursion */
  if (search_quadrant_fn != NULL &&
      !search_quadrant_fn (p4est, which_tree, quadrant, local_num, NULL)) {
    return;
  }

  /* query callback for all points and return if none remain */
  sc_array_init (&child_actives, sizeof (size_t));
  for (zz = 0; zz < actives->elem_count; ++zz) {
    pz = (size_t *) sc_array_index (actives, zz);
    is_match = search_point_fn (p4est, which_tree, quadrant, local_num,
                                sc_array_index (points, *pz));
    if (!is_leaf && is_match) {
      qz = (size_t *) sc_array_push (&child_actives);
      *qz = *pz;
    }
  }
  if (child_actives.elem_count == 0)
    return;

  /* leaf situation has returned above */
  P4EST_ASSERT (!is_leaf);

  /* split quadrant array and run recursion */
  p4est_split_array (quadrants, (int) quadrant->level, split);
  p4est_quadrant_childrenv (quadrant, children);
  for (i = 0; i < P4EST_CHILDREN; ++i) {
    if (split[i] < split[i + 1]) {
      sc_array_init_view (&child_quadrants, quadrants,
                          split[i], split[i + 1] - split[i]);
      p4est_search_recursion (p4est, which_tree, &children[i],
                              search_quadrant_fn, search_point_fn,
                              &child_quadrants, points, &child_actives);
      sc_array_reset (&child_quadrants);
    }
  }
  sc_array_reset (&child_actives);
}

void
p4est_search (p4est_t * p4est, p4est_search_query_t search_quadrant_fn,
              p4est_search_query_t search_point_fn, sc_array_t * points)
{
  p4est_topidx_t      jt;
  p4est_tree_t       *tree;
  p4est_quadrant_t    root;
  p4est_quadrant_t    temp, temp2;
  p4est_quadrant_t   *f, *l;
  uint64_t            midx;
  sc_array_t          actives;
  sc_array_t         *tquadrants;
  size_t              zz, *pz;

  for (jt = p4est->first_local_tree; jt <= p4est->last_local_tree; ++jt) {

    /* start recursion with root quadrant */
    p4est_quadrant_set_morton (&root, 0, 0);
    f = NULL;
    l = NULL;
    if (jt == p4est->first_local_tree) {
      f = &p4est->global_first_position[p4est->mpirank];
      p4est_quadrant_last_descendant (&root, &temp, P4EST_QMAXLEVEL);
      l = &temp;
    }
    if (jt == p4est->last_local_tree) {
      if (!f) {
        p4est_quadrant_first_descendant (&root, &temp, P4EST_QMAXLEVEL);
        f = &temp;
      }
      l = &p4est->global_first_position[p4est->mpirank + 1];
      if (l->p.which_tree == jt) {
        midx = p4est_quadrant_linear_id (l, P4EST_QMAXLEVEL);
        p4est_quadrant_set_morton (&temp2, P4EST_QMAXLEVEL, midx - 1);
        l = &temp2;
      }
      else {
        p4est_quadrant_last_descendant (&root, &temp2, P4EST_QMAXLEVEL);
        l = &temp2;
      }
    }
    if (f != NULL) {
      P4EST_ASSERT (l != NULL);
      p4est_nearest_common_ancestor (f, l, &root);
    }

    /* grab complete tree quadrant array */
    tree = p4est_tree_array_index (p4est->trees, jt);
    tquadrants = &tree->quadrants;

    /* mark all points as active */
    sc_array_init_size (&actives, sizeof (size_t), points->elem_count);
    for (zz = 0; zz < points->elem_count; ++zz) {
      pz = (size_t *) sc_array_index (&actives, zz);
      *pz = zz;
    }

    p4est_search_recursion (p4est, jt, &root, search_quadrant_fn,
                            search_point_fn, tquadrants, points, &actives);
    sc_array_reset (&actives);
  }
}
