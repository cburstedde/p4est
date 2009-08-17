/*
  This file is part of p4est.
  p4est is a C library to manage a parallel collection of quadtrees and/or
  octrees.

  Copyright (C) 2009 Carsten Burstedde, Lucas Wilcox,
                     Toby Isaac.

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
#include <p8est_iterate.h>
#else
#include <p4est_algorithms.h>
#include <p4est_bits.h>
#include <p4est_iterate.h>
#endif

#define P4EST_ITER_STRIDE (P4EST_CHILDREN + 1)
typedef struct p4est_iter_tier
{
  p4est_quadrant_t   *key;
  size_t              array[P4EST_ITER_STRIDE];
}
p4est_iter_tier_t;

typedef struct p4est_iter_tier_ring
{
  int                 next;
  sc_array_t          tiers;
}
p4est_iter_tier_ring_t;

static void
p4est_iter_tier_update (sc_array_t * view, int level, size_t * next_tier,
                        size_t shift)
{
  int                 i;
  p4est_split_array (view, level, next_tier);
  for (i = 0; i < P4EST_ITER_STRIDE; i++) {
    next_tier[i] += shift;
  }
}

static void
p4est_iter_tier_insert (sc_array_t * view, int level, size_t * next_tier,
                        size_t shift, sc_array_t * tier_rings,
                        p4est_quadrant_t * q)
{
  int                 i, limit;
  p4est_iter_tier_ring_t *ring;
  p4est_iter_tier_t  *tier;
  p4est_quadrant_t   *key;

  if (q == NULL) {
    for (i = 0; i < P4EST_ITER_STRIDE; i++) {
      next_tier[i] = shift;
    }
    return;
  }

  if (level >= (int) tier_rings->elem_count) {
    p4est_iter_tier_update (view, level, next_tier, shift);
    return;
  }
  ring = sc_array_index_int (tier_rings, level);

  limit = (int) ring->tiers.elem_count;
  for (i = 0; i < limit; i++) {
    tier = sc_array_index_int (&(ring->tiers), i);
    key = tier->key;
    if (key == NULL) {
      if (ring->next != i) {
        P4EST_VERBOSEF ("level %d i %d next %d limit %d\n",
                        level, i, ring->next, limit);
      }
      P4EST_ASSERT (ring->next == i);
      p4est_iter_tier_update (view, level, next_tier, shift);
      memcpy (tier->array, next_tier, P4EST_ITER_STRIDE * sizeof (size_t));
      tier->key = q;
      ring->next++;
      ring->next %= limit;
      return;
    }
    if (q == key) {
      memcpy (next_tier, tier->array, P4EST_ITER_STRIDE * sizeof (size_t));
      return;
    }
  }

  /* if the tier wasn't already computed, compute it */
  p4est_iter_tier_update (view, level, next_tier, shift);
  /* we always rewrite over the oldest created tier */
  tier = sc_array_index_int (&(ring->tiers), ring->next++);
  memcpy (tier->array, next_tier, P4EST_ITER_STRIDE * sizeof (size_t));
  tier->key = q;
  ring->next %= limit;
}

typedef struct p4est_iter_corner_args
{
  int                 level;    /* the level of the initial search areas */
  int                 num_sides;
  int                *start_idx2;       /* for each side, the ancestor_id at
                                           level of the initial search area */
  int                *quad_idx2;        /* an indexing variable used in
                                           corner_iterate: passed as an
                                           argument to avoid using alloc/free
                                           on each call */
  sc_array_t        **quadrants;        /* the arrays, two for each side (one
                                           local, one ghost), that contain the
                                           quadrants in each search area */
  size_t            **index;    /* for each sidetype, the indices in quadrants
                                   that form the bounds of the heirarchical
                                   search areas */
  size_t             *first_index;      /* an indexing variable used in
                                           corner_iterate: passed as an
                                           argument to avoid using alloc/free
                                           on each call */
  size_t             *count;    /* a counting variable used in
                                   corner_iterate: passed as an argument to
                                   avoid using alloc/free on each call */
  p4est_quadrant_t  **test;     /* a testing variable used in
                                   corner_iterate: passed as an argument to
                                   avoid using alloc/free on each call */
  p4est_iter_corner_info_t *info;
}
p4est_iter_corner_args_t;

int
p4est_quadrant_compare_contains (const void *a, const void *b)
{
  const p4est_quadrant_t *q = (p4est_quadrant_t *) a;
  const p4est_quadrant_t *r = (p4est_quadrant_t *) b;
  int8_t              level = (q->level < r->level) ? q->level : r->level;
  p4est_qcoord_t      mask =
    ((p4est_qcoord_t) - 1) << (P4EST_MAXLEVEL - level);

  if (((q->x ^ r->x) & mask) || ((q->y ^ r->y) & mask)
#ifdef P4_TO_P8
      || ((q->z ^ r->z) & mask)
#endif
    ) {
    return p4est_quadrant_compare (a, b);
  }

  return 0;
}

static void
p4est_corner_iterate (p4est_iter_corner_args_t * args, void *user_data,
                      p4est_iter_corner_t iter_corner)
{
  const int           local = 0;
  const int           ghost = 1;

  int                 side, sidetype;
  int                 level = args->level;
  int                 num_sides = args->num_sides;
  const int          *start_idx2 = args->start_idx2;
  int                *quad_idx2 = args->quad_idx2;
  int                 this_corner;
  sc_array_t        **quadrants = args->quadrants;
  size_t            **index = args->index;
  size_t             *first_index = args->first_index;
  size_t             *count = args->count;
  p4est_quadrant_t  **test = args->test;
  p4est_quadrant_t    temp;
  p4est_qcoord_t      mask =
    ((p4est_qcoord_t) - 1) << (P4EST_MAXLEVEL - level);
  sc_array_t          test_view;
  p4est_iter_corner_info_t *info = args->info;
  p4est_iter_corner_side_t *cside;
  ssize_t             temp_idx;
  int                 level_idx2;
  int                 type;
  bool                has_local;

  /* level_idx2 moves us to the correct set of bounds within the index arrays
   * for the level: it is a set of bounds because it includes all children at
   * this level */
  level_idx2 = level * P4EST_ITER_STRIDE;

  for (side = 0; side < num_sides; side++) {
    /* start_idx2 gives the ancestor id at level for the search area on this
     * side, so quad_idx2[side] now gives the correct location in
     * index[sidetype] of the bounds of the search area */
    quad_idx2[side] = level_idx2 + start_idx2[side];

    /* get the location in quadrants[sidetype] of the first quadrant in the
     * search area, and the count of quadrants in the search area, and
     * initialize tests to NULL */
    for (type = local; type <= ghost; type++) {
      sidetype = side * 2 + type;
      first_index[sidetype] = index[sidetype][quad_idx2[side]];
      count[sidetype] = (index[sidetype][quad_idx2[side] + 1] -
                         first_index[sidetype]);
      test[sidetype] = NULL;
    }
  }

  /* corner_iterate only runs if there is a chance of a local quadrant touching
   * the desired corner */
  for (side = 0; side < num_sides; side++) {
    if (count[side * 2 + local]) {
      break;
    }
  }
  if (side == num_sides) {
    return;
  }

  has_local = false;
  for (side = 0; side < num_sides; side++) {

    cside = sc_array_index_int (&(info->sides), side);
    cside->quad = NULL;
    cside->is_local = false;
    cside->quadid = -1;
    this_corner = cside->corner;
    for (type = local; type <= ghost; type++) {
      sidetype = side * 2 + type;
      /* if we already found something locally, there's no need to search the
       * ghost layer */
      if (test[side * 2 + local] != NULL) {
        continue;
      }

      /* for this sidetype, we must find the most likely candidate in the
       * search area for touching the desired corner */
      if (count[sidetype]) {
        /* get a candidate */
        if (count[sidetype] == 1) {
          test[sidetype] = sc_array_index (quadrants[sidetype],
                                           first_index[sidetype]);
          temp_idx = 0;
        }
        else {
          switch (this_corner) {
          case (P4EST_CHILDREN - 1):
            P4EST_ASSERT (first_index[sidetype] + count[sidetype] - 1 <
                          quadrants[sidetype]->elem_count);
            test[sidetype] = sc_array_index (quadrants[sidetype],
                                             first_index[sidetype] +
                                             count[sidetype] - 1);
            temp_idx = ((ssize_t) count[sidetype]) - 1;
            break;
          default:
            P4EST_ASSERT (first_index[sidetype] <
                          quadrants[sidetype]->elem_count);
            test[sidetype] = sc_array_index (quadrants[sidetype],
                                             first_index[sidetype]);
            temp_idx = 0;
            break;
          }
        }
        /* create the smallest quadrant in the appropriate corner */
        temp = *(test[sidetype]);
        temp.x &= mask;
        temp.y &= mask;
#ifdef P4_TO_P8
        temp.z &= mask;
#endif
        temp.level = P4EST_QMAXLEVEL;
        P4EST_ASSERT (p4est_quadrant_is_valid (&temp));
        temp.x += (this_corner % 2) ?
          P4EST_QUADRANT_LEN (level) -
          P4EST_QUADRANT_LEN (P4EST_QMAXLEVEL) : 0;
        temp.y += ((this_corner % 4) / 2) ?
          P4EST_QUADRANT_LEN (level) -
          P4EST_QUADRANT_LEN (P4EST_QMAXLEVEL) : 0;
#ifdef P4_TO_P8
        temp.z += (this_corner / 4) ?
          P4EST_QUADRANT_LEN (level) -
          P4EST_QUADRANT_LEN (P4EST_QMAXLEVEL) : 0;
#endif
        P4EST_ASSERT (p4est_quadrant_is_valid (&temp));
        /* we do not have to search if there is one quadrant, or if we are in
         * the first or last corner */
        if (count[sidetype] == 1 || this_corner == 0 ||
            this_corner == P4EST_CHILDREN - 1) {
          /* if test[sidetype] does not contain temp */
          if ((!p4est_quadrant_is_equal (test[sidetype], &temp)) &&
              (!p4est_quadrant_is_ancestor (test[sidetype], &temp))) {
            test[sidetype] = NULL;
          }
        }
        else {
          /* we search for the last quadrant before the temp */
          sc_array_init_view (&test_view, quadrants[sidetype],
                              first_index[sidetype], count[sidetype]);
          temp_idx = sc_array_bsearch (&test_view, &temp,
                                       p4est_quadrant_compare_contains);
          /* if there is no quadrant containing temp, then no quad in the
           * search area can touch the corner */
          if (temp_idx == -1) {
            test[sidetype] = NULL;
          }
          else {
            P4EST_ASSERT ((size_t) temp_idx < test_view.elem_count);
            test[sidetype] = sc_array_index_ssize_t (&test_view, temp_idx);
          }
        }
        /* if we have found the right quadrant for this side of the corner */
        if (test[sidetype] != NULL) {
          P4EST_ASSERT (p4est_quadrant_compare_contains
                        (test[sidetype], &temp) == 0);
          P4EST_ASSERT (temp_idx >= 0 && (size_t) temp_idx < count[sidetype]);
          temp_idx += first_index[sidetype];
          cside->quad = test[sidetype];
          cside->is_local = (type == local);
          cside->quadid = (p4est_locidx_t) temp_idx;
          if (type == local) {
            has_local = true;
          }
        }
      }
    }
  }
  /* if none of the quads touching the corner is local, nothing is done */
  if (!has_local) {
    return;
  }

  /* run the callback */
  iter_corner (info, user_data);
}

#ifdef P4_TO_P8
typedef struct p8est_iter_edge_args
{
  int                 num_sides;
  int                 level;    /* the level of the initial search areas */
  int                *start_idx2;       /* for each side, the ancestor_id at
                                           level of the initial search area */
  int                *level_num;        /* an array that keeps track of which
                                           branch we take at each step in the
                                           heirarchical search areas */
  sc_array_t         *common_corners[2];        /* for each side of the edge,
                                                   there are two corners that
                                                   touch the edge */
  sc_array_t        **quadrants;        /* the arrays, two for each side (one
                                           local, one ghost), that contain the
                                           quadrants in each search area */
  size_t            **index;    /* for each sidetype, the indices in quadrants
                                   that form the bounds of the heirarchical
                                   search areas */
  size_t             *first_index;      /* an indexing variable used in
                                           edge_iterate: passed as an
                                           argument to avoid using alloc/free
                                           on each call */
  p8est_quadrant_t  **test;     /* a testing variable used in
                                   edge_iterate: passed as an argument to
                                   avoid using alloc/free on each call */
  size_t             *count;    /* a counting variable used in
                                   edge_iterate: passed as an argument to
                                   avoid using alloc/free on each call */
  int                *test_level;       /* a testing variable used in
                                           edge_iterate: passed as an argument
                                           to avoid using alloc/free on each 
                                           call */
  int                *quad_idx2;        /* an indexing variable used in
                                           edge_iterate: passed as an
                                           argument to avoid using alloc/free
                                           on each call */
  bool               *refine;   /* a testing variable used in edge_iterate:
                                   passed as an argument to avoid using
                                   alloc/free on each call */
  p4est_iter_corner_args_t *corner_args;
  sc_array_t         *tier_rings;
  p8est_iter_edge_info_t *info;
}
p8est_iter_edge_args_t;

static void
p8est_edge_iterate (p8est_iter_edge_args_t * args, void *user_data,
                    p8est_iter_edge_t iter_edge,
                    p8est_iter_corner_t iter_corner)
{
  const int           local = 0;
  const int           ghost = 1;

  int                 num_sides = args->num_sides;
  int                 start_level = args->level;
  int                *start_idx2 = args->start_idx2;
  int                *level_num = args->level_num;
  sc_array_t        **quadrants = args->quadrants;
  size_t            **index = args->index;
  size_t             *first_index = args->first_index;
  sc_array_t        **common_corners = args->common_corners;
  p8est_quadrant_t  **test = args->test;
  size_t             *count = args->count;
  int                *test_level = args->test_level;
  int                *quad_idx2 = args->quad_idx2;
  bool               *refine = args->refine;
  int                *temp_int, *temp_int2;
  int                 i;
  int                 level;
  int                 side;
  int                 type;
  int                 sidetype;
  int                 level_idx2;
  p8est_iter_edge_info_t *info = args->info;
  p8est_iter_edge_side_t *eside;
  p8est_quadrant_t  **quads;
  p4est_locidx_t     *quadids;
  bool               *is_local;
  int                 child_corner;
  bool                has_local;
  sc_array_t          test_view;
  bool                all_empty, stop_refine;
  int                 temp_idx2;
  p4est_iter_corner_args_t *corner_args = args->corner_args;
  sc_array_t         *corner_sides = &(args->corner_args->info->sides);
  p4est_iter_corner_side_t *cside;
  sc_array_t         *tier_rings = args->tier_rings;

  /* level_idx2 moves us to the correct set of bounds within the index arrays
   * for the level: it is a set of bounds because it includes all children at
   * this level */
  level = start_level;
  level_idx2 = level * P4EST_ITER_STRIDE;
  for (side = 0; side < num_sides; side++) {

    /* start_idx2 gives the ancestor id at level for the search area on this
     * side, so quad_idx2[side] now gives the correct location in
     * index[sidetype] of the bounds of the search area */
    quad_idx2[side] = level_idx2 + start_idx2[side];

    /* get the location in quadrants[sidetype] of the first quadrant in the
     * search area, and the count of quadrants in the search area */
    for (type = local; type <= ghost; type++) {
      sidetype = side * 2 + type;
      first_index[sidetype] = index[sidetype][quad_idx2[side]];
      count[sidetype] = (index[sidetype][quad_idx2[side] + 1] -
                         first_index[sidetype]);
    }
  }

  /* edge_iterate only runs if there is a chance of a local quadrant touching
   * the desired edge */
  for (side = 0; side < num_sides; side++) {
    if (count[side * 2 + local]) {
      break;
    }
  }
  if (side == num_sides) {
    return;
  }

  /* we think of the search tree as being rooted at start_level, so we can
   * think the branch number at start_level as 0, even if it actually is not */
  level_num[start_level] = 0;

  for (;;) {
    /* for each sidetype, get the first quadrant in that sidetype search area
     */
    for (side = 0; side < num_sides; side++) {
      for (type = local; type <= ghost; type++) {
        sidetype = side * 2 + type;
        if (count[sidetype]) {
          P4EST_ASSERT (first_index[sidetype] <
                        quadrants[sidetype]->elem_count);
          test[sidetype] = sc_array_index (quadrants[sidetype],
                                           first_index[sidetype]);
          test_level[sidetype] = (int) test[sidetype]->level;
        }
        else {
          test[sidetype] = NULL;
          test_level[sidetype] = -1;
        }
      }
      /* initially assume that every side needs to be refined */
      refine[side] = true;
    }
    /* initially assume that we are going to have to refine our search areas */
    stop_refine = false;
    has_local = false;
    for (side = 0; side < num_sides; side++) {
      for (type = local; type <= ghost; type++) {
        sidetype = side * 2 + type;
        /* if the candidate from sidetype is the same size as the search area,
         * then we do not refine this side */
        if (test_level[sidetype] == level) {
          if (iter_edge != NULL) {
            refine[side] = false;
            /* by the two to one condition, we do not need to continue the
             * search beyond the possibility of neighbors to this quad being
             * one size smaller */
            stop_refine = true;
            eside = sc_array_index_int (&(info->sides), side);
            eside->is_hanging = false;
            eside->is.full.quad = test[sidetype];
            eside->is.full.is_local = (type == local);
            eside->is.full.quadid = first_index[sidetype];
            has_local = (has_local || (type == local));
          }
          /* if there is no edge callback (i.e., we are just running
           * edge_iterate to find corners), then we're done with this branch */
          else {
            level_num[level]++;
            goto change_search_area;
          }
        }
      }
    }
    /* if no side needs to be refined, then we run the iter_edge and proceed to
     * the next search area on this level */
    for (side = 0; side < num_sides; side++) {
      if (refine[side]) {
        break;
      }
    }
    if (side == num_sides) {
      iter_edge (info, user_data);
      level_num[level]++;
      goto change_search_area;
    }

    /* at this point, some sides need to be refined, so we take the search area
     * and split it up, taking the indices for the refined search areas and
     * placing them on the next tier in index[sidetype] */
    for (side = 0; side < num_sides; side++) {
      if (refine[side]) {
        quad_idx2[side] = level_idx2 + P4EST_ITER_STRIDE;
        for (type = local; type <= ghost; type++) {
          sidetype = side * 2 + type;
          sc_array_init_view (&test_view, quadrants[sidetype],
                              first_index[sidetype], count[sidetype]);
          p4est_iter_tier_insert (&test_view, level, index[sidetype] +
                                  quad_idx2[side], first_index[sidetype],
                                  tier_rings, test[sidetype]);
        }
      }
    }

    /* if one of sides was not refined, then it is time to run iter_edge */
    if (stop_refine) {
      /* for each corner of the common edge, college all of the quadrants that
       * touch that corner and run iter_edge with that collection */
      for (side = 0; side < num_sides; side++) {
        if (refine[side]) {
          eside = sc_array_index_int (&(info->sides), side);
          eside->is_hanging = true;
          quads = eside->is.hanging.quad;
          is_local = eside->is.hanging.is_local;
          quadids = eside->is.hanging.quadid;
          for (i = 0; i < 2; i++) {
            temp_int = sc_array_index_int (common_corners[i], side);
            temp_int2 = sc_array_index_int (common_corners[1 - i], side);
            if (*temp_int < *temp_int2) {
              P4EST_ASSERT (p8est_edge_corners[eside->edge][0] == *temp_int);
              child_corner = 0;
            }
            else {
              P4EST_ASSERT (p8est_edge_corners[eside->edge][1] == *temp_int);
              child_corner = 1;
            }
            quads[child_corner] = NULL;
            is_local[child_corner] = false;
            quadids[child_corner] = -1;

            quad_idx2[side] = level_idx2 + P4EST_ITER_STRIDE + *temp_int;
            for (type = local; type <= ghost; type++) {
              sidetype = side * 2 + type;
              first_index[sidetype] = index[sidetype][quad_idx2[side]];
              count[sidetype] =
                (size_t) index[sidetype][quad_idx2[side] + 1] -
                first_index[sidetype];
              /* if the search area is non-empty, by the two to one condition
               * it must contain exactly one quadrant, which we add to the
               * collection */
              if (count[sidetype]) {
                P4EST_ASSERT (count[sidetype] == 1);
                P4EST_ASSERT (first_index[sidetype] <
                              quadrants[sidetype]->elem_count);
                quads[child_corner] = sc_array_index (quadrants[sidetype],
                                                      first_index[sidetype]);
                P4EST_ASSERT ((int) quads[child_corner]->level == level + 1);
                is_local[child_corner] = (type == local);
                quadids[child_corner] =
                  (p4est_locidx_t) first_index[sidetype];
                has_local = (has_local || (type == local));
              }
            }
          }
        }
      }
      if (has_local) {
        iter_edge (info, user_data);
      }
      /* we proceed to the next search area (i.e., branch) on this level */
      level_num[level]++;
      goto change_search_area;
    }
    /* if every side needed to be refined, then we descend along this branch to
     * this next level and search there */
    level_num[++level] = 0;
    level_idx2 += P4EST_ITER_STRIDE;

  change_search_area:
    /* if we tried to advance the search area on start_level, we've completed
     * the search */
    if (level_num[start_level] > 0) {
      break;
    }
    /* if we have tried to advance the search area past two branches, that
     * means that we have completed all of the branches on this level */
    if (level_num[level] == 2) {
      /* if we have a corner callback, we need to run it on the corner between
       * the edge branches on this level */
      if (iter_corner != NULL) {
        corner_args->level = level;
        /* find the correct corner ids for the corner search areas, and copy
         * the index information necessary to run corner_iterate over to the
         * additional index arrays */
        for (i = 0; i < 2 * num_sides; i++) {
          cside = sc_array_index_int (corner_sides, i);
          temp_int = sc_array_index_int (common_corners[1 - i / num_sides],
                                         i % num_sides);
          cside->corner = *temp_int;
          temp_int = sc_array_index_int (common_corners[i / num_sides],
                                         i % num_sides);
          start_idx2[i] = *temp_int;
          temp_idx2 = level_idx2 + start_idx2[i];
          index[i * 2 + local][temp_idx2] =
            index[(i % num_sides) * 2 + local][temp_idx2];
          index[i * 2 + local][temp_idx2 + 1] =
            index[(i % num_sides) * 2 + local][temp_idx2 + 1];
          index[i * 2 + ghost][temp_idx2] =
            index[(i % num_sides) * 2 + ghost][temp_idx2];
          index[i * 2 + ghost][temp_idx2 + 1] =
            index[(i % num_sides) * 2 + ghost][temp_idx2 + 1];
        }
        p4est_corner_iterate (corner_args, user_data, iter_corner);
      }
      /* now that we're done on this level, go up a level and over a branch */
      level_num[--level]++;
      level_idx2 -= P4EST_ITER_STRIDE;
      goto change_search_area;
    }

    /* at this point, we need to initialize the bounds of the search areas
     * for this new branch */
    all_empty = true;
    for (side = 0; side < num_sides; side++) {
      P4EST_ASSERT ((size_t) side <
                    common_corners[level_num[level]]->elem_count);
      temp_int = sc_array_index_int (common_corners[level_num[level]], side);
      quad_idx2[side] = level_idx2 + *temp_int;
      for (type = local; type <= ghost; type++) {
        sidetype = side * 2 + type;
        first_index[sidetype] = index[sidetype][quad_idx2[side]];
        count[sidetype] = (index[sidetype][quad_idx2[side] + 1] -
                           first_index[sidetype]);
        if (type == local && count[sidetype]) {
          all_empty = false;
        }
      }
    }
    /* if there are no local quadrants in any of the search areas, we're done
     * with this search area and proceed to the next branch on this level */
    if (all_empty) {
      level_num[level]++;
      goto change_search_area;
    }
  }
}
#endif

typedef struct p4est_iter_face_args
{
  int                 level;    /* the level of the initial search areas */
  int                *start_idx2;       /* for each side, the ancestor_id at
                                           level of the initial search area */
  int                *level_num;        /* an array that keeps track of which
                                           branch we take at each step in the
                                           heirarchical search areas */

  int                *num_to_child;     /* when a search branch is refined,
                                           num_to_child says which child id
                                           corresponds to the branch number for
                                           each side of the face. e.g. Suppose
                                           face[left] = 1, face[right] = 0, and
                                           orientation = 0 in 3D. The child ids
                                           of the descendants of the current
                                           search area that touch face[left]
                                           are 1, 3, 5, 7, and given
                                           face[right] and the orientation, the
                                           descendants that are opposite them
                                           are 0, 2, 4, 6, respectively:
                                           therefore num_to_child =
                                           { 1, 3, 5, 7, 0, 2, 4, 6} */

  sc_array_t        **quadrants;        /* the arrays, two for each side (one
                                           local, one ghost), that contain the
                                           quadrants in each search area */
  size_t            **index;    /* for each sidetype, the indices in quadrants
                                   that form the bounds of the heirarchical
                                   search areas */
  size_t             *first_index;      /* an indexing variable used in
                                           face_iterate: passed as an
                                           argument to avoid using alloc/free
                                           on each call */
  p4est_quadrant_t  **test;     /* a testing variable used in
                                   face_iterate: passed as an argument to
                                   avoid using alloc/free on each call */
  size_t             *count;    /* a counting variable used in
                                   face_iterate: passed as an argument to
                                   avoid using alloc/free on each call */
  int                *test_level;       /* a testing variable used in
                                           edge_iterate: passed as an argument
                                           to avoid using alloc/free on each 
                                           call */
  int                *quad_idx2;        /* an indexing variable used in
                                           edge_iterate: passed as an
                                           argument to avoid using alloc/free
                                           on each call */
  bool               *refine;   /* a testing variable used in edge_iterate:
                                   passed as an argument to avoid using
                                   alloc/free on each call */
  bool                outside_face;     /* indicates if we are at a tree
                                           boundary without a neighbor across
                                           the face */
#ifdef P4_TO_P8
  p8est_iter_edge_args_t *edge_args;
#endif
  p4est_iter_corner_args_t *corner_args;
  sc_array_t         *tier_rings;
  p4est_iter_face_info_t *info;
}
p4est_iter_face_args_t;

static void
p4est_face_iterate (p4est_iter_face_args_t * args, void *user_data,
                    p4est_iter_face_t iter_face,
#ifdef P4_TO_P8
                    p8est_iter_edge_t iter_edge,
#endif
                    p4est_iter_corner_t iter_corner)
{

  const int           left = 0;
  const int           right = 1;
  const int           local = 0;
  const int           ghost = 1;
  const int           ntc_str = P4EST_CHILDREN / 2;

  int                 start_level = args->level;
  int                *start_idx2 = args->start_idx2;
  int                *level_num = args->level_num;
  sc_array_t        **quadrants = args->quadrants;
  size_t            **index = args->index;
  size_t             *first_index = args->first_index;
  int                *num_to_child = args->num_to_child;
  p4est_quadrant_t  **test = args->test;
  size_t             *count = args->count;
  int                *test_level = args->test_level;
  int                *quad_idx2 = args->quad_idx2;
  bool               *refine = args->refine;
  int                 limit;

  int                 i, j;
  int                 level;
  int                 side, n_side;
  int                 type, n_type;
  int                 sidetype, nsidentype;
  int                 level_idx2;
  p4est_iter_face_info_t *info = args->info;
  p4est_iter_face_side_t *fside;
  p4est_quadrant_t  **quads;
  p4est_locidx_t     *quadids;
  bool               *is_local;
  int                 child_corner;
  bool                has_local;
  sc_array_t          test_view;
  int                 temp_idx2;
  p4est_iter_corner_args_t *corner_args = args->corner_args;
  sc_array_t         *tier_rings = args->tier_rings;
  sc_array_t         *corner_sides = &(args->corner_args->info->sides);
  p4est_iter_corner_side_t *cside;
#ifdef P4_TO_P8
  int                *c;
  int                 v1, v2, true_dir, dir, k;
  p8est_iter_edge_args_t *edge_args = args->edge_args;
  sc_array_t         *common_corners[2];

  sc_array_t         *edge_sides = &(args->edge_args->info->sides);
  p8est_iter_edge_side_t *eside;

  common_corners[0] = edge_args->common_corners[0];
  common_corners[1] = edge_args->common_corners[1];
#endif

  /* if we are at an outside face, then there is no right half to our search
   * that needs to be coordinated with the left half */
  limit = args->outside_face ? left : right;

  /* level_idx2 moves us to the correct set of bounds within the index arrays
   * for the level: it is a set of bounds because it includes all children at
   * this level */
  level = start_level;
  level_idx2 = level * P4EST_ITER_STRIDE;

  for (side = left; side <= limit; side++) {

    /* start_idx2 gives the ancestor id at level for the search area on this
     * side, so quad_idx2[side] now gives the correct location in
     * index[sidetype] of the bounds of the search area */
    quad_idx2[side] = level_idx2 + start_idx2[side];

    /* get the location in quadrants[sidetype] of the first quadrant in the
     * search area, and the count of quadrants in the search area */
    for (type = local; type <= ghost; type++) {
      sidetype = side * 2 + type;
      first_index[sidetype] = index[sidetype][quad_idx2[side]];
      count[sidetype] = (index[sidetype][quad_idx2[side] + 1] -
                         first_index[sidetype]);
    }
  }

  /* face_iterate only runs if there is a chance of a local quadrant touching
   * the desired face */
  if (!args->outside_face) {
    if (!count[left * 2 + local] && !count[right * 2 + local]) {
      return;
    }
  }
  else {
    if (!count[left * 2 + local]) {
      return;
    }
  }

  /* we think of the search tree as being rooted at start_level, so we can
   * think the branch number at start_level as 0, even if it actually is not */
  level_num[start_level] = 0;
  for (;;) {
    /* for each sidetype, get the first quadrant in that sidetype search area
     */
    for (side = left; side <= limit; side++) {
      for (type = local; type <= ghost; type++) {
        sidetype = side * 2 + type;
        if (count[sidetype]) {
          P4EST_ASSERT (first_index[sidetype] <
                        quadrants[sidetype]->elem_count);
          test[sidetype] = sc_array_index (quadrants[sidetype],
                                           first_index[sidetype]);
          test_level[sidetype] = (int) test[sidetype]->level;
        }
        else {
          test[sidetype] = NULL;
          test_level[sidetype] = -1;
        }
      }
    }
    /* initially assume that each side needs to be refined */
    refine[left] = refine[right] = true;
    has_local = false;

    /* get a candidate from each sidetype */
    for (side = left; side <= limit; side++) {
      for (type = local; type <= ghost; type++) {
        sidetype = side * 2 + type;
        /* if the candidate from sidetype is the same size as the search area,
         * then we do not refine this side */
        if (test_level[sidetype] == level) {
          if (iter_face != NULL) {
            P4EST_ASSERT (count[sidetype] == 1);
            P4EST_ASSERT (count[side * 2 + (type ^ 1)] == 0);
            fside = sc_array_index_int (&(info->sides), side);
            fside->is_hanging = false;
            fside->is.full.quad = test[sidetype];
            fside->is.full.quadid = (p4est_locidx_t) first_index[sidetype];
            refine[side] = false;
            has_local = (type == local);
            fside->is.full.is_local = has_local;
            if (!args->outside_face) {
              n_side = side ^ 1;
              for (n_type = type; n_type <= ghost; n_type++) {
                nsidentype = n_side * 2 + n_type;
                /* if the quadrant opposite our candidate is the same size,
                 * then we run iter_face and proceed to the next branch on this
                 * level */
                if ((n_type > type || n_side > side) &&
                    test_level[nsidentype] == level) {
                  P4EST_ASSERT (count[nsidentype] == 1);
                  P4EST_ASSERT (count[n_side * 2 + (n_type ^ 1)] == 0);
                  P4EST_ASSERT (!(type == ghost && n_type == ghost));
                  fside = sc_array_index_int (&(info->sides), n_side);
                  fside->is_hanging = false;
                  fside->is.full.quad = test[nsidentype];
                  fside->is.full.quadid =
                    (p4est_locidx_t) first_index[nsidentype];
                  fside->is.full.is_local = (n_type == local);
                  iter_face (info, user_data);
                  level_num[level]++;
                  goto change_search_area;
                }
              }
#ifdef P4EST_DEBUG
              if (type == local) {
                P4EST_ASSERT (count[n_side * 2 + type] > 0 ||
                              count[n_side * 2 + (type ^ 1)] > 0);
              }
#endif
            }
            else {
              iter_face (info, user_data);
              level_num[level]++;
              goto change_search_area;
            }
          }
          /* if there is no face callback, i.e. we are running face_iterate only
           * to find the edges/corners that live on faces, then we are done once
           * we find a side that does not need to be refined */
          else {
            level_num[level]++;
            goto change_search_area;
          }
        }
      }
    }

    /* if a side needs to be refined, we take the search area and split it up,
     * taking the indices for the refined search areas and placing them on the
     * next tier in index[sidetype] */
    for (side = left; side <= limit; side++) {
      if (refine[side]) {
        quad_idx2[side] = level_idx2 + P4EST_ITER_STRIDE;
        for (type = local; type <= ghost; type++) {
          sidetype = side * 2 + type;
          sc_array_init_view (&test_view, quadrants[sidetype],
                              first_index[sidetype], count[sidetype]);
          p4est_iter_tier_insert (&test_view, level, index[sidetype] +
                                  quad_idx2[side], first_index[sidetype],
                                  tier_rings, test[sidetype]);
        }
      }
    }

    for (side = left; side <= limit; side++) {
      /* if this side was not refined, then we need to run iter_face with the
       * candidate from this side and each of its hanging neighbors */
      if (!refine[side]) {
        n_side = side ^ 1;
        fside = sc_array_index_int (&(info->sides), n_side);
        fside->is_hanging = true;
        quads = fside->is.hanging.quad;
        quadids = fside->is.hanging.quadid;
        is_local = fside->is.hanging.is_local;
        for (i = 0; i < P4EST_CHILDREN / 2; i++) {
          child_corner = num_to_child[n_side * ntc_str + i];
#ifndef P4_TO_P8
          if (child_corner < num_to_child[n_side * ntc_str + (1 - i)]) {
            child_corner = 0;
          }
          else {
            child_corner = 1;
          }
#else
          child_corner = p8est_corner_face_corners[child_corner][fside->face];
#endif
          quads[child_corner] = NULL;
          quadids[child_corner] = -1;
          is_local[child_corner] = false;
          quad_idx2[n_side] = level_idx2 + P4EST_ITER_STRIDE +
            num_to_child[n_side * ntc_str + i];
          for (n_type = local; n_type <= ghost; n_type++) {
            nsidentype = n_side * 2 + n_type;
            first_index[nsidentype] = index[nsidentype][quad_idx2[n_side]];
            count[nsidentype] = index[nsidentype][quad_idx2[n_side] + 1] -
              first_index[nsidentype];
            /* if the search area is non-empty, by the two to one condition
             * it must contain exactly one quadrant: if one of the two types
             * is local, we run iter_face */
            if (count[nsidentype]) {
              P4EST_ASSERT (count[nsidentype] == 1);
              P4EST_ASSERT (index[n_side * 2 + (n_type ^ 1)]
                            [quad_idx2[n_side]]
                            == index[n_side * 2 + (n_type ^ 1)]
                            [quad_idx2[n_side] + 1]);
              quads[child_corner] = sc_array_index (quadrants[nsidentype],
                                                    first_index[nsidentype]);
              quadids[child_corner] =
                (p4est_locidx_t) first_index[nsidentype];
              is_local[child_corner] = (n_type == local);
              if (n_type == local) {
                has_local = true;
              }
            }
          }
        }
        if (has_local) {
          iter_face (info, user_data);
        }
        /* once we've run iter_face for the hanging faces, we
         * proceed to the next branch on this level */
        level_num[level]++;
        goto change_search_area;
      }
    }
    /* if we refined both sides, we descend to the next level from this branch
     * and continue searching there */
    level_num[++level] = 0;
    level_idx2 += P4EST_ITER_STRIDE;

  change_search_area:
    /* if we tried to advance the search area on start_level, we've completed
     * the search */
    if (level_num[start_level] > 0) {
      break;
    }

    /* if we have tried to advance the search area past the number of
     * descendants, that means that we have completed all of the branches on
     * this level */
    if (level_num[level] == P4EST_CHILDREN / 2) {
#ifdef P4_TO_P8
      /* if we have an edge callback, we need to run it on all of the edges
       * between the face branches on this level */
      if (iter_edge != NULL || iter_corner != NULL) {
        edge_args->level = level;
        if (!args->outside_face) {
          for (dir = 0; dir < 2; dir++) {

            /* we have to set up the common corners on either side of the
             * desired edge */
            for (side = 0; side < 2; side++) {
              for (j = 0; j < 4; j++) {
                k = j >> 1;
                c = sc_array_index_int (common_corners[side], j);
                if (dir == 0) {
                  *c = num_to_child[(j % 2) * ntc_str + (2 * (1 - k) + side)];
                }
                else {
                  *c = num_to_child[(j % 2) * ntc_str + ((1 - k) + 2 * side)];
                }
              }
            }

            /* we have to set up the edge id for all four search areas that
             * touch the desired edge */
            for (j = 0; j < 4; j++) {
              c = sc_array_index_int (common_corners[0], j);
              v1 = *c;
              c = sc_array_index_int (common_corners[1], j);
              v2 = *c;
              true_dir = v2 - v1;
              true_dir = (true_dir >= 0) ? true_dir : -true_dir;
              eside = sc_array_index_int (edge_sides, j);
              eside->edge = p8est_corner_edges[v1][true_dir >> 1];
              eside->orientation = (v1 < v2) ? 0 : 1;
            }

            /* copy the index information necessary to run edge_iterate over to
             * the additional index arrays */
            for (side = 0; side < 2; side++) {
              for (j = 0; j < 4; j++) {
                k = j >> 1;
                if (dir == 0) {
                  start_idx2[j] =
                    num_to_child[(j % 2) * ntc_str + side + 2 * k];
                }
                else {
                  start_idx2[j] =
                    num_to_child[(j % 2) * ntc_str + 2 * side + k];
                }
                temp_idx2 = level_idx2 + start_idx2[j];
                index[j * 2 + local][temp_idx2] = index[(j % 2) * 2 + local]
                  [temp_idx2];
                index[j * 2 + local][temp_idx2 + 1] =
                  index[(j % 2) * 2 + local]
                  [temp_idx2 + 1];
                index[j * 2 + ghost][temp_idx2] = index[(j % 2) * 2 + ghost]
                  [temp_idx2];
                index[j * 2 + ghost][temp_idx2 + 1] =
                  index[(j % 2) * 2 + ghost]
                  [temp_idx2 + 1];
              }
              p8est_edge_iterate (edge_args, user_data, iter_edge,
                                  iter_corner);
            }
          }
        }
        /* if we are on on outside face, there are only two quadrants to an
         * edge instead of four, and we have to set up the arrays
         * differently */
        else {
          for (dir = 0; dir < 2; dir++) {

            /* we have to set up the common corners on either side of the
             * desired edge */
            for (side = 0; side < 2; side++) {
              for (j = 0; j < 2; j++) {
                c = sc_array_index_int (common_corners[side], j);
                if (dir == 0) {
                  *c = num_to_child[(2 * (1 - j) + side)];
                }
                else {
                  *c = num_to_child[((1 - j) + 2 * side)];
                }
              }
            }

            /* we have to set up the edge id for both search areas that
             * touch the desired edge */
            for (j = 0; j < 2; j++) {
              c = sc_array_index_int (common_corners[0], j);
              v1 = *c;
              c = sc_array_index_int (common_corners[1], j);
              v2 = *c;
              true_dir = v2 - v1;
              true_dir = (true_dir >= 0) ? true_dir : -true_dir;
              eside = sc_array_index_int (edge_sides, j);
              eside->edge = p8est_corner_edges[v1][true_dir >> 1];
              eside->orientation = (v1 < v2) ? 0 : 1;
            }

            /* copy the index information necessary to run edge_iterate over to
             * the additional index arrays */
            for (side = 0; side < 2; side++) {
              for (j = 0; j < 2; j++) {
                if (dir == 0) {
                  start_idx2[j] = num_to_child[side + 2 * j];
                }
                else {
                  start_idx2[j] = num_to_child[2 * side + j];
                }
                temp_idx2 = level_idx2 + start_idx2[j];
                index[j * 2 + local][temp_idx2] = index[local][temp_idx2];
                index[j * 2 + local][temp_idx2 + 1] = index[local]
                  [temp_idx2 + 1];
                index[j * 2 + ghost][temp_idx2] = index[ghost][temp_idx2];
                index[j * 2 + ghost][temp_idx2 + 1] = index[ghost]
                  [temp_idx2 + 1];
              }
              p8est_edge_iterate (edge_args, user_data, iter_edge,
                                  iter_corner);
            }
          }
        }
      }
#endif
      /* if we have a corner callback, we need to run it on the corner between
       * the face branches on this level */
      if (iter_corner != NULL) {
        corner_args->level = level;
        /* find the correct corner ids for the corner search areas, and copy
         * the index information necessary to run corner_iterate over to the
         * additional index arrays */
        if (!args->outside_face) {
          for (j = 0; j < P4EST_CHILDREN; j++) {
            cside = sc_array_index_int (corner_sides, j);
            cside->corner = num_to_child
              [(j % 2) * ntc_str + P4EST_CHILDREN / 2 - 1 - (j / 2)];
            start_idx2[j] = num_to_child[(j % 2) * ntc_str + (j / 2)];
            temp_idx2 = level_idx2 + start_idx2[j];
            index[j * 2 + local][temp_idx2] = index[(j % 2) * 2 + local]
              [temp_idx2];
            index[j * 2 + local][temp_idx2 + 1] = index[(j % 2) * 2 + local]
              [temp_idx2 + 1];
            index[j * 2 + ghost][temp_idx2] = index[(j % 2) * 2 + ghost]
              [temp_idx2];
            index[j * 2 + ghost][temp_idx2 + 1] = index[(j % 2) * 2 + ghost]
              [temp_idx2 + 1];
          }
        }
        else {
          for (j = 0; j < P4EST_CHILDREN / 2; j++) {
            cside = sc_array_index_int (corner_sides, j);
            cside->corner = num_to_child[P4EST_CHILDREN / 2 - 1 - j];
            start_idx2[j] = num_to_child[j];
            temp_idx2 = level_idx2 + start_idx2[j];
            index[j * 2 + local][temp_idx2] = index[local][temp_idx2];
            index[j * 2 + local][temp_idx2 + 1] = index[local]
              [temp_idx2 + 1];
            index[j * 2 + ghost][temp_idx2] = index[ghost][temp_idx2];
            index[j * 2 + ghost][temp_idx2 + 1] = index[ghost]
              [temp_idx2 + 1];
          }
        }
        p4est_corner_iterate (corner_args, user_data, iter_corner);
      }

      /* now that we're done on this level, go up a level and over a branch */
      level_num[--level]++;
      level_idx2 -= P4EST_ITER_STRIDE;
      goto change_search_area;
    }

    /* at this point, we need to initialize the bounds of the search areas
     * for this new branch */
    for (side = left; side <= limit; side++) {
      quad_idx2[side] =
        level_idx2 + num_to_child[side * ntc_str + level_num[level]];
    }
    for (side = left; side <= limit; side++) {
      for (type = local; type <= ghost; type++) {
        sidetype = side * 2 + type;
        first_index[sidetype] = index[sidetype][quad_idx2[side]];
        count[sidetype] = (index[sidetype][quad_idx2[side] + 1] -
                           first_index[sidetype]);
      }
    }

    /* if there are no local quadrants in either of the search areas, we're
     * done with this search area and proceed to the next branch on this
     * level */
    if (!args->outside_face) {
      if (!count[left * 2 + local] && !count[right * 2 + local]) {
        level_num[level]++;
        goto change_search_area;
      }
    }
    else {
      if (!count[left * 2 + local]) {
        level_num[level]++;
        goto change_search_area;
      }
    }
  }
}

void static
p4est_volume_iterate (p4est_t * p4est, sc_array_t * ghost_layer,
                      void *user_data, p4est_iter_volume_t iter_volume)
{
  p4est_topidx_t      t;
  p4est_topidx_t      first_local_tree = p4est->first_local_tree;
  p4est_topidx_t      last_local_tree = p4est->last_local_tree;
  sc_array_t         *trees = p4est->trees;
  p4est_tree_t       *tree;
  size_t              si, n_quads;
  sc_array_t         *quadrants;
  p4est_iter_volume_info_t info;

  info.p4est = p4est;
  info.ghost_layer = ghost_layer;

  for (t = first_local_tree; t <= last_local_tree; t++) {
    info.treeid = t;
    tree = p4est_array_index_topidx (trees, t);
    quadrants = &(tree->quadrants);
    n_quads = quadrants->elem_count;
    for (si = 0; si < n_quads; si++) {
      info.quad = sc_array_index (quadrants, si);
      info.quadid = si;
      iter_volume (&info, user_data);
    }
  }
}

static              size_t
p4est_split_array_which_tree (sc_array_t * array, size_t index, void *data)
{
  p4est_quadrant_t   *q = sc_array_index (array, index);
  return ((size_t) q->p.which_tree);
}

static void
p4est_split_ghost_layer_by_tree (sc_array_t * ghosts, size_t * tree_offsets,
                                 p4est_topidx_t num_trees)
{
  sc_array_t          view;
  sc_array_init_data (&view, tree_offsets, sizeof (size_t), num_trees + 1);
  sc_array_split (ghosts, &view, (size_t) num_trees,
                  p4est_split_array_which_tree, NULL);
}

void
p4est_iterate (p4est_t * p4est, sc_array_t * ghost_layer, void *user_data,
               p4est_iter_volume_t iter_volume, p4est_iter_face_t iter_face,
#ifdef P4_TO_P8
               p8est_iter_edge_t iter_edge,
#endif
               p4est_iter_corner_t iter_corner)
{
  /* constants */
  const int           left = 0;
  const int           right = 1;
  const int           local = 0;
  const int           ghost = 1;
  const int           ntc_str = P4EST_CHILDREN / 2;

  /* iterators */
  int                 i, j, k, side, type;
  int                 left_side, right_side;
  p4est_topidx_t      t, c;
#ifdef P4_TO_P8
  p4est_topidx_t      e;
  p4est_topidx_t      jt, kt;
#endif

  /* counts */
  size_t              num_ghosts, num_local;

  /* trees and quadrant arrays */
  sc_array_t         *trees = p4est->trees;
  p4est_tree_t       *tree;
  sc_array_t        **quadrants;
  size_t             *tree_first_ghost;

  /* topology */
  p4est_connectivity_t *conn = p4est->connectivity;
  const int8_t       *ttf = conn->tree_to_face;
  const p4est_topidx_t *ttt = conn->tree_to_tree;
  p4est_topidx_t      nt;
  size_t              global_num_trees;
  p4est_topidx_t      this_corner;
#ifndef P4_TO_P8
  const p4est_topidx_t *ttc = conn->tree_to_vertex;
  p4est_topidx_t      num_corners = conn->num_vertices;
  const p4est_topidx_t *ctt_offset = conn->vtt_offset;
  p4est_corner_transform_t *ct;
  sc_array_t          ctransforms, *cta = &ctransforms;
#else
  const p4est_topidx_t *ttc = conn->tree_to_corner;
  p4est_topidx_t      num_corners = conn->num_corners;
  const p4est_topidx_t *ctt_offset = conn->ctt_offset;
  const p4est_topidx_t *corner_to_tree = conn->corner_to_tree;
  const int8_t       *corner_to_corner = conn->corner_to_corner;

  const p4est_topidx_t *tte = conn->tree_to_edge;
  p4est_topidx_t      num_edges = conn->num_edges;
  const p4est_topidx_t *ett_offset = conn->ett_offset;
  const p4est_topidx_t *edge_to_tree = conn->edge_to_tree;
  p4est_topidx_t      this_edge;
  const int8_t       *edge_to_edge = conn->edge_to_edge;
#endif

  /* allocation size variables */
  int                 alloc_size;
  int                 max_corner_size;
  int                 corner_size;
#ifdef P4_TO_P8
  int                 max_edge_size;
  int                 edge_size;
#endif

  /* search loop variables */
  p4est_quadrant_t  **test;
  int                 level;
  int                *test_level;
  int                *level_num;
  int                *start_idx2;
  size_t            **index;
  size_t             *this_index;
  size_t             *first_index;
  size_t             *count;
  int                *quad_idx2;
  bool               *refine;
  int                 level_idx2;
  sc_array_t          test_view;

  /* face variables */
  int                 face[2];
#ifndef P4_TO_P8
  int                 rface[2];
#endif
  int                 num_to_child[P4EST_CHILDREN];
  p4est_iter_face_args_t face_args;
  int                 orientation;
  p4est_iter_face_info_t face_info;
  p4est_iter_face_side_t face_sides[2], *fside;

  /* edge variables */
#ifdef P4_TO_P8
  p8est_iter_edge_args_t edge_args;
  sc_array_t          common_corners[2];
  int                *e_common_corners[2], *e_edges;
  int                 perm, corner0, corner1, ne, m, l, dir;
  p4est_topidx_t     *e_treeids;
  p8est_iter_edge_info_t edge_info;
  p8est_iter_edge_side_t *edge_sides, *eside;
#endif

  /* corner variables */
  int                 nc;
  int                *c_corners;
  p4est_topidx_t     *c_treeids;
  p4est_iter_corner_args_t corner_args;
  p4est_iter_corner_info_t corner_info;
  p4est_iter_corner_side_t *corner_sides, *cside;

  /* volume info */
  p4est_iter_volume_info_t info;

  /* tier_ring variables */
  int                 tier_ring_max;
  int                 tier_level_max;
  sc_array_t          tier_rings;
  p4est_iter_tier_ring_t *ring;
  p4est_iter_tier_t  *tier;

  /* multiuse variable */
  int                 modulus;

  P4EST_ASSERT (p4est_is_valid (p4est));
  P4EST_ASSERT (ghost_layer != NULL);
  P4EST_ASSERT (sc_array_is_sorted
                (ghost_layer, p4est_quadrant_compare_piggy));

  if (p4est->first_local_tree < 0) {
    return;
  }

  /* simple loop if there is only a volume callback */
  if (iter_face == NULL && iter_corner == NULL
#ifdef P4_TO_P8
      && iter_edge == NULL
#endif
    ) {
    p4est_volume_iterate (p4est, ghost_layer, user_data, iter_volume);
    return;
  }

  /** alloc_size is the number of index arrays that are needed in the program.
   * at minimum we need two for each side of the face iterator: one for local,
   * one for ghost */
  alloc_size = 4;
  /** in the absence of strange corners (or strange edges), P4EST_CHILDREN is
   * the most quadrants that can meet at a corner */
  max_corner_size = P4EST_CHILDREN;
#ifdef P4_TO_P8
  /** if there are no strange edges between trees, then at most 4 quadrants
   * meet at an edge */
  max_edge_size = 4;
  if (iter_edge != NULL || iter_corner != NULL) {
    for (e = 0; e < num_edges; e++) {
      edge_size = (int) (ett_offset[e + 1] - ett_offset[e]);
      max_edge_size = (edge_size > max_edge_size) ? edge_size : max_edge_size;
    }
    /** we need to have two index arrays for every side of the edge iterator:
     * one for local, one for ghost */
    alloc_size = (2 * max_edge_size > alloc_size) ?
      2 * max_edge_size : alloc_size;
    /** even if there are no strange corners, for a corner that is in the
     * middle of a strange edge, there will be two quadrants that meet at the
     * corner for every quadrant that meets at the edge */
    max_corner_size = (max_edge_size * 2 > max_corner_size) ?
      max_edge_size * 2 : max_corner_size;
  }
#endif

  if (iter_corner != NULL) {
    for (c = 0; c < num_corners; c++) {
      corner_size = (int) (ctt_offset[c + 1] - ctt_offset[c]);
      max_corner_size = (corner_size > max_corner_size) ? corner_size :
        max_corner_size;
    }
    /** Similar to edges, we need to arrays for every quadrant that meets at a
     * corner */
    alloc_size = (2 * max_corner_size > alloc_size) ?
      2 * max_corner_size : alloc_size;
  }

  /** initialize arrays that keep track of where we are in the search */
  index = P4EST_ALLOC (size_t *, alloc_size);
  quadrants = P4EST_ALLOC (sc_array_t *, alloc_size);
  start_idx2 = P4EST_ALLOC (int, alloc_size / 2);
  first_index = P4EST_ALLOC (size_t, alloc_size);
  test = P4EST_ALLOC (p4est_quadrant_t *, alloc_size);
  count = P4EST_ALLOC (size_t, alloc_size);
  test_level = P4EST_ALLOC (int, alloc_size);
  quad_idx2 = P4EST_ALLOC (int, alloc_size / 2);
  refine = P4EST_ALLOC (bool, alloc_size / 2);
  for (i = 0; i < alloc_size; i++) {
    index[i] = P4EST_ALLOC (size_t, P4EST_MAXLEVEL * P4EST_ITER_STRIDE);
    if (i % 2)
      quadrants[i] = ghost_layer;
  }
  level_num = P4EST_ALLOC (int, P4EST_MAXLEVEL);

  /* initialize the tier rings */
  tier_ring_max = (p4est->mpisize == 1 ? P4EST_CHILDREN : 2 * P4EST_CHILDREN);
  tier_level_max = P4EST_QMAXLEVEL - 1;
  sc_array_init (&tier_rings, sizeof (p4est_iter_tier_ring_t));
  sc_array_resize (&tier_rings, (size_t) tier_level_max);
  for (i = 0; i < tier_level_max; i++) {
    ring = sc_array_index_int (&tier_rings, i);
    ring->next = 0;
    sc_array_init (&(ring->tiers), sizeof (p4est_iter_tier_t));
    sc_array_resize (&(ring->tiers), (size_t) tier_ring_max);
    for (j = 0; j < tier_ring_max; j++) {
      tier = sc_array_index_int (&(ring->tiers), j);
      tier->key = NULL;
    }
  }

  /** Divide the ghost_layer by p.which_tree */
  global_num_trees = trees->elem_count;
  tree_first_ghost = P4EST_ALLOC (size_t, global_num_trees + 1);
  p4est_split_ghost_layer_by_tree (ghost_layer, tree_first_ghost,
                                   global_num_trees);

  /** set up the arguments passed to face_iterate that are invariant */
  face_args.level_num = level_num;
  face_args.index = index;
  face_args.quadrants = quadrants;
  face_args.num_to_child = num_to_child;
  face_args.start_idx2 = start_idx2;
  face_args.first_index = first_index;
  face_args.test = test;
  face_args.count = count;
  face_args.test_level = test_level;
  face_args.quad_idx2 = quad_idx2;
  face_args.refine = refine;
#ifdef P4_TO_P8
  face_args.edge_args = &edge_args;
#endif
  face_args.corner_args = &corner_args;
  face_args.tier_rings = &tier_rings;
  face_info.p4est = p4est;
  face_info.ghost_layer = ghost_layer;
  sc_array_init_data (&face_info.sides, face_sides, sizeof (*face_sides), 2);
  face_args.info = &face_info;

  /** set up the arguments passed to corner_iterate that are invariant */
  if (iter_corner != NULL) {
    c_corners = P4EST_ALLOC (int, max_corner_size);
    c_treeids = P4EST_ALLOC (p4est_topidx_t, max_corner_size);
    corner_sides = P4EST_ALLOC (p4est_iter_corner_side_t, max_corner_size);
    corner_args.quadrants = quadrants;
    corner_args.index = index;
    corner_args.start_idx2 = start_idx2;
    corner_args.first_index = first_index;
    corner_args.test = test;
    corner_args.count = count;
    corner_args.quad_idx2 = quad_idx2;
    sc_array_init_data (&corner_info.sides, corner_sides,
                        sizeof (*corner_sides), (size_t) max_corner_size);
    corner_info.p4est = p4est;
    corner_info.ghost_layer = ghost_layer;
    corner_args.info = &corner_info;
  }

#ifdef P4_TO_P8
  /** set up the arguments passed to edge_iterate that are invariant */
  if (iter_edge != NULL || iter_corner != NULL) {
    e_common_corners[0] = P4EST_ALLOC (int, max_edge_size);
    e_common_corners[1] = P4EST_ALLOC (int, max_edge_size);
    e_edges = P4EST_ALLOC (int, max_edge_size);
    e_treeids = P4EST_ALLOC (p4est_topidx_t, max_edge_size);
    edge_sides = P4EST_ALLOC (p8est_iter_edge_side_t, max_edge_size);
    sc_array_init_data (&(common_corners[0]), e_common_corners[0],
                        sizeof (*(e_common_corners[0])),
                        (size_t) max_edge_size);
    sc_array_init_data (&(common_corners[1]), e_common_corners[1],
                        sizeof (*(e_common_corners[1])),
                        (size_t) max_edge_size);
    edge_args.level_num = level_num;
    edge_args.quadrants = quadrants;
    edge_args.index = index;
    edge_args.start_idx2 = start_idx2;
    edge_args.first_index = first_index;
    edge_args.test = test;
    edge_args.count = count;
    edge_args.test_level = test_level;
    edge_args.quad_idx2 = quad_idx2;
    edge_args.refine = refine;
    edge_args.common_corners[0] = &(common_corners[0]);
    edge_args.common_corners[1] = &(common_corners[1]);
    edge_args.corner_args = &corner_args;
    edge_args.tier_rings = &tier_rings;
    sc_array_init_data (&edge_info.sides, edge_sides, sizeof (*edge_sides),
                        (size_t) max_edge_size);
    edge_info.p4est = p4est;
    edge_info.ghost_layer = ghost_layer;
    edge_args.info = &edge_info;
  }
#endif

  /** set up the arguments passed to iter_volume that are invariant */
  info.p4est = p4est;
  info.ghost_layer = ghost_layer;

  /** we have to loop over all trees and not just local trees because of the
   * ghost layer */
  for (t = 0; t < global_num_trees; t++) {

    tree = p4est_array_index_topidx (trees, t);
    quadrants[local] = &(tree->quadrants);

    P4EST_ASSERT (p4est_tree_is_sorted (tree));
    P4EST_ASSERT (tree->maxlevel <= P4EST_QMAXLEVEL);

    /* set up the index information for level 0 */
    this_index = index[left * 2 + local];
    this_index[0] = 0;
    num_local = quadrants[local]->elem_count;
    this_index[1] = (p4est_locidx_t) num_local;
    this_index = index[left * 2 + ghost];
    this_index[0] = tree_first_ghost[t];
    this_index[1] = tree_first_ghost[t + 1];
    num_ghosts = (this_index[1] - this_index[0]);

    if (!num_local && !num_ghosts) {
      continue;
    }

    /* since we are intra-tree, both treeids are the current tree,
     * orientation is always 0 */
    face_args.outside_face = false;
    face_info.orientation = 0;
    sc_array_resize (&face_info.sides, 2);
    for (side = left; side <= right; side++) {
      face_args.quadrants[side * 2 + local] = quadrants[local];
      fside = sc_array_index_int (&face_info.sides, side);
      fside->treeid = t;
    }

#ifdef P4_TO_P8
    /* since we are intra-tree, all treeids and quadrant arrays correspond to
     * the current tree */
    if (iter_edge != NULL || iter_corner != NULL) {
      edge_args.num_sides = 4;
      sc_array_resize (&(common_corners[0]), 4);
      sc_array_resize (&(common_corners[1]), 4);
      sc_array_resize (&edge_info.sides, 4);
      for (side = 0; side < 4; side++) {
        e_treeids[side] = t;
        edge_args.quadrants[side * 2 + local] = quadrants[local];
        eside = sc_array_index_int (&edge_info.sides, side);
        eside->treeid = t;
        eside->orientation = 0;
      }
    }
#endif

    /* since we are intra-tree, all treeids and quadrant arrays correspond to
     * the current tree */
    if (iter_corner != NULL) {
      corner_args.num_sides = P4EST_CHILDREN;
      sc_array_resize (&corner_info.sides, P4EST_CHILDREN);
      for (i = 0; i < P4EST_CHILDREN; i++) {
        c_treeids[i] = t;
        quadrants[2 * i] = quadrants[local];
        cside = sc_array_index_int (&corner_info.sides, i);
        cside->treeid = t;
      }
    }

    info.treeid = t;

    level = 0;
    level_idx2 = 0;
    level_num[0] = 0;

    /* there is no need to iterate within the tree if there are no local
     * quadrants */
    if (num_local) {

      /* if we try to advance the branch number at level 0, then we have
       * completed the search through the tree */
      while (level_num[0] == 0) {
        /* if we have advanced the branch number to P4EST_CHILDREN, we have
         * completed the search at this level, and we now set up the face_args
         * for the faces between the quad branches at this level */
        if (level_num[level] == P4EST_CHILDREN) {
          /* for each direction */
          for (left_side = 0, right_side = 1; right_side < 2 * P4EST_DIM;
               left_side += 2, right_side += 2) {
            fside = sc_array_index (&face_info.sides, 0);
            fside->face = right_side;
            fside = sc_array_index (&face_info.sides, 1);
            fside->face = left_side;
            /* set up the num_to_child array
             * (see face_args_t for a description) */
            for (i = 0; i < P4EST_CHILDREN / 2; i++) {
#ifndef P4_TO_P8
              num_to_child[left * ntc_str + i] = p4est_face_corners
                [p4est_zface_to_rface[right_side]][i];
              num_to_child[right * ntc_str + i] = p4est_face_corners
                [p4est_zface_to_rface[left_side]][i];
#else
              num_to_child[left * ntc_str + i] =
                p8est_face_corners[right_side][i];
              num_to_child[right * ntc_str + i] =
                p8est_face_corners[left_side][i];
#endif
            }
            face_args.level = level;
            for (i = 0; i < P4EST_CHILDREN / 2; i++) {
#ifndef P4_TO_P8
              start_idx2[left] =
                p4est_face_corners[p4est_zface_to_rface[left_side]][i];
              start_idx2[right] =
                p4est_face_corners[p4est_zface_to_rface[right_side]][i];
#else
              start_idx2[left] = p8est_face_corners[left_side][i];
              start_idx2[right] = p8est_face_corners[right_side][i];
#endif
              /** initialize index[right] for use in face_iterate */
              quad_idx2[0] = level_idx2 + start_idx2[right];
              index[right * 2 + local][quad_idx2[0]] =
                index[left * 2 + local][quad_idx2[0]];
              index[right * 2 + local][quad_idx2[0] + 1] =
                index[left * 2 + local][quad_idx2[0] + 1];
              index[right * 2 + ghost][quad_idx2[0]] =
                index[left * 2 + ghost][quad_idx2[0]];
              index[right * 2 + ghost][quad_idx2[0] + 1] =
                index[left * 2 + ghost][quad_idx2[0] + 1];

              p4est_face_iterate (&face_args, user_data, iter_face,
#ifdef P4_TO_P8
                                  iter_edge,
#endif
                                  iter_corner);
            }
          }
#ifdef P4_TO_P8
          /* if there is an edge or a corner callback, we need to use
           * edge_iterate, so we set up the common corners and edge ids
           * for all of the edges between the quad search areas */
          if (iter_edge != NULL || iter_corner != NULL) {
            edge_args.level = level;
            for (dir = 0; dir < P4EST_DIM; dir++) {
              for (side = 0; side < 2; side++) {
                for (j = 0; j < 4; j++) {
                  e_common_corners[side][j] =
                    p8est_face_corners[dir * 2 + side][3 - j];
                }
              }
              for (j = 0; j < 4; j++) {
                eside = sc_array_index_int (&edge_info.sides, j);
                eside->edge = 4 * dir + (3 - j);
              }
              for (side = 0; side < 2; side++) {
                for (j = 0; j < 4; j++) {
                  start_idx2[j] = p8est_face_corners[dir * 2 + side][j];
                  if (j != 0) {
                    quad_idx2[0] = level_idx2 + start_idx2[j];
                    index[j * 2 + local][quad_idx2[0]] =
                      index[local][quad_idx2[0]];
                    index[j * 2 + local][quad_idx2[0] + 1] =
                      index[local][quad_idx2[0] + 1];
                    index[j * 2 + ghost][quad_idx2[0]] =
                      index[ghost][quad_idx2[0]];
                    index[j * 2 + ghost][quad_idx2[0] + 1] =
                      index[ghost][quad_idx2[0] + 1];
                  }
                }
                p8est_edge_iterate (&edge_args, user_data, iter_edge,
                                    iter_corner);
              }
            }
          }
#endif
          /* if there is a corner callback, we need to call corner_iterate on
           * the corner in the middle of the quad search areas */
          if (iter_corner != NULL) {
            corner_args.level = level;
            for (i = 0; i < P4EST_CHILDREN; i++) {
              cside = sc_array_index_int (&corner_info.sides, i);
              cside->corner = P4EST_CHILDREN - i - 1;
              start_idx2[i] = i;
              quad_idx2[0] = level_idx2 + start_idx2[i];
              index[i * 2 + local][quad_idx2[0]] = index[local][quad_idx2[0]];
              index[i * 2 + local][quad_idx2[0] + 1] =
                index[local][quad_idx2[0] + 1];
              index[i * 2 + ghost][quad_idx2[0]] = index[ghost][quad_idx2[0]];
              index[i * 2 + ghost][quad_idx2[0] + 1] =
                index[ghost][quad_idx2[0] + 1];
            }
            p4est_corner_iterate (&corner_args, user_data, iter_corner);
          }
          /* we are done at the level, so we go up a level and over a branch */
          level_num[--level]++;
          level_idx2 -= P4EST_ITER_STRIDE;
          continue;
        }

        /* quad_idx[0] now gives the location in index[type] of the bounds
         * of the current search area, from which we get the first quad
         * and the count */
        quad_idx2[0] = level_idx2 + level_num[level];
        for (type = local; type <= ghost; type++) {
          this_index = index[type];
          first_index[type] = this_index[quad_idx2[0]];
          count[type] = (this_index[quad_idx2[0] + 1] - first_index[type]);
        }
        /* if there are no local quadrants, we are done with this search area,
         * and we advance to the next branch at this level */
        if (!count[local]) {
          level_num[level]++;
          continue;
        }

        for (type = local; type <= ghost; type++) {
          if (count[type]) {
            P4EST_ASSERT (first_index[type] < quadrants[type]->elem_count);
            test[type] = sc_array_index (quadrants[type], first_index[type]);
            test_level[type] = (int) test[type]->level;
            /* if the test quad form this type is the size of the search area,
             * we run the quad callback */
            if (test_level[type] == level) {
              if (type == local) {
                info.quad = test[type];
                info.quadid = (p4est_locidx_t) first_index[type];
                if (iter_volume != NULL)
                  iter_volume (&info, user_data);
              }
              break;
            }
          }
          else {
            test[type] = NULL;
            test_level[type] = -1;
          }
        }
        /* if we ran a quad callback, then we are done with this search
         * area, and advance to the next branch on this level */
        if (type <= ghost) {
          level_num[level]++;
          continue;
        }

        /* otherwise, we need to refine the search areas, and place the
         * bounds of the refined search areas on the next tier in index[type]*/
        quad_idx2[0] = level_idx2 + P4EST_ITER_STRIDE;
        for (type = local; type <= ghost; type++) {
          sc_array_init_view (&test_view, quadrants[type],
                              first_index[type], count[type]);
          this_index = index[type];
          p4est_iter_tier_insert (&test_view, level,
                                  this_index + quad_idx2[0],
                                  first_index[type], &tier_rings, test[type]);
        }

        /* We descend to the next level and continue the search there */
        level_num[++level] = 0;
        level_idx2 += P4EST_ITER_STRIDE;
      }
    }

    /* Now we need to run face_iterate on the faces between trees */
    for (i = 0; i < 2 * P4EST_DIM; i++) {
      /* the level of the initial search area has to be the whole tree, i.e.
       * level 0 */
      face_args.level = 0;
      face_args.start_idx2[left] = face_args.start_idx2[right] = 0;
      sc_array_resize (&face_info.sides, 2);

      /* we find the corresponding tree nt and face[right]  across from
       * face[left] = i in tree t, and we initialize num_to_child
       * to agree with this face/orientation combo */
#ifndef P4_TO_P8
      rface[left] = p4est_zface_to_rface[i];
#endif
      for (j = 0; j < P4EST_CHILDREN / 2; j++) {
#ifndef P4_TO_P8
        num_to_child[left * ntc_str + j] = p4est_face_corners[rface[left]][j];
#else
        num_to_child[left * ntc_str + j] = p4est_face_corners[i][j];
#endif
      }
      face[left] = i;
      fside = sc_array_index (&face_info.sides, 0);
      fside->treeid = t;
      fside->face = i;
#ifndef P4_TO_P8
      nt = ttt[t * 2 * P4EST_DIM + rface[left]];
      rface[right] = (int) ttf[t * 2 * P4EST_DIM + rface[left]];
      face_info.orientation = rface[right] / 4;
      rface[right] %= 4;
      P4EST_ASSERT (ttt[nt * 2 * P4EST_DIM + rface[right]] == t);
      face[right] = p4est_rface_to_zface[rface[right]];
      fside = sc_array_index (&face_info.sides, 1);
      fside->treeid = nt;
      fside->face = p4est_rface_to_zface[rface[right]];
#else
      nt = ttt[t * 2 * P4EST_DIM + face[left]];
      face[right] = (int) ttf[t * 2 * P4EST_DIM + face[left]];
      face[right] %= 6;
      fside = sc_array_index (&face_info.sides, 1);
      fside->treeid = nt;
      fside->face = (int) ttf[t * 2 * P4EST_DIM + face[left]];
      face_info.orientation = fside->face / 6;
      fside->face %= 6;
#endif

      /* we only want to call face_iterate once per face */
      if ((t > nt) || (t == nt && face[left] > face[right])) {
        face_args.outside_face = false;
#ifndef P4_TO_P8
        modulus = face[left] ^ face[right];
        modulus = ((modulus & 2) >> 1) ^ (modulus & 1) ^ 1;
        if ((modulus ^ face_info.orientation) == 0) {
          num_to_child[right * ntc_str] = p4est_face_corners[rface[right]][0];
          num_to_child[right * ntc_str + 1] =
            p4est_face_corners[rface[right]][1];
        }
        else {
          num_to_child[right * ntc_str] = p4est_face_corners[rface[right]][1];
          num_to_child[right * ntc_str + 1] =
            p4est_face_corners[rface[right]][0];
        }
#else
        modulus = p8est_face_permutation_refs[face[left]][face[right]];
        modulus = p8est_face_permutation_sets[modulus][face_info.orientation];
        for (j = 0; j < P4EST_CHILDREN / 2; j++) {
          k = p8est_face_permutations[modulus][j];
          num_to_child[right * ntc_str + j] =
            p8est_face_corners[face[right]][k];
        }
#endif

        /* we need to set up the quadrants arrays and the index arrays for the 
         * right side of the face */
        level_num[0] = 0;
        P4EST_ASSERT (nt < trees->elem_count);
        tree = p4est_array_index_topidx (trees, nt);
        face_args.quadrants[left * 2 + local] = quadrants[local];
        face_args.quadrants[right * 2 + local] = &(tree->quadrants);
        index[right * 2 + local][0] = 0;
        index[right * 2 + local][1] = num_local = tree->quadrants.elem_count;
        index[right * 2 + ghost][0] = tree_first_ghost[nt];
        index[right * 2 + ghost][1] = tree_first_ghost[nt + 1];

#ifdef P4_TO_P8
        /* we set the arrays in edge_args to the appropriate size and
         * make sure that the quadrants arrays are pointed at the appropriate
         * arrays, and that the tree ids reflect which side of the
         * face each edge-side is on*/
        if (iter_edge != NULL || iter_corner != NULL) {
          edge_args.num_sides = 4;
          sc_array_resize (&(common_corners[0]), 4);
          sc_array_resize (&(common_corners[1]), 4);
          sc_array_resize (&edge_info.sides, 4);
          for (j = 0; j < 4; j++) {
            edge_args.quadrants[j * 2 + local] =
              (j % 2) ? face_args.quadrants[right * 2 + local] :
              face_args.quadrants[left * 2 + local];
            eside = sc_array_index_int (&edge_info.sides, j);
            eside->treeid = (j % 2) ? nt : t;
          }
        }
#endif

        /* we set the arrays in edge_args the appropriate size and
         * make sure that the quadrants arrays are pointed at the appropriate
         * arrays, and that the tree ids reflect which side of the
         * face each corner-side is on*/
        if (iter_corner != NULL) {
          corner_args.num_sides = P4EST_CHILDREN;
          sc_array_resize (&corner_info.sides, P4EST_CHILDREN);
          for (j = 0; j < P4EST_CHILDREN; j++) {
            cside = sc_array_index_int (&corner_info.sides, j);
            if (j % 2) {
              cside->treeid = nt;
              corner_args.quadrants[j * 2 + local] =
                face_args.quadrants[right * 2 + local];
            }
            else {
              cside->treeid = t;
              corner_args.quadrants[j * 2 + local] =
                face_args.quadrants[left * 2 + local];
            }
          }
        }

        p4est_face_iterate (&face_args, user_data, iter_face,
#ifdef P4_TO_P8
                            iter_edge,
#endif
                            iter_corner);
      }
      /* if we are an outside face, the edge/corner args need to be 
       * initialized differently */
      else if (t == nt && face[left] == face[right]) {
        sc_array_resize (&face_info.sides, 1);
        face_args.outside_face = true;
        level_num[0] = 0;
#ifdef P4_TO_P8
        if (iter_edge != NULL || iter_corner != NULL) {
          edge_args.num_sides = 2;
          sc_array_resize (&(common_corners[0]), 2);
          sc_array_resize (&(common_corners[1]), 2);
          sc_array_resize (&edge_info.sides, 2);
          for (j = 0; j < 2; j++) {
            eside = sc_array_index_int (&edge_info.sides, j);
            eside->treeid = t;
            edge_args.quadrants[j * 2 + local] =
              face_args.quadrants[left * 2 + local];
          }
        }
#endif

        if (iter_corner != NULL) {
          corner_args.num_sides = P4EST_CHILDREN / 2;
          sc_array_resize (&corner_info.sides, P4EST_CHILDREN / 2);
          for (j = 0; j < P4EST_CHILDREN / 2; j++) {
            cside = sc_array_index_int (&corner_info.sides, j);
            cside->treeid = t;
            quadrants[2 * j] = quadrants[local];
          }
        }

        p4est_face_iterate (&face_args, user_data, iter_face,
#ifdef P4_TO_P8
                            iter_edge,
#endif
                            iter_corner);
      }
    }

    /* if there is an edge or a corner callback, we need to run
     * edge_iterate on the edges between trees */
#ifdef P4_TO_P8
    if (iter_edge != NULL || iter_corner != NULL) {
      for (i = 0; i < 12; i++) {
        edge_args.level = 0;
        e_edges[0] = i;
        /* find out if there is a global edge corresponding to this edge */
        if (num_edges > 0) {
          this_edge = tte[t * 12 + i];
        }
        else {
          this_edge = -1;
        }
        modulus = 0;

        /* if there is a global edge, modulus = 0 if they have the same
         * orientation, = 1 if opposite */
        if (this_edge > -1) {
          for (jt = ett_offset[this_edge];
               jt < ett_offset[this_edge + 1]; jt++) {
            if (edge_to_tree[jt] == t) {
              modulus = (int) edge_to_edge[jt];
              if (modulus % 12 == i) {
                modulus /= 12;
                break;
              }
            }
          }
        }

        /* set the treeid and common corners that go with edge */
        edge_args.num_sides = 1;
        e_treeids[0] = t;
        if (modulus == 0) {
          e_common_corners[0][0] = p8est_edge_corners[i][0];
          e_common_corners[1][0] = p8est_edge_corners[i][1];
        }
        else {
          e_common_corners[0][0] = p8est_edge_corners[i][1];
          e_common_corners[1][0] = p8est_edge_corners[i][0];
        }

        /* for each face touching the edge, add the tree across that face,
         * since it must also touch the edge */
        for (side = left; side <= right; side++) {
          face[0] = p8est_edge_faces[i][side];
          nt = ttt[t * 6 + face[0]];
          face[1] = (int) ttf[t * 6 + face[0]];
          orientation = face[1] / 6;
          face[1] %= 6;
          /* figure out which common corners from the neighbor tree
           * correspond to the common corners from tree t*/
          if (nt != t || face[1] != face[0]) {
            perm = p8est_face_permutation_refs[face[0]][face[1]];
            perm = p8est_face_permutation_sets[perm][orientation];
            corner0 = p8est_face_permutations[perm]
              [p8est_edge_face_corners[i][face[0]][0]];
            corner1 = p8est_face_permutations[perm]
              [p8est_edge_face_corners[i][face[0]][1]];
            k = edge_args.num_sides;
            if (modulus == 0) {
              e_common_corners[0][k] = p8est_face_corners[face[1]][corner0];
              e_common_corners[1][k] = p8est_face_corners[face[1]][corner1];
            }
            else {
              e_common_corners[0][k] = p8est_face_corners[face[1]][corner1];
              e_common_corners[1][k] = p8est_face_corners[face[1]][corner0];
            }

            /* figure out the edge between the common corners */
            ne = e_common_corners[1][k] - e_common_corners[0][k];
            ne = (ne > 0) ? ne : -ne;
            ne = p8est_corner_edges[e_common_corners[0][k]][ne >> 1];
            e_edges[k] = ne;
            e_treeids[k] = nt;
            edge_args.num_sides++;
          }
        }

        /* if there is a global edge, add all of the edge neighbors */
        if (this_edge > -1) {
          for (jt = ett_offset[this_edge];
               jt < ett_offset[this_edge + 1]; jt++) {
            nt = edge_to_tree[jt];
            ne = (int) edge_to_edge[jt];
            orientation = ne / 12;
            ne %= 12;
            if (orientation == 0) {
              corner0 = p8est_edge_corners[ne][0];
              corner1 = p8est_edge_corners[ne][1];
            }
            else {
              corner0 = p8est_edge_corners[ne][1];
              corner1 = p8est_edge_corners[ne][0];
            }
            /* only add if the neighbor is not the original and
             * is not a face neighbor */
            for (k = 0; k < edge_args.num_sides; k++) {
              if (nt == e_treeids[k]) {
                if (corner0 == e_common_corners[0][k] &&
                    corner1 == e_common_corners[1][k]) {
                  break;
                }
              }
            }
            if (k == edge_args.num_sides) {
              e_treeids[k] = nt;
              e_edges[k] = ne;
              e_common_corners[0][k] = corner0;
              e_common_corners[1][k] = corner1;
              edge_args.num_sides++;
            }
          }
        }
        /* we only want to run edge_iterate on this edge once,
         * so we only run if this is the highest tree id and
         * the highest edge on that tree in the list */
        for (j = 1; j < edge_args.num_sides; j++) {
          if (e_treeids[j] > t) {
            break;
          }
          if (e_treeids[j] == t && e_edges[j] > i) {
            break;
          }
        }
        /* set up the quadrants arrays to point to the right arrays and the
         * index arrays to the correct initial bounds */
        if (j == edge_args.num_sides) {
          sc_array_resize (&edge_info.sides, (size_t) edge_args.num_sides);
          for (j = 0; j < edge_args.num_sides; j++) {
            edge_args.start_idx2[j] = 0;
            index[j * 2 + local][0] = 0;
            tree = p4est_array_index_topidx (trees, e_treeids[j]);
            quadrants[j * 2 + local] = &(tree->quadrants);
            index[j * 2 + local][1] = quadrants[j * 2 + local]->elem_count;
            index[j * 2 + ghost][0] = tree_first_ghost[e_treeids[j]];
            index[j * 2 + ghost][1] = tree_first_ghost[e_treeids[j] + 1];
            eside = sc_array_index_int (&edge_info.sides, j);
            eside->treeid = e_treeids[j];
            eside->edge = e_edges[j];
            if (e_common_corners[0][j] < e_common_corners[1][j]) {
              eside->orientation = 0;
            }
            else {
              eside->orientation = 1;
            }
          }
          sc_array_resize (&(common_corners[0]),
                           (size_t) edge_args.num_sides);
          sc_array_resize (&(common_corners[1]),
                           (size_t) edge_args.num_sides);

          /* if there is a corner callback, initialize the quadrants arrays to
           * point to the correct arrays and make sure the arrays
           * are the correct sizes */
          if (iter_corner != NULL) {
            corner_args.num_sides = edge_args.num_sides * 2;
            sc_array_resize (&corner_info.sides,
                             (size_t) corner_args.num_sides);
            for (j = 0; j < corner_args.num_sides; j++) {
              eside = sc_array_index_int (&edge_info.sides,
                                          j % edge_args.num_sides);
              cside = sc_array_index_int (&corner_info.sides, j);
              cside->treeid = eside->treeid;
              quadrants[2 * j] = quadrants[2 * (j % edge_args.num_sides)];
            }
          }

          p8est_edge_iterate (&edge_args, user_data, iter_edge, iter_corner);
        }
      }
    }
#endif

    if (iter_corner != NULL) {
#ifndef P4_TO_P8
      sc_array_init (cta, sizeof (p4est_corner_transform_t));
#endif

      for (i = 0; i < P4EST_CHILDREN; i++) {
        corner_args.level = 0;
        c_treeids[0] = t;
        c_corners[0] = i;
        corner_args.num_sides = 1;

        /* for every face that touches corner i, add the tree neighbor across
         * that face to the search */
        for (j = 0; j < P4EST_DIM; j++) {
          face[0] = p4est_corner_faces[i][j];
          nt = ttt[t * P4EST_DIM * 2 + face[0]];
          face[1] = (int) ttf[t * P4EST_DIM * 2 + face[0]];
          orientation = face[1] / (P4EST_DIM * 2);
          face[1] = face[1] % (P4EST_DIM * 2);
          if (nt == t && face[0] == face[1]) {
            continue;
          }
          c_treeids[corner_args.num_sides] = nt;
          /* find the correct corner in nt corresponding to corner i in tree t
           */
          for (k = 0; k < P4EST_CHILDREN / 2; k++) {
            if (p4est_face_corners[face[0]][k] == i) {
#ifndef P4_TO_P8
              rface[0] = face[0];
              rface[1] = face[1];
              face[0] = p4est_rface_to_zface[rface[0]];
              face[1] = p4est_rface_to_zface[rface[1]];
              modulus = face[0] ^ face[1];
              modulus = ((modulus & 2) >> 1) ^ (modulus & 1) ^ 1;
              face[0] = rface[0];
              face[1] = rface[1];
              c_corners[corner_args.num_sides] =
                p4est_face_corners[face[1]][modulus ^ orientation ^ k];
#else
              modulus = p8est_face_permutation_refs[face[0]][face[1]];
              modulus = p8est_face_permutation_sets[modulus][orientation];
              c_corners[corner_args.num_sides] = p8est_face_corners[face[1]]
                [p8est_face_permutations[modulus][k]];
#endif
            }
          }
          corner_args.num_sides++;
        }

#ifdef P4_TO_P8
        /* for every edge that touches corner i, add the tree neighbors
         * across this edge to the search */
        if (num_edges > 0) {
          for (j = 0; j < P4EST_DIM; j++) {
            m = p8est_corner_edges[i][j];
            this_edge = tte[t * 12 + m];
            if (this_edge == -1) {
              continue;
            }
            modulus = 0;
            /* figure out which corner of the global edge corresponds
             * to corner i on tree t */
            for (kt = ett_offset[this_edge];
                 kt < ett_offset[this_edge + 1]; kt++) {
              if (edge_to_tree[kt] == t) {
                modulus = (int) edge_to_edge[kt];
                if (modulus % 12 == m) {
                  modulus /= 12;
                  break;
                }
              }
            }
            if (p8est_edge_corners[m][modulus] == i) {
              modulus = 0;
            }
            else {
              modulus = 1;
            }

            /* find the correct corner in nt corresponding to corner i in tree
             * t */
            for (kt = ett_offset[this_edge];
                 kt < ett_offset[this_edge + 1]; kt++) {
              nt = edge_to_tree[kt];
              ne = (int) edge_to_edge[kt];
              orientation = ne / 12;
              ne %= 12;
              nc = p8est_edge_corners[ne][orientation ^ modulus];
              /* only add tree/corners that haven't been added already */
              for (l = 0; l < corner_args.num_sides; l++) {
                if (nt == c_treeids[l] && nc == c_corners[l]) {
                  break;
                }
              }
              if (l == corner_args.num_sides) {
                c_treeids[l] = nt;
                c_corners[l] = nc;
                corner_args.num_sides++;
              }
            }
          }
        }
#endif

        if (num_corners > 0) {
#ifndef P4_TO_P8
          this_corner = ttc[t * P4EST_CHILDREN + p4est_corner_to_zorder[i]];
#else
          this_corner = ttc[t * P4EST_CHILDREN + i];
#endif
        }
        else {
          this_corner = -1;
        }

        /* if this corner is a global corner, add its corner neighbors */
        if (this_corner != -1) {
#ifndef P4_TO_P8
          p4est_find_corner_transform (conn, t, p4est_corner_to_zorder[i],
                                       cta);
          for (j = 0; j < (int) cta->elem_count; j++) {
            ct = sc_array_index_int (cta, j);
            nt = ct->ntree;
            nc = (int) ct->ncorner;
            /* only add tree/corners that haven't been added yet */
            for (k = 0; k < corner_args.num_sides; k++) {
              if (nt == c_treeids[k] && nc == c_corners[k]) {
                break;
              }
            }
            if (k == corner_args.num_sides) {
              c_treeids[k] = nt;
              c_corners[k] = nc;
              corner_args.num_sides++;
            }
          }
#else
          for (jt = ctt_offset[this_corner];
               jt < ctt_offset[this_corner + 1]; jt++) {
            nt = corner_to_tree[jt];
            nc = (int) corner_to_corner[jt];
            /* only add tree/corners that haven't been added yet */
            for (k = 0; k < corner_args.num_sides; k++) {
              if (nt == c_treeids[k] && nc == c_corners[k]) {
                break;
              }
            }
            if (k == corner_args.num_sides) {
              c_treeids[k] = nt;
              c_corners[k] = nc;
              corner_args.num_sides++;
            }
          }
#endif
        }

        /* we only want to run corner_iterate once, so we only run if
         * t is the greatest tree and i is the greatest corner on that tree */
        for (j = 0; j < corner_args.num_sides; j++) {
          if (c_treeids[j] > t) {
            break;
          }
          if (c_treeids[j] == t && c_corners[j] > i) {
            break;
          }
        }
        /* set up the quadrants arrays to point to the correct arrays,
         * and resize all of the corner_args arrays to the correct size */
        if (j == corner_args.num_sides) {
          sc_array_resize (&corner_info.sides,
                           (size_t) corner_args.num_sides);
          for (j = 0; j < corner_args.num_sides; j++) {
            cside = sc_array_index_int (&corner_info.sides, j);
            cside->treeid = c_treeids[j];
            cside->corner = c_corners[j];
            corner_args.start_idx2[j] = 0;
            index[j * 2 + local][0] = 0;
            tree = p4est_array_index_topidx (trees, c_treeids[j]);
            quadrants[j * 2 + local] = &(tree->quadrants);
            index[j * 2 + local][1] = quadrants[j * 2 + local]->elem_count;
            index[j * 2 + ghost][0] = tree_first_ghost[c_treeids[j]];
            index[j * 2 + ghost][1] = tree_first_ghost[c_treeids[j] + 1];
          }

          p4est_corner_iterate (&corner_args, user_data, iter_corner);
        }
      }
#ifndef P4_TO_P8
      sc_array_reset (cta);
#endif
    }

  }

  /* clean up */

  if (iter_corner != NULL) {
    P4EST_FREE (c_corners);
    P4EST_FREE (c_treeids);
    P4EST_FREE (corner_sides);
  }

#ifdef P4_TO_P8
  if (iter_edge != NULL || iter_corner != NULL) {
    P4EST_FREE (e_common_corners[0]);
    P4EST_FREE (e_common_corners[1]);
    P4EST_FREE (e_treeids);
    P4EST_FREE (e_edges);
    P4EST_FREE (edge_sides);
  }
#endif

  P4EST_FREE (tree_first_ghost);
  P4EST_FREE (level_num);
  for (i = 0; i < alloc_size; i++) {
    P4EST_FREE (index[i]);
  }
  P4EST_FREE (index);
  P4EST_FREE (refine);
  P4EST_FREE (quad_idx2);
  P4EST_FREE (test_level);
  P4EST_FREE (count);
  P4EST_FREE (test);
  P4EST_FREE (first_index);
  P4EST_FREE (start_idx2);
  P4EST_FREE (quadrants);

  for (i = 0; i < tier_level_max; i++) {
    ring = sc_array_index_int (&tier_rings, i);
    sc_array_reset (&(ring->tiers));
  }
  sc_array_reset (&tier_rings);
}
