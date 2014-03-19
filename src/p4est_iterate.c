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
#include <p8est_bits.h>
#include <p8est_iterate.h>
#include <p8est_search.h>
#else
#include <p4est_algorithms.h>
#include <p4est_bits.h>
#include <p4est_iterate.h>
#include <p4est_search.h>
#endif

/* tier ring functions:
 *
 * A tier saves the indices created  by an array split.  In all of the iterate
 * loops, those indices are uniquely determined by the level of the split and
 * the first quadrant in the section that's split.  For each level, there is a
 * set number of tiers saved in a ring.  To get the indices from a split,
 * we search the ring of the level of the split for the same first quadrant: if
 * it is found, the indices are copied from the saved tier; otherwise, they are
 * computed, and saved in the ring in place of the oldest tier in the ring */
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

static sc_array_t  *
p4est_iter_tier_rings_new (int num_procs)
{
  int                 i, j;
  int                 tier_ring_max;
  int                 tier_level_max;
  sc_array_t         *tier_rings;
  p4est_iter_tier_ring_t *ring;
  p4est_iter_tier_t  *tier;

  tier_rings = sc_array_new (sizeof (p4est_iter_tier_ring_t));
  tier_ring_max = (num_procs == 1 ? P4EST_CHILDREN : 2 * P4EST_CHILDREN);
  tier_level_max = P4EST_QMAXLEVEL;
  sc_array_resize (tier_rings, (size_t) tier_level_max);
  for (i = 0; i < tier_level_max; i++) {
    ring = (p4est_iter_tier_ring_t *) sc_array_index_int (tier_rings, i);
    ring->next = 0;
    sc_array_init (&(ring->tiers), sizeof (p4est_iter_tier_t));
    sc_array_resize (&(ring->tiers), (size_t) tier_ring_max);
    for (j = 0; j < tier_ring_max; j++) {
      tier = (p4est_iter_tier_t *) sc_array_index_int (&(ring->tiers), j);
      tier->key = NULL;
    }
  }

  return tier_rings;
}

static void
p4est_iter_tier_rings_destroy (sc_array_t * tier_rings)
{
  size_t              zz;
  p4est_iter_tier_ring_t *ring;

  for (zz = 0; zz < tier_rings->elem_count; zz++) {
    ring = (p4est_iter_tier_ring_t *) sc_array_index (tier_rings, zz);
    sc_array_reset (&(ring->tiers));
  }
  sc_array_destroy (tier_rings);
}

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
  ring = (p4est_iter_tier_ring_t *) sc_array_index_int (tier_rings, level);

  limit = (int) ring->tiers.elem_count;
  for (i = 0; i < limit; i++) {
    tier = (p4est_iter_tier_t *) sc_array_index_int (&(ring->tiers), i);
    key = tier->key;
    if (key == NULL) {
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
  tier = (p4est_iter_tier_t *) sc_array_index_int (&(ring->tiers),
                                                   ring->next++);
  memcpy (tier->array, next_tier, P4EST_ITER_STRIDE * sizeof (size_t));
  tier->key = q;
  ring->next %= limit;
}

/* loop arg functions */
typedef struct p4est_iter_loop_args
{
  int                 alloc_size;       /* large enough to accomodate strange
                                           corners/edges between trees */
#ifdef P4_TO_P8
  int8_t              loop_edge;        /* should edge_iterate be run */
#endif
  int8_t              loop_corner;      /* should corner_iterate by run */

  int                 level;
  int                *level_num;        /* an array that keeps track of which
                                           branch we take at each step in the
                                           heirarchical search areas */
  int                *quad_idx2;        /* an indexing variable used in
                                           the iterate functions: passed as an
                                           argument to avoid using alloc/free
                                           on each call */
  sc_array_t        **quadrants;        /* the arrays, two for each side (one
                                           local, one ghost), that contain the
                                           quadrants in each search area */
  size_t            **index;    /* for each sidetype, the indices in quadrants
                                   that form the bounds of the heirarchical
                                   search areas */
  size_t             *first_index;      /* an indexing variable used in the
                                           iterate functions: passed as an
                                           argument to avoid using alloc/free
                                           on each call */
  size_t             *count;    /* a counting variable used in the iterate
                                   functions: passed as an argument to
                                   avoid using alloc/free on each call */
  p4est_quadrant_t  **test;     /* a testing variable used in the iterate
                                   functions: passed as an argument to
                                   avoid using alloc/free on each call */
  int                *test_level;       /* a testing variable used in
                                           the iterate functions:: passed as an
                                           argument to avoid using alloc/free
                                           on each call */
  int8_t             *refine;   /* a testing variable used in the iterate
                                   functions: passed as an argument to avoid
                                   using alloc/free on each call */
  sc_array_t         *tier_rings;
}
p4est_iter_loop_args_t;

static p4est_iter_loop_args_t *
p4est_iter_loop_args_new (p4est_connectivity_t * conn,
#ifdef P4_TO_P8
                          p8est_iter_edge_t iter_edge,
#endif
                          p4est_iter_corner_t iter_corner,
                          p4est_ghost_t * ghost_layer, int num_procs)
{
  int                 i;
  p4est_topidx_t      c;
  int                 alloc_size;
  p4est_topidx_t      num_corners = conn->num_corners;
  const p4est_topidx_t *ctt_offset = conn->ctt_offset;
#ifdef P4_TO_P8
  p4est_topidx_t      e;
  int                 max_edge_size;
  int                 edge_size;
  const p4est_topidx_t *ett_offset = conn->ett_offset;
  p4est_topidx_t      num_edges = conn->num_edges;
#endif
  int                 max_corner_size;
  int                 corner_size;
  p4est_iter_loop_args_t *loop_args;

  loop_args = P4EST_ALLOC (p4est_iter_loop_args_t, 1);

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
  loop_args->alloc_size = alloc_size;
  loop_args->level_num = P4EST_ALLOC (int, (P4EST_QMAXLEVEL + 1));
  loop_args->quad_idx2 = P4EST_ALLOC (int, alloc_size / 2);
  loop_args->quadrants = P4EST_ALLOC (sc_array_t *, alloc_size);
  loop_args->index = P4EST_ALLOC (size_t *, alloc_size);
  for (i = 0; i < alloc_size; i++) {
    loop_args->index[i] =
      P4EST_ALLOC (size_t, (P4EST_QMAXLEVEL + 1) * P4EST_ITER_STRIDE);
    if (i & 1) {
      loop_args->quadrants[i] = &(ghost_layer->ghosts);
    }
  }
  loop_args->first_index = P4EST_ALLOC (size_t, alloc_size);
  loop_args->count = P4EST_ALLOC (size_t, alloc_size);
  loop_args->test = P4EST_ALLOC (p4est_quadrant_t *, alloc_size);
  loop_args->test_level = P4EST_ALLOC (int, alloc_size);
  loop_args->refine = P4EST_ALLOC (int8_t, alloc_size / 2);

  loop_args->tier_rings = p4est_iter_tier_rings_new (num_procs);

#ifdef P4_TO_P8
  loop_args->loop_edge = ((iter_corner != NULL) || (iter_edge != NULL));
#endif
  loop_args->loop_corner = (iter_corner != NULL);

  return loop_args;
}

static void
p4est_iter_loop_args_destroy (p4est_iter_loop_args_t * loop_args)
{
  int                 i;
  int                 alloc_size = loop_args->alloc_size;

  P4EST_FREE (loop_args->level_num);
  P4EST_FREE (loop_args->quad_idx2);
  P4EST_FREE (loop_args->quadrants);
  for (i = 0; i < alloc_size; i++) {
    P4EST_FREE (loop_args->index[i]);
  }
  P4EST_FREE (loop_args->index);
  P4EST_FREE (loop_args->first_index);
  P4EST_FREE (loop_args->count);
  P4EST_FREE (loop_args->test);
  P4EST_FREE (loop_args->test_level);
  P4EST_FREE (loop_args->refine);
  p4est_iter_tier_rings_destroy (loop_args->tier_rings);
  P4EST_FREE (loop_args);
}

/* initialize loop_args for a volume search */
static void
p4est_iter_init_loop_volume (p4est_iter_loop_args_t * loop_args,
                             p4est_topidx_t t, p4est_t * p4est,
                             p4est_ghost_t * ghost_layer)
{
  const int           left = 0;
  const int           right = 1;
  const int           local = 0;
  const int           ghost = 1;
  int                 i;
  sc_array_t         *trees = p4est->trees;
  p4est_tree_t       *tree;
  sc_array_t         *local_quads;
  size_t              first_ghost_quad =
    (size_t) ghost_layer->tree_offsets[t];
  size_t              stop_ghost_quad =
    (size_t) ghost_layer->tree_offsets[t + 1];

  tree = p4est_tree_array_index (trees, t);
  local_quads = &(tree->quadrants);

  loop_args->level = 0;
  loop_args->level_num[0] = 0;

  for (i = left; i <= right; i++) {
    loop_args->index[i * 2 + local][0] = 0;
    loop_args->index[i * 2 + local][1] = local_quads->elem_count;
    loop_args->index[i * 2 + ghost][0] = first_ghost_quad;
    loop_args->index[i * 2 + ghost][1] = stop_ghost_quad;
  }
  for (i = 0; i < 4; i++) {
    loop_args->quadrants[i] = (i & 1) ? &(ghost_layer->ghosts) : local_quads;
  }
#ifdef P4_TO_P8
  if (loop_args->loop_edge) {
    for (; i < 8; i++) {
      loop_args->quadrants[i] =
        (i & 1) ? &(ghost_layer->ghosts) : local_quads;
    }
  }
#endif
  if (loop_args->loop_corner) {
    for (; i < 2 * P4EST_CHILDREN; i++) {
      loop_args->quadrants[i] =
        (i & 1) ? &(ghost_layer->ghosts) : local_quads;
    }
  }
}

/* initialize loop_args for a face search between trees */
static void
p4est_iter_init_loop_face (p4est_iter_loop_args_t * loop_args,
                           p4est_topidx_t t, p4est_topidx_t nt,
                           p4est_t * p4est, p4est_ghost_t * ghost_layer)
{
  const int           left = 0;
  const int           right = 1;
  const int           local = 0;
  const int           ghost = 1;
  int                 i;
  sc_array_t         *trees = p4est->trees;
  p4est_tree_t       *tree;
  sc_array_t         *left_local_quads;
  sc_array_t         *right_local_quads;
  size_t              left_first_ghost =
    (size_t) ghost_layer->tree_offsets[t];
  size_t              left_stop_ghost =
    (size_t) ghost_layer->tree_offsets[t + 1];
  size_t              right_first_ghost =
    (size_t) ghost_layer->tree_offsets[nt];
  size_t              right_stop_ghost =
    (size_t) ghost_layer->tree_offsets[nt + 1];

  tree = p4est_tree_array_index (trees, t);
  left_local_quads = &(tree->quadrants);
  tree = p4est_tree_array_index (trees, nt);
  right_local_quads = &(tree->quadrants);

  loop_args->level = 0;
  loop_args->level_num[0] = 0;

  loop_args->index[left * 2 + local][0] = 0;
  loop_args->index[left * 2 + local][1] = left_local_quads->elem_count;
  loop_args->index[left * 2 + ghost][0] = left_first_ghost;
  loop_args->index[left * 2 + ghost][1] = left_stop_ghost;

  loop_args->index[right * 2 + local][0] = 0;
  loop_args->index[right * 2 + local][1] = right_local_quads->elem_count;
  loop_args->index[right * 2 + ghost][0] = right_first_ghost;
  loop_args->index[right * 2 + ghost][1] = right_stop_ghost;

  loop_args->quadrants[left * 2 + local] = left_local_quads;
  loop_args->quadrants[left * 2 + ghost] = &(ghost_layer->ghosts);
  loop_args->quadrants[right * 2 + local] = right_local_quads;
  loop_args->quadrants[right * 2 + ghost] = &(ghost_layer->ghosts);
  i = 4;
#ifdef P4_TO_P8
  if (loop_args->loop_edge) {
    for (; i < 8; i++) {
      loop_args->quadrants[i] = (i & 1) ? &(ghost_layer->ghosts) :
        ((i >> 1) & 1) ? right_local_quads : left_local_quads;
    }
  }
#endif
  if (loop_args->loop_corner) {
    for (; i < 2 * P4EST_CHILDREN; i++) {
      loop_args->quadrants[i] = (i & 1) ? &(ghost_layer->ghosts) :
        ((i >> 1) & 1) ? right_local_quads : left_local_quads;
    }
  }
}

/* initialize loop_args for a face search on a boundary face */
static void
p4est_iter_init_loop_outside_face (p4est_iter_loop_args_t * loop_args,
                                   p4est_topidx_t t, p4est_t * p4est,
                                   p4est_ghost_t * ghost_layer)
{
  const int           local = 0;
  const int           ghost = 1;
  int                 i;
  sc_array_t         *trees = p4est->trees;
  p4est_tree_t       *tree;
  sc_array_t         *local_quads;
  size_t              first_ghost_quad =
    (size_t) ghost_layer->tree_offsets[t];
  size_t              stop_ghost_quad =
    (size_t) ghost_layer->tree_offsets[t + 1];

  tree = p4est_tree_array_index (trees, t);
  local_quads = &(tree->quadrants);

  loop_args->level = 0;
  loop_args->level_num[0] = 0;

  loop_args->index[local][0] = 0;
  loop_args->index[local][1] = local_quads->elem_count;
  loop_args->index[ghost][0] = first_ghost_quad;
  loop_args->index[ghost][1] = stop_ghost_quad;

  loop_args->quadrants[local] = local_quads;
  loop_args->quadrants[ghost] = &(ghost_layer->ghosts);
  i = 2;
#ifdef P4_TO_P8
  if (loop_args->loop_edge) {
    for (; i < 4; i++) {
      loop_args->quadrants[i] =
        (i & 1) ? &(ghost_layer->ghosts) : local_quads;
    }
  }
#endif
  if (loop_args->loop_corner) {
    for (; i < P4EST_CHILDREN; i++) {
      loop_args->quadrants[i] =
        (i & 1) ? &(ghost_layer->ghosts) : local_quads;
    }
  }
}

/* initialize loop_args for an edge search between trees */
#ifdef P4_TO_P8
static void
p8est_iter_init_loop_edge (p4est_iter_loop_args_t * loop_args,
                           p8est_t * p8est, p4est_ghost_t * ghost_layer,
                           p8est_iter_edge_info_t * info)
{
  const int           local = 0;
  const int           ghost = 1;

  size_t              zz;
  size_t              limit = info->sides.elem_count;
  p8est_iter_edge_side_t *side;
  p4est_topidx_t      t;
  sc_array_t         *trees = p8est->trees;
  p8est_tree_t       *tree;
  sc_array_t         *local_quads;

  loop_args->level = 0;
  loop_args->level_num[0] = 0;

  for (zz = 0; zz < limit; zz++) {
    side = p8est_iter_eside_array_index (&info->sides, zz);
    t = side->treeid;
    tree = p4est_tree_array_index (trees, t);
    local_quads = &(tree->quadrants);
    loop_args->index[zz * 2 + local][0] = 0;
    loop_args->index[zz * 2 + local][1] = local_quads->elem_count;
    loop_args->index[zz * 2 + ghost][0] =
      (size_t) ghost_layer->tree_offsets[t];
    loop_args->index[zz * 2 + ghost][1] =
      (size_t) ghost_layer->tree_offsets[t + 1];
    loop_args->quadrants[zz * 2 + local] = local_quads;
    loop_args->quadrants[zz * 2 + ghost] = &(ghost_layer->ghosts);
    if (loop_args->loop_corner) {
      loop_args->quadrants[(limit + zz) * 2 + local] = local_quads;
      loop_args->quadrants[(limit + zz) * 2 + ghost] = &(ghost_layer->ghosts);
    }
  }
}
#endif

/* initialize loop_args for a corner search between trees */
static void
p4est_iter_init_loop_corner (p4est_iter_loop_args_t * loop_args,
                             p4est_t * p4est, p4est_ghost_t * ghost_layer,
                             p4est_iter_corner_info_t * info)
{
  const int           local = 0;
  const int           ghost = 1;

  size_t              zz;
  size_t              limit = info->sides.elem_count;
  p4est_iter_corner_side_t *side;
  p4est_topidx_t      t;
  sc_array_t         *trees = p4est->trees;
  p4est_tree_t       *tree;
  sc_array_t         *local_quads;

  loop_args->level = 0;
  loop_args->level_num[0] = 0;

  for (zz = 0; zz < limit; zz++) {
    side = p4est_iter_cside_array_index (&info->sides, zz);
    t = side->treeid;
    tree = p4est_tree_array_index (trees, t);
    local_quads = &(tree->quadrants);
    loop_args->index[zz * 2 + local][0] = 0;
    loop_args->index[zz * 2 + local][1] = local_quads->elem_count;
    loop_args->index[zz * 2 + ghost][0] =
      (size_t) ghost_layer->tree_offsets[t];
    loop_args->index[zz * 2 + ghost][1] =
      (size_t) ghost_layer->tree_offsets[t + 1];
    loop_args->quadrants[zz * 2 + local] = local_quads;
    loop_args->quadrants[zz * 2 + ghost] = &(ghost_layer->ghosts);
  }
}

/* When one iterate loop calls another, e.g. volume_iterate calls face_iterate,
 * the initial bounds for the new sides of the search need to be initialized.
 * The whole heirarchy doesn't need to be copied, just the most recent bounds
 * from the correct starting sections (start_idx2).
 */
static void
p4est_iter_copy_indices (p4est_iter_loop_args_t * loop_args,
                         const int *start_idx2, int old_num, int factor)
{
  const int           local = 0;
  const int           ghost = 1;
  int                 r, s, t;
  size_t            **zindex = loop_args->index;
  int                 i, j, idx2 = loop_args->level * P4EST_ITER_STRIDE;

  for (r = old_num, j = 1; j < factor; j++) {
    for (s = 0; s < old_num; s++, r++) {
      i = idx2 + start_idx2[r];
      for (t = local; t <= ghost; t++) {
        zindex[r * 2 + t][i] = zindex[s * 2 + t][i];
        zindex[r * 2 + t][i + 1] = zindex[s * 2 + t][i + 1];
      }
    }
  }
}

/* corner iterate function */
typedef struct p4est_iter_corner_args
{
  int                 num_sides;
  int                *start_idx2;
  int                 remote;
  p4est_iter_loop_args_t *loop_args;
  p4est_iter_corner_info_t info;
}
p4est_iter_corner_args_t;

static int
cside_compare (const void *a, const void *b)
{
  p4est_locidx_t      diff;
  const p4est_iter_corner_side_t *A = (const p4est_iter_corner_side_t *) a;
  const p4est_iter_corner_side_t *B = (const p4est_iter_corner_side_t *) b;

  diff = (A->treeid - B->treeid);
  if (diff) {
    return (int) diff;
  }
  return (A->corner - B->corner);
}

/* initialize corner args for a corner between trees */
static void
p4est_iter_init_corner (p4est_iter_corner_args_t * args,
                        p4est_t * p4est, p4est_ghost_t * ghost_layer,
                        p4est_iter_loop_args_t * loop_args, p4est_topidx_t t,
                        int c)
{
  p4est_topidx_t      ti;
  int                 i, j, k;
  int                 f, nf, o;
  int                 c2, nc;
  int                 count = 0;
  p4est_topidx_t      nt;
  p4est_connectivity_t *conn = p4est->connectivity;
  p4est_topidx_t     *ttt = conn->tree_to_tree;
  int8_t             *ttf = conn->tree_to_face;
  p4est_topidx_t     *ttc = conn->tree_to_corner;
  p4est_topidx_t     *ctt_offset = conn->ctt_offset;
  p4est_topidx_t     *ctt = conn->corner_to_tree;
  int8_t             *ctc = conn->corner_to_corner;
  p4est_topidx_t      corner = ttc != NULL ? ttc[t * P4EST_CHILDREN + c] : -1;
#ifdef P4_TO_P8
  int                 ref, set, nc2, orig_o;
  int                 e, ne, l;
  p4est_topidx_t     *tte = conn->tree_to_edge;
  p4est_topidx_t     *ett_offset = conn->ett_offset;
  p4est_topidx_t     *ett = conn->edge_to_tree;
  int8_t             *ete = conn->edge_to_edge;
  p4est_topidx_t      edge;
#endif
  p4est_iter_corner_info_t *info = &(args->info);
  p4est_iter_corner_side_t *cside;
  int                *start_idx2;

  info->p4est = p4est;
  info->ghost_layer = ghost_layer;
  info->tree_boundary = (int8_t) P4EST_CONNECT_CORNER;
  sc_array_init (&(info->sides), sizeof (p4est_iter_corner_side_t));
  start_idx2 = args->start_idx2 =
    P4EST_ALLOC (int, loop_args->alloc_size / 2);
  args->loop_args = loop_args;

  if (corner >= 0) {
    for (ti = ctt_offset[corner]; ti < ctt_offset[corner + 1]; ti++) {
      nt = ctt[ti];
      nc = (int) ctc[ti];
      cside = (p4est_iter_corner_side_t *) sc_array_push (&(info->sides));
      cside->corner = (int8_t) nc;
      cside->treeid = nt;
      start_idx2[count++] = 0;
      for (j = 0; j < P4EST_DIM; j++) {
        cside->faces[j] = -1;
#ifdef P4_TO_P8
        cside->edges[j] = -1;
#endif
      }
    }
  }
  else {
    cside = (p4est_iter_corner_side_t *) sc_array_push (&(info->sides));
    cside->corner = (int8_t) c;
    cside->treeid = t;
    for (j = 0; j < P4EST_DIM; j++) {
      cside->faces[j] = -1;
#ifdef P4_TO_P8
      cside->edges[j] = -1;
#endif
    }
    start_idx2[count++] = 0;
    for (i = 0; i < P4EST_DIM; i++) {
      f = p4est_corner_faces[c][i];
      c2 = p4est_corner_face_corners[c][f];
      nt = ttt[t * P4EST_FACES + f];
      nf = (int) ttf[t * P4EST_FACES + f];
      o = nf / P4EST_FACES;
      nf %= P4EST_FACES;
      if (nt == t && nf == f) {
        continue;
      }
#ifndef P4_TO_P8
      nc = p4est_face_corners[nf][(o == 0) ? c2 : (1 - c2)];
#else
      ref = p8est_face_permutation_refs[f][nf];
      set = p8est_face_permutation_sets[ref][o];
      nc2 = p8est_face_permutations[set][c2];
      nc = p8est_face_corners[nf][nc2];
#endif
      for (j = 0; j < count; j++) {
        cside = p4est_iter_cside_array_index_int (&(info->sides), j);
        if (cside->treeid == nt && (int) cside->corner == nc) {
          break;
        }
      }
      if (j < count) {
        continue;
      }
      cside = (p4est_iter_corner_side_t *) sc_array_push (&(info->sides));
      cside->corner = (int8_t) nc;
      cside->treeid = nt;
      for (j = 0; j < P4EST_DIM; j++) {
        cside->faces[j] = -1;
#ifdef P4_TO_P8
        cside->edges[j] = -1;
#endif
      }
      start_idx2[count++] = 0;
    }
#ifdef P4_TO_P8
    for (i = 0; i < 3; i++) {
      e = p8est_corner_edges[c][i];
      c2 = (p8est_edge_corners[e][0] == c) ? 0 : 1;
      edge = (tte != NULL) ? tte[t * 12 + e] : -1;
      if (edge >= 0) {
        orig_o = -1;
        for (ti = ett_offset[edge]; ti < ett_offset[edge + 1]; ti++) {
          nt = ett[ti];
          ne = (int) ete[ti];
          o = ne / 12;
          ne %= 12;
          if (nt == t && ne == e) {
            orig_o = o;
            break;
          }
        }
        P4EST_ASSERT (orig_o >= 0);
        P4EST_ASSERT (ti < ett_offset[edge + 1]);
        for (ti = ett_offset[edge]; ti < ett_offset[edge + 1]; ti++) {
          nt = ett[ti];
          ne = (int) ete[ti];
          o = ne / 12;
          ne %= 12;
          if (nt == t && ne == e) {
            continue;
          }
          nc = p8est_edge_corners[ne][(o == orig_o) ? c2 : (1 - c2)];
          for (j = 0; j < count; j++) {
            cside = p4est_iter_cside_array_index_int (&(info->sides), j);
            if (cside->treeid == nt && (int) cside->corner == nc) {
              break;
            }
          }
          if (j < count) {
            continue;
          }
          cside = (p4est_iter_corner_side_t *) sc_array_push (&(info->sides));
          cside->corner = (int8_t) nc;
          cside->treeid = nt;
          for (j = 0; j < P4EST_DIM; j++) {
            cside->faces[j] = -1;
#ifdef P4_TO_P8
            cside->edges[j] = -1;
#endif
          }
          start_idx2[count++] = 0;
        }
      }
    }
#endif
  }

  sc_array_sort (&(info->sides), cside_compare);

  {
    int                 faces_count = 0;
#ifdef P4_TO_P8
    int                 edges_count = 0;
#endif
    for (i = 0; i < count; i++) {
      cside = p4est_iter_cside_array_index_int (&(info->sides), i);
      nt = cside->treeid;
      nc = cside->corner;
      for (j = 0; j < P4EST_DIM; j++) {
        if (cside->faces[j] == -1) {
          int                 nnc;
          p4est_topidx_t      nnt;
          p4est_iter_corner_side_t *cside2;

          cside->faces[j] = faces_count;

          f = p4est_corner_faces[nc][j];
          nnt = ttt[P4EST_FACES * nt + f];
          c2 = p4est_corner_face_corners[nc][f];
          nf = (int) ttf[P4EST_FACES * nt + f];
          o = nf / P4EST_FACES;

          nf %= P4EST_FACES;
          if (nnt == nt && nf == f) {
            faces_count++;
            continue;
          }

#ifndef P4_TO_P8
          nnc = p4est_face_corners[nf][(o == 0) ? c2 : (1 - c2)];
#else
          ref = p8est_face_permutation_refs[f][nf];
          set = p8est_face_permutation_sets[ref][o];
          nc2 = p8est_face_permutations[set][c2];
          nnc = p8est_face_corners[nf][nc2];
#endif

          for (k = 0; k < count; k++) {
            if (k == i) {
              continue;
            }
            cside2 = p4est_iter_cside_array_index_int (&(info->sides), k);
            if (cside2->treeid == nnt && (int) cside2->corner == nnc) {
              cside2->faces[nf / 2] = faces_count;
            }
          }
          faces_count++;
        }
#ifdef P4_TO_P8
        if (cside->edges[j] == -1) {
          p4est_topidx_t      nnt;

          cside->edges[j] = edges_count;
          e = p8est_corner_edges[nc][j];
          c2 = (p8est_edge_corners[e][0] == nc) ? 0 : 1;
          edge = (tte != NULL) ? tte[nt * 12 + e] : -1;
          if (edge >= 0) {
            orig_o = -1;
            for (ti = ett_offset[edge]; ti < ett_offset[edge + 1]; ti++) {
              nnt = ett[ti];
              ne = (int) ete[ti];
              o = ne / 12;
              ne %= 12;
              if (nnt == nt && ne == e) {
                orig_o = o;
                break;
              }
            }
            P4EST_ASSERT (orig_o >= 0);
            P4EST_ASSERT (ti < ett_offset[edge + 1]);
            for (ti = ett_offset[edge]; ti < ett_offset[edge + 1]; ti++) {
              int                 nnc;

              nnt = ett[ti];
              ne = (int) ete[ti];
              o = ne / 12;
              ne %= 12;
              if (nnt == nt && ne == e) {
                continue;
              }
              nnc = p8est_edge_corners[ne][(o == orig_o) ? c2 : (1 - c2)];
              for (k = 0; k < count; k++) {
                p4est_iter_corner_side_t *cside2;

                cside2 = p4est_iter_cside_array_index_int (&(info->sides), k);
                if (cside2->treeid == nnt && (int) cside2->corner == nnc) {
                  cside2->edges[ne / 4] = edges_count;
                  break;
                }
              }
              P4EST_ASSERT (k < count);
            }
          }
          else {
            for (l = 0; l < 2; l++) {
              int                 nnc, e2;
              p4est_iter_corner_side_t *cside2;

              f = p8est_edge_faces[e][l];
              nnt = ttt[P4EST_FACES * nt + f];
              c2 = p4est_corner_face_corners[nc][f];
              nf = (int) ttf[P4EST_FACES * nt + f];
              o = nf / P4EST_FACES;

              nf %= P4EST_FACES;
              if (nnt == nt && nf == f) {
                continue;
              }
              ref = p8est_face_permutation_refs[f][nf];
              set = p8est_face_permutation_sets[ref][o];
              nc2 = p8est_face_permutations[set][c2];
              nnc = p8est_face_corners[nf][nc2];
              c2 = p8est_edge_corners[e][0] == nc ? p8est_edge_corners[e][1] :
                p8est_edge_corners[e][0];
              c2 = p8est_corner_face_corners[c2][f];
              nc2 = p8est_face_permutations[set][c2];
              nc2 = p8est_face_corners[nf][nc2];
              e2 = p8est_child_corner_edges[nnc][nc2];
              P4EST_ASSERT (e2 >= 0);
              for (k = 0; k < count; k++) {
                if (k == i) {
                  continue;
                }
                cside2 = p4est_iter_cside_array_index_int (&(info->sides), k);
                if (cside2->treeid == nnt && (int) cside2->corner == nnc) {
                  cside2->edges[e2 / 4] = edges_count;

                }
              }
            }
          }
          edges_count++;
        }
#endif
      }
    }
  }

  args->num_sides = count;
  p4est_iter_init_loop_corner (loop_args, p4est, ghost_layer, info);
}

static void
p4est_iter_reset_corner (p4est_iter_corner_args_t * args)
{
  sc_array_reset (&(args->info.sides));
  P4EST_FREE (args->start_idx2);
}

static void
p4est_corner_iterate (p4est_iter_corner_args_t * args, void *user_data,
                      p4est_iter_corner_t iter_corner)
{
  const int           local = 0;
  const int           ghost = 1;

  int                 side, st;

  p4est_iter_loop_args_t *loop_args = args->loop_args;
  int                 level = loop_args->level;
  int                 num_sides = args->num_sides;
  const int          *start_idx2 = args->start_idx2;
  int                *quad_idx2 = loop_args->quad_idx2;
  int                 this_corner;
  sc_array_t        **quadrants = loop_args->quadrants;
  size_t            **zindex = loop_args->index;
  size_t             *first_index = loop_args->first_index;
  size_t             *count = loop_args->count;
  p4est_quadrant_t  **test = loop_args->test;
  p4est_quadrant_t    temp;
  p4est_qcoord_t      mask =
    ((p4est_qcoord_t) - 1) << (P4EST_MAXLEVEL - level);
  sc_array_t          test_view;
  p4est_iter_corner_info_t *info = &(args->info);
  p4est_iter_corner_side_t *cside;
  ssize_t             temp_idx;
  int                 level_idx2;
  int                 type;
  int8_t              has_local;
  int                 lmax = 0;

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
      st = side * 2 + type;
      first_index[st] = zindex[st][quad_idx2[side]];
      count[st] = (zindex[st][quad_idx2[side] + 1] - first_index[st]);
      test[st] = NULL;
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

  has_local = 0;
  for (side = 0; side < num_sides; side++) {

    cside = p4est_iter_cside_array_index_int (&info->sides, side);
    cside->quad = NULL;
    cside->is_ghost = 1;
    cside->quadid = -1;
    this_corner = (int) cside->corner;
    for (type = local; type <= ghost; type++) {
      st = side * 2 + type;
      /* if we already found something locally, there's no need to search the
       * ghost layer */
      if (test[side * 2 + local] != NULL) {
        continue;
      }

      /* for this sidetype, we must find the most likely candidate in the
       * search area for touching the desired corner */
      if (count[st]) {
        /* get a candidate */
        if (count[st] == 1) {
          test[st] = p4est_quadrant_array_index (quadrants[st],
                                                 first_index[st]);
          temp_idx = 0;
        }
        else {
          switch (this_corner) {
          case (P4EST_CHILDREN - 1):
            test[st] = p4est_quadrant_array_index (quadrants[st],
                                                   first_index[st] +
                                                   count[st] - 1);
            temp_idx = ((ssize_t) count[st]) - 1;
            break;
          default:
            P4EST_ASSERT (first_index[st] < quadrants[st]->elem_count);
            test[st] = p4est_quadrant_array_index (quadrants[st],
                                                   first_index[st]);
            temp_idx = 0;
            break;
          }
        }
        /* create the smallest quadrant in the appropriate corner */
        temp = *(test[st]);
        temp.level = (int8_t) level;
        temp.x &= mask;
        temp.y &= mask;
#ifdef P4_TO_P8
        temp.z &= mask;
#endif
        p4est_quadrant_corner_descendant (&temp, &temp, this_corner,
                                          P4EST_QMAXLEVEL);
        P4EST_ASSERT (p4est_quadrant_is_valid (&temp));
        /* we do not have to search if there is one quadrant, or if we are in
         * the first or last corner */
        if (count[st] == 1 || this_corner == 0 ||
            this_corner == P4EST_CHILDREN - 1) {
          /* if test[sidetype] does not contain temp */
          if ((!p4est_quadrant_is_equal (test[st], &temp)) &&
              (!p4est_quadrant_is_ancestor (test[st], &temp))) {
            test[st] = NULL;
          }
        }
        else {
          /* we search for the quadrant containing temp */
          sc_array_init_view (&test_view, quadrants[st],
                              first_index[st], count[st]);
          temp_idx = sc_array_bsearch (&test_view, &temp,
                                       p4est_quadrant_disjoint);
          /* if there is no quadrant containing temp, then no quad in the
           * search area can touch the corner */
          if (temp_idx == -1) {
            test[st] = NULL;
          }
          else {
            P4EST_ASSERT (temp_idx >= 0);
            test[st] = p4est_quadrant_array_index (&test_view,
                                                   (size_t) temp_idx);
          }
        }
        /* if we have found the right quadrant for this side of the corner */
        if (test[st] != NULL) {
          P4EST_ASSERT (p4est_quadrant_overlaps (test[st], &temp));
          P4EST_ASSERT (temp_idx >= 0 && (size_t) temp_idx < count[st]);
          temp_idx += first_index[st];
          cside->quad = test[st];
          cside->is_ghost = (type == ghost);
          cside->quadid = (p4est_locidx_t) temp_idx;
          if (type == local) {
            has_local = 1;
          }
          lmax = SC_MAX (lmax, test[st]->level);
        }
      }
    }
  }
  if (has_local) {
    /* always run if a local quadrant touches the corner */
    iter_corner (info, user_data);
  }
  else if (args->remote) {
    /* otherwise, if are supposed to run the callback if the corner is
     * remotely touched, we have to determine if the corner is remotely
     * touched */
    p4est_t            *p4est = info->p4est;
    p4est_quadrant_t   *fq = &(p4est->global_first_position[p4est->mpirank]);
    p4est_quadrant_t   *lq =
      &(p4est->global_first_position[p4est->mpirank + 1]);
    p4est_topidx_t      flt = fq->p.which_tree;
    p4est_topidx_t      llt = lq->p.which_tree;

    for (side = 0; side < num_sides; side++) {
      cside = p4est_iter_cside_array_index_int (&info->sides, side);
      p4est_quadrant_t   *q = cside->quad;
      int                 dir;
      int                 l = -1;

      if (q == NULL || (l = q->level) == lmax) {
        continue;
      }
      for (dir = 0; dir < P4EST_DIM; dir++) {
        int                 f = cside->faces[dir];
        int                 side2;
        int                 ht;

        for (side2 = 0; side2 < num_sides; side2++) {
          p4est_iter_corner_side_t *cside2;
          int                 dir2;

          if (side2 == side) {
            continue;
          }

          cside2 = p4est_iter_cside_array_index_int (&info->sides, side2);

          ht = cside2->treeid;

          for (dir2 = 0; dir2 < P4EST_DIM; dir2++) {
            if (cside2->faces[dir2] == f && cside2->quad) {
              p4est_quadrant_t   *q2 = cside2->quad;
              int                 i;
              int                 nf =
                p4est_corner_faces[cside2->corner][dir2];

              if (q2->level < l + 1) {
                continue;
              }

              for (i = 0; i < P4EST_HALF; i++) {
                int                 nc = p4est_face_corners[nf][i];
                p4est_quadrant_t    h;
                if (nc == cside2->corner) {
                  continue;
                }
                p4est_quadrant_sibling (q2, &h, nc);
                if ((ht > flt ||
                     (ht == flt && p4est_quadrant_disjoint (fq, &h) <= 0))
                    &&
                    (ht < llt ||
                     (ht == llt && p4est_quadrant_disjoint (&h, lq) < 0))) {
                  iter_corner (info, user_data);
                  return;
                }
              }
            }
          }
        }
      }
#ifdef P4_TO_P8
      for (dir = 0; dir < P4EST_DIM; dir++) {
        int                 e = cside->edges[dir];
        int                 side2;
        int                 ht;

        for (side2 = 0; side2 < num_sides; side2++) {
          p4est_iter_corner_side_t *cside2;
          int                 dir2;

          if (side2 == side) {
            continue;
          }
          cside2 = p4est_iter_cside_array_index_int (&info->sides, side2);
          ht = cside2->treeid;
          for (dir2 = 0; dir2 < P4EST_DIM; dir2++) {
            if (cside2->edges[dir2] == e && cside2->quad) {
              p4est_quadrant_t    h;
              p4est_quadrant_t   *q2 = cside2->quad;
              int                 ne =
                p8est_corner_edges[cside2->corner][dir2];

              if (q2->level < l + 1) {
                continue;
              }

              if (p8est_edge_corners[ne][0] == cside2->corner) {
                p4est_quadrant_sibling (q2, &h, p8est_edge_corners[ne][1]);
              }
              else {
                p4est_quadrant_sibling (q2, &h, p8est_edge_corners[ne][0]);
              }
              if ((ht > flt ||
                   (ht == flt && p4est_quadrant_disjoint (fq, &h) <= 0))
                  &&
                  (ht < llt ||
                   (ht == llt && p4est_quadrant_disjoint (&h, lq) < 0))) {
                /* run the callback */
                iter_corner (info, user_data);
                return;
              }
            }
          }
        }
      }
#endif
    }
  }
}

/* edge iterate functions */
#ifdef P4_TO_P8
typedef struct p8est_iter_edge_args
{
  int                 num_sides;
  int                *start_idx2;
  sc_array_t          common_corners[2];        /* for each side of the edge,
                                                   there are two corners that
                                                   touch the edge */
  p4est_iter_loop_args_t *loop_args;
  p4est_iter_corner_args_t corner_args;
  p8est_iter_edge_info_t info;
  int                 remote;
}
p8est_iter_edge_args_t;

/* given valid edge arguments, setup the corner arguments for a corner search
 * that is called for the corner where two adjacent, colinear edges meet */
static void
p8est_iter_init_corner_from_edge (p4est_iter_corner_args_t * args,
                                  p8est_iter_edge_args_t * edge_args)
{
  int                 j, k, dir, ndir1, ndir2;;
  p8est_iter_corner_info_t *info = &(args->info);
  p8est_iter_edge_side_t *eside;
  p8est_iter_corner_side_t *cside;
  sc_array_t         *common_corners = edge_args->common_corners;
  int                *c_start_idx2;
  int                 n_edge_faces = 0;

  info->p4est = edge_args->info.p4est;
  info->ghost_layer = edge_args->info.ghost_layer;
  info->tree_boundary = edge_args->info.tree_boundary;
  sc_array_init (&(info->sides), sizeof (p4est_iter_corner_side_t));
  args->loop_args = edge_args->loop_args;
  args->num_sides = edge_args->num_sides * 2;
  c_start_idx2 = args->start_idx2 = P4EST_ALLOC (int, args->num_sides);
  sc_array_resize (&(info->sides), (size_t) args->num_sides);

  for (j = 0; j < edge_args->num_sides; j++) {
    eside = p8est_iter_eside_array_index_int (&edge_args->info.sides, j);
    n_edge_faces = SC_MAX (n_edge_faces, eside->faces[0] + 1);
    n_edge_faces = SC_MAX (n_edge_faces, eside->faces[1] + 1);
  }
  for (j = 0; j < args->num_sides; j++) {
    k = j % edge_args->num_sides;
    eside = p8est_iter_eside_array_index_int (&edge_args->info.sides, k);
    cside = p4est_iter_cside_array_index_int (&info->sides, j);
    cside->treeid = eside->treeid;
    dir = eside->edge / 4;
    cside->faces[dir] = 2 * n_edge_faces + k;
    ndir1 = SC_MIN (((dir + 1) % 3), ((dir + 2) % 3));
    ndir2 = SC_MAX (((dir + 1) % 3), ((dir + 2) % 3));
    cside->edges[ndir1] = eside->faces[1];
    cside->edges[ndir2] = eside->faces[0];
    if (j == k) {
      cside->corner = (int8_t) * ((int *) sc_array_index_int
                                  (&(common_corners[1]), k));
      c_start_idx2[j] = *((int *) sc_array_index_int (&(common_corners[0]),
                                                      k));
      cside->edges[dir] = n_edge_faces;
      cside->faces[ndir1] = eside->faces[0];
      cside->faces[ndir2] = eside->faces[1];
    }
    else {
      cside->corner = (int8_t) * ((int *) sc_array_index_int
                                  (&(common_corners[0]), k));
      c_start_idx2[j] = *((int *) sc_array_index_int (&(common_corners[1]),
                                                      k));
      cside->edges[dir] = n_edge_faces + 1;
      cside->faces[ndir1] = eside->faces[0] + n_edge_faces;
      cside->faces[ndir2] = eside->faces[1] + n_edge_faces;
    }
  }
  args->remote = edge_args->remote;
}

static int
eside_compare (const void *a, const void *b)
{
  p4est_locidx_t      diff;
  const p8est_iter_edge_side_t *A = (const p8est_iter_edge_side_t *) a;
  const p8est_iter_edge_side_t *B = (const p8est_iter_edge_side_t *) b;

  diff = (A->treeid - B->treeid);
  if (diff) {
    return (int) diff;
  }
  return (A->edge - B->edge);
}

/* initialize edge args for an edge between trees */
static void
p8est_iter_init_edge (p8est_iter_edge_args_t * args, p8est_t * p8est,
                      p4est_ghost_t * ghost_layer,
                      p4est_iter_loop_args_t * loop_args, p4est_topidx_t t,
                      int e)
{
  p4est_topidx_t      ti;
  int                 i, j, k;
  int                 f, nf, o, ref, set;
  int                 ne;
  int                 c0, c1, nc0, nc1, *cc;
  int                 count = 0;
  p4est_topidx_t      nt;
  p8est_connectivity_t *conn = p8est->connectivity;
  p8est_iter_edge_info_t *info = &(args->info);
  p8est_iter_edge_side_t *eside;
  int                *start_idx2;
  p4est_topidx_t     *ttt = conn->tree_to_tree;
  int8_t             *ttf = conn->tree_to_face;
  p4est_topidx_t     *tte = conn->tree_to_edge;
  p4est_topidx_t     *ett_offset = conn->ett_offset;
  p4est_topidx_t     *ett = conn->edge_to_tree;
  int8_t             *ete = conn->edge_to_edge;
  p4est_topidx_t      edge = (tte != NULL) ? tte[t * 12 + e] : -1;
  sc_array_t         *common_corners = args->common_corners;

  info->p4est = p8est;
  info->ghost_layer = ghost_layer;
  info->tree_boundary = (int8_t) P8EST_CONNECT_EDGE;
  start_idx2 = args->start_idx2 =
    P4EST_ALLOC (int, loop_args->alloc_size / 2);
  sc_array_init (&(info->sides), sizeof (p8est_iter_edge_side_t));
  sc_array_init (&(args->common_corners[0]), sizeof (int));
  sc_array_init (&(args->common_corners[1]), sizeof (int));
  args->loop_args = loop_args;

  if (edge >= 0) {
    for (ti = ett_offset[edge]; ti < ett_offset[edge + 1]; ti++) {
      nt = ett[ti];
      ne = (int) ete[ti];
      o = ne / 12;
      ne %= 12;
      eside = (p8est_iter_edge_side_t *) sc_array_push (&(info->sides));
      eside->orientation = (int8_t) o;
      eside->edge = (int8_t) ne;
      eside->treeid = nt;
      start_idx2[count++] = 0;
      eside->faces[0] = -1;
      eside->faces[1] = -1;
    }
  }
  else {
    eside = (p8est_iter_edge_side_t *) sc_array_push (&(info->sides));
    eside->edge = (int8_t) e;
    eside->treeid = t;
    eside->orientation = 0;
    start_idx2[count++] = 0;
    eside->faces[0] = -1;
    eside->faces[1] = -1;
    for (i = 0; i < 2; i++) {
      f = p8est_edge_faces[e][i];
      nt = ttt[t * P4EST_FACES + f];
      nf = (int) ttf[t * P4EST_FACES + f];
      o = nf / P4EST_FACES;
      nf %= P4EST_FACES;
      if (nt == t && nf == f) {
        continue;
      }
      c0 = p8est_edge_corners[e][0];
      c1 = p8est_edge_corners[e][1];
      c0 = p8est_corner_face_corners[c0][f];
      c1 = p8est_corner_face_corners[c1][f];
      ref = p8est_face_permutation_refs[f][nf];
      set = p8est_face_permutation_sets[ref][o];
      nc0 = p8est_face_permutations[set][c0];
      nc1 = p8est_face_permutations[set][c1];
      nc0 = p8est_face_corners[nf][nc0];
      nc1 = p8est_face_corners[nf][nc1];
      ne = p8est_child_corner_edges[nc0][nc1];
      for (j = 0; j < count; j++) {
        eside = p8est_iter_eside_array_index_int (&(info->sides), j);
        if (eside->treeid == nt && (int) eside->edge == ne) {
          break;
        }
      }
      if (j < count) {
        continue;
      }
      eside = (p8est_iter_edge_side_t *) sc_array_push (&(info->sides));
      eside->orientation = (int8_t) ((nc0 < nc1) ? 0 : 1);
      eside->edge = (int8_t) ne;
      eside->treeid = nt;
      start_idx2[count++] = 0;
      eside->faces[0] = -1;
      eside->faces[1] = -1;
    }
  }

  sc_array_sort (&(info->sides), eside_compare);
  for (i = 0; i < count; i++) {
    eside = p8est_iter_eside_array_index_int (&(info->sides), i);
    if (!i) {
      o = eside->orientation;
      eside->orientation = 0;
    }
    else {
      eside->orientation ^= o;
    }
    cc = (int *) sc_array_push (&(common_corners[0]));
    *cc = p8est_edge_corners[eside->edge][eside->orientation];
    cc = (int *) sc_array_push (&(common_corners[1]));
    *cc = p8est_edge_corners[eside->edge][1 - eside->orientation];
  }

  args->num_sides = count;

  {
    int                 faces_count = 0;
    for (i = 0; i < count; i++) {
      eside = p8est_iter_eside_array_index_int (&(info->sides), i);
      nt = eside->treeid;
      ne = eside->edge;
      for (j = 0; j < 2; j++) {
        if (eside->faces[j] == -1) {
          int                 nne;
          p4est_topidx_t      nnt;

          eside->faces[j] = faces_count;

          f = p8est_edge_faces[ne][j];
          nnt = ttt[P4EST_FACES * nt + f];
          nf = (int) ttf[nt * P4EST_FACES + f];
          o = nf / P4EST_FACES;
          nf %= P4EST_FACES;
          if (nnt == nt && nf == f) {
            faces_count++;
            continue;
          }
          c0 = p8est_edge_corners[ne][0];
          c1 = p8est_edge_corners[ne][1];
          c0 = p8est_corner_face_corners[c0][f];
          c1 = p8est_corner_face_corners[c1][f];
          ref = p8est_face_permutation_refs[f][nf];
          set = p8est_face_permutation_sets[ref][o];
          nc0 = p8est_face_permutations[set][c0];
          nc1 = p8est_face_permutations[set][c1];
          nc0 = p8est_face_corners[nf][nc0];
          nc1 = p8est_face_corners[nf][nc1];
          nne = p8est_child_corner_edges[nc0][nc1];

          for (k = 0; k < count; k++) {
            p8est_iter_edge_side_t *eside2;

            if (k == i) {
              continue;
            }
            eside2 = p8est_iter_eside_array_index_int (&(info->sides), k);

            if (eside2->treeid == nnt && (int) eside2->edge == nne) {
              if (p8est_edge_faces[nne][0] == nf) {
                eside2->faces[0] = faces_count;
              }
              else {
                P4EST_ASSERT (p8est_edge_faces[nne][1] == nf);
                eside2->faces[1] = faces_count;
              }
            }
          }
          faces_count++;
        }
      }
    }
  }

  if (loop_args->loop_corner) {
    p8est_iter_init_corner_from_edge (&(args->corner_args), args);
  }
  p8est_iter_init_loop_edge (loop_args, p8est, ghost_layer, info);
}

static void
p8est_iter_reset_edge (p8est_iter_edge_args_t * args)
{
  if (args->loop_args->loop_corner) {
    p4est_iter_reset_corner (&args->corner_args);
  }
  sc_array_reset (&(args->common_corners[0]));
  sc_array_reset (&(args->common_corners[1]));
  sc_array_reset (&(args->info.sides));
  P4EST_FREE (args->start_idx2);
}

static void
p8est_edge_iterate (p8est_iter_edge_args_t * args, void *user_data,
                    p8est_iter_edge_t iter_edge,
                    p8est_iter_corner_t iter_corner)
{
  const int           local = 0;
  const int           ghost = 1;

  p4est_iter_loop_args_t *loop_args = args->loop_args;
  int                 num_sides = args->num_sides;
  int                 start_level = loop_args->level;
  int                *start_idx2 = args->start_idx2;
  int                *level_num = loop_args->level_num;
  sc_array_t        **quadrants = loop_args->quadrants;
  size_t            **zindex = loop_args->index;
  size_t             *first_index = loop_args->first_index;
  sc_array_t         *common_corners = args->common_corners;
  p8est_quadrant_t  **test = loop_args->test;
  size_t             *count = loop_args->count;
  int                *test_level = loop_args->test_level;
  int                *quad_idx2 = loop_args->quad_idx2;
  int8_t             *refine = loop_args->refine;
  int                *temp_int, *temp_int2;
  int                 i;
  int                *Level = &(loop_args->level);
  int                 side;
  int                 type;
  int                 st;
  int                 level_idx2;
  p8est_iter_edge_info_t *info = &(args->info);
  p8est_iter_edge_side_t *eside;
  p8est_quadrant_t  **quads;
  p4est_locidx_t     *quadids;
  int8_t             *is_ghost;
  int                 child_corner;
  int8_t              has_local;
  sc_array_t          test_view;
  int8_t              all_empty, stop_refine;
  p4est_iter_corner_args_t *corner_args = &(args->corner_args);
  sc_array_t         *tier_rings = loop_args->tier_rings;

  /* level_idx2 moves us to the correct set of bounds within the index arrays
   * for the level: it is a set of bounds because it includes all children at
   * this level */
  level_idx2 = start_level * P4EST_ITER_STRIDE;
  for (side = 0; side < num_sides; side++) {

    /* start_idx2 gives the ancestor id at level for the search area on this
     * side, so quad_idx2[side] now gives the correct location in
     * index[sidetype] of the bounds of the search area */
    quad_idx2[side] = level_idx2 + start_idx2[side];

    /* get the location in quadrants[sidetype] of the first quadrant in the
     * search area, and the count of quadrants in the search area */
    for (type = local; type <= ghost; type++) {
      st = side * 2 + type;
      first_index[st] = zindex[st][quad_idx2[side]];
      count[st] = (zindex[st][quad_idx2[side] + 1] - first_index[st]);
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
        st = side * 2 + type;
        if (count[st]) {
          test[st] = p4est_quadrant_array_index (quadrants[st],
                                                 first_index[st]);
          test_level[st] = (int) test[st]->level;
        }
        else {
          test[st] = NULL;
          test_level[st] = -1;
        }
      }
      /* initially assume that every side needs to be refined */
      refine[side] = 1;
    }
    /* initially assume that we are going to have to refine our search areas */
    stop_refine = 0;
    has_local = 0;
    for (side = 0; side < num_sides; side++) {
      for (type = local; type <= ghost; type++) {
        st = side * 2 + type;
        /* if the candidate from sidetype is the same size as the search area,
         * then we do not refine this side */
        if (test_level[st] == *Level) {
          if (!stop_refine) {
            stop_refine = 1;
            /* we are not going to recur on the next level, instead moving to
             * the next branch on this level */
            level_num[*Level]++;
            /* if there is no edge callback (i.e., we are just running
             * edge_iterate to find corners), then we're done with this branch */
            if (iter_edge == NULL) {
              goto change_search_area;
            }
          }
          /* this side is full-sized, it does not need to be refined */
          refine[side] = 0;
          eside = p8est_iter_eside_array_index_int (&info->sides, side);
          eside->is_hanging = 0;
          eside->is.full.quad = test[st];
          eside->is.full.is_ghost = (type == ghost);
          eside->is.full.quadid = first_index[st];
          has_local = (has_local || (type == local));
        }
      }
    }
    for (side = 0; side < num_sides; side++) {
      if (refine[side]) {
        if (stop_refine && count[side * 2 + local] == 0 &&
            count[side * 2 + ghost] == 0) {
          /* if a side is empty, the appropriate quadrant(s) are missing from the
           * ghost layer, so we fill in with NULL */
          eside = p8est_iter_eside_array_index_int (&info->sides, side);
          eside->is_hanging = 0;
          eside->is.full.quad = NULL;
          eside->is.full.is_ghost = 1;
          eside->is.full.quadid = -1;
          refine[side] = 0;
        }
      }
      if (refine[side]) {
        /* at this point, the side needs to be refined, so we take the search area
         * and split it up, taking the indices for the refined search areas and
         * placing them on the next tier in index[sidetype] */
        quad_idx2[side] = level_idx2 + P4EST_ITER_STRIDE;
        for (type = local; type <= ghost; type++) {
          st = side * 2 + type;
          sc_array_init_view (&test_view, quadrants[st],
                              first_index[st], count[st]);
          p4est_iter_tier_insert (&test_view, *Level, zindex[st] +
                                  quad_idx2[side], first_index[st],
                                  tier_rings, test[st]);
        }
        if (stop_refine) {
          /* fill the side info */
          eside = p8est_iter_eside_array_index_int (&(info->sides), side);
          eside->is_hanging = 1;
          quads = eside->is.hanging.quad;
          is_ghost = eside->is.hanging.is_ghost;
          quadids = eside->is.hanging.quadid;
          for (i = 0; i < 2; i++) {
            /* esides expects things in z-order, not the search order */
            temp_int =
              (int *) sc_array_index_int (&(common_corners[i]), side);
            temp_int2 =
              (int *) sc_array_index_int (&(common_corners[1 - i]), side);
            if (*temp_int < *temp_int2) {
              P4EST_ASSERT (p8est_edge_corners[eside->edge][0] == *temp_int);
              child_corner = 0;
            }
            else {
              P4EST_ASSERT (p8est_edge_corners[eside->edge][1] == *temp_int);
              child_corner = 1;
            }
            quads[child_corner] = NULL;
            is_ghost[child_corner] = 1;
            quadids[child_corner] = -1;

            /* get the location in index[sidetype] of the index for the hanging
             * quadrant */
            quad_idx2[side] = level_idx2 + P4EST_ITER_STRIDE + *temp_int;
            for (type = local; type <= ghost; type++) {
              st = side * 2 + type;
              first_index[st] = zindex[st][quad_idx2[side]];
              count[st] =
                (size_t) zindex[st][quad_idx2[side] + 1] - first_index[st];
              /* if the search area is non-empty, by the two to one condition
               * it must contain exactly quadrant half as large as the search
               * level, which we add to the collection */
              if (count[st]) {
                quads[child_corner] = p4est_quadrant_array_index
                  (quadrants[st], first_index[st]);
                P4EST_ASSERT ((int) quads[child_corner]->level == *Level + 1);
                is_ghost[child_corner] = (type == ghost);
                quadids[child_corner] = (p4est_locidx_t) first_index[st];
                has_local = (has_local || (type == local));
              }
            }
          }
        }                       /* if (stop_refine) */
      }                         /* if (refine[side]) */
    }                           /* for (side = 0; side < num_sides; side++) */
    if (stop_refine) {
      if (has_local) {
        /* if there is a local quadrant, we run the callback */
        iter_edge (info, user_data);
      }
      else if (args->remote) {
        /* if we are supposed to run on remotely touched edges, determine if
         * the edge is remotely touched first */
        p4est_t            *p4est = info->p4est;
        p4est_quadrant_t   *fq =
          &(p4est->global_first_position[p4est->mpirank]);
        p4est_quadrant_t   *lq =
          &(p4est->global_first_position[p4est->mpirank + 1]);
        p4est_topidx_t      flt = fq->p.which_tree;
        p4est_topidx_t      llt = lq->p.which_tree;

        for (side = 0; side < num_sides; side++) {
          int                 dir;
          p4est_topidx_t      ht;

          eside = p8est_iter_eside_array_index_int (&info->sides, side);
          ht = eside->treeid;
          if (!eside->is_hanging || eside->is.hanging.quad[0] == NULL ||
              eside->is.hanging.quad[1] == NULL) {
            continue;
          }
          for (dir = 0; dir < 2; dir++) {
            int                 f = eside->faces[dir];
            int                 side2;

            for (side2 = 0; side2 < num_sides; side2++) {
              p8est_iter_edge_side_t *eside2;
              int                 dir2;

              if (side2 == side) {
                continue;
              }
              eside2 = p8est_iter_eside_array_index_int (&info->sides, side2);
              if (eside2->is_hanging || eside2->is.full.quad == NULL) {
                continue;
              }
              for (dir2 = 0; dir2 < 2; dir2++) {
                if (eside2->faces[dir2] == f) {
                  p4est_quadrant_t    h;
                  int                 c1, c2, fc[2];
                  int                 face;

                  face = p8est_edge_faces[eside->edge][dir];
                  c1 = p8est_edge_corners[eside->edge][0];
                  c2 = p8est_edge_corners[eside->edge][1];
                  fc[0] = p8est_corner_face_corners[c1][face];
                  fc[1] = p8est_corner_face_corners[c2][face];
                  fc[0] = p8est_face_corners[face][fc[0] ^ 3];
                  fc[1] = p8est_face_corners[face][fc[1] ^ 3];
                  for (i = 0; i < 2; i++) {
                    p4est_quadrant_sibling (eside->is.hanging.quad[0], &h,
                                            fc[i]);
                    if ((ht > flt
                         || (ht == flt
                             && p4est_quadrant_disjoint (fq, &h) <= 0))
                        && (ht < llt
                            || (ht == llt
                                && p4est_quadrant_disjoint (&h, lq) < 0))) {
                      iter_edge (info, user_data);
                      /* this goto is to avoid a 5-level break */
                      goto change_search_area;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }                           /* if (stop_refine) */
    else {
      /* if every side needed to be refined, then we descend along this branch to
       * the next level and search there */
      level_num[++(*Level)] = 0;
      level_idx2 += P4EST_ITER_STRIDE;
    }

  change_search_area:

    for (;;) {
      /* if we tried to advance the search area on start_level, we've completed
       * the search */
      if (level_num[start_level] > 0) {
        P4EST_ASSERT (*Level == start_level);
        return;
      }
      /* if we have tried to advance the search area past two branches, that
       * means that we have completed all of the branches on this level */
      if (level_num[*Level] == 2) {
        /* if we have a corner callback, we need to run it on the corner between
         * the edge branches on this level */
        if (loop_args->loop_corner) {
          P4EST_ASSERT (corner_args->num_sides == 2 * num_sides);
          p4est_iter_copy_indices (loop_args, corner_args->start_idx2,
                                   num_sides, 2);
          p4est_corner_iterate (corner_args, user_data, iter_corner);
        }
        /* now that we're done on this level, go up a level and over a branch */
        level_num[--(*Level)]++;
        level_idx2 -= P4EST_ITER_STRIDE;
      }
      else {
        /* at this point, we need to initialize the bounds of the search areas
         * for this new branch */
        all_empty = 1;
        for (side = 0; side < num_sides; side++) {
          temp_int = (int *) sc_array_index_int
            (&(common_corners[level_num[*Level]]), side);
          quad_idx2[side] = level_idx2 + *temp_int;
          for (type = local; type <= ghost; type++) {
            st = side * 2 + type;
            first_index[st] = zindex[st][quad_idx2[side]];
            count[st] = (zindex[st][quad_idx2[side] + 1] - first_index[st]);
            if (type == local && count[st]) {
              all_empty = 0;
            }
          }
        }
        if (all_empty) {
          /* if there are no local quadrants in any of the search areas, we're done
           * with this search area and proceed to the next branch on this level */
          level_num[*Level]++;
        }
        else {
          /* otherwise we are done changing the search area */
          break;
        }
      }
    }                           /* change_search_area for (;;) */
  }                             /* outer for (;;) */
}
#endif

/* face iterate functions */
typedef struct p4est_iter_face_args
{
  p4est_iter_loop_args_t *loop_args;
  int                 start_idx2[2];
  /* when a search branch is refined,
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
  int                 num_to_child[P4EST_CHILDREN];
  int8_t              outside_face;     /* indicates if we are at a tree
                                           boundary without a neighbor across
                                           the face */
#ifdef P4_TO_P8
  p8est_iter_edge_args_t edge_args[2][2];
#endif
  p4est_iter_corner_args_t corner_args;
  p4est_iter_face_info_t info;
  int                 remote;
}
p4est_iter_face_args_t;

/* given valid face arguments, setup corner arguments for a corner search that
 * is called for a corner where four adjacent, coplaner faces meet.
 */
static void
p4est_iter_init_corner_from_face (p4est_iter_corner_args_t * args,
                                  p4est_iter_face_args_t * face_args)
{
  const int           ntc_str = P4EST_CHILDREN / 2;
  int                 j, k;
  p4est_iter_corner_info_t *info = &(args->info);
  p4est_iter_face_side_t *fside;
  p4est_iter_corner_side_t *cside;
  int                *num_to_child = face_args->num_to_child;
  int                *c_start_idx2;
  int                 limit = face_args->outside_face ? 1 : 2;
  int                 count = 0;
  int                 dir;
#ifdef P4_TO_P8
  int                 ndir1, ndir2;
#endif

  info->p4est = face_args->info.p4est;
  info->ghost_layer = face_args->info.ghost_layer;
  info->tree_boundary = face_args->info.tree_boundary;
  sc_array_init (&(info->sides), sizeof (p4est_iter_corner_side_t));
  args->num_sides = limit * ntc_str;
  sc_array_resize (&(info->sides), (size_t) args->num_sides);
  c_start_idx2 = args->start_idx2 = P4EST_ALLOC (int, args->num_sides);
  args->loop_args = face_args->loop_args;

  for (j = 0; j < ntc_str; j++) {
    for (k = 0; k < limit; k++) {
      fside = p4est_iter_fside_array_index_int (&face_args->info.sides, k);
      cside = p4est_iter_cside_array_index_int (&info->sides, count);
      cside->treeid = fside->treeid;

      dir = fside->face / 2;
      cside->corner = (int8_t) num_to_child[k * ntc_str + (ntc_str - 1 - j)];
      c_start_idx2[count++] = num_to_child[k * ntc_str + j];
      cside->faces[dir] = j;
#ifndef P4_TO_P8
      cside->faces[dir ^ 1] = ntc_str + k;
#else
      cside->edges[dir] = ntc_str + k;
      ndir1 = SC_MIN (((dir + 1) % 3), ((dir + 2) % 3));
      ndir2 = SC_MAX (((dir + 1) % 3), ((dir + 2) % 3));
      if (!k) {
        cside->faces[ndir1] = ntc_str + ntc_str * k + (j >> 1);
        cside->faces[ndir2] = ntc_str + ntc_str * k + 2 + (j & 1);
        cside->edges[ndir1] = (j & 1);
        cside->edges[ndir2] = 2 + (j >> 1);
      }
      else {
        int                 c0 =
          p8est_corner_face_corners[num_to_child[ntc_str]][fside->face];
        int                 c1 =
          p8est_corner_face_corners[num_to_child[ntc_str + 1]][fside->face];

        if (c0 == (c1 ^ 1)) {
          cside->faces[ndir1] = ntc_str + ntc_str * k + (j >> 1);
          cside->faces[ndir2] = ntc_str + ntc_str * k + 2 + (j & 1);
          cside->edges[ndir1] = (j & 1);
          cside->edges[ndir2] = 2 + (j >> 1);
        }
        else {
          P4EST_ASSERT (c0 == (c1 ^ 2));
          cside->faces[ndir2] = ntc_str + ntc_str * k + (j >> 1);
          cside->faces[ndir1] = ntc_str + ntc_str * k + 2 + (j & 1);
          cside->edges[ndir2] = (j & 1);
          cside->edges[ndir1] = 2 + (j >> 1);
        }
      }
#endif
    }
  }
  args->remote = face_args->remote;
}

/* given valid face arguments, setup edge arguments for an edge search that is
 * called for an edge where two adjacent, coplaner faces meet: each plane
 * has two directions, and each direction can begin from one of two potential
 * starting sides.
 */
#ifdef P4_TO_P8
static void
p8est_iter_init_edge_from_face (p8est_iter_edge_args_t * args,
                                p4est_iter_face_args_t * face_args,
                                int dir, int side)
{
  const int           ntc_str = P4EST_CHILDREN / 2;
  int                 j, k;
  int                 c0, c1, *cc;
  p8est_iter_edge_info_t *info = &(args->info);
  p8est_iter_face_side_t *fside;
  p8est_iter_edge_side_t *eside;
  int                *num_to_child = face_args->num_to_child;
  int                *e_start_idx2;
  int                 limit = face_args->outside_face ? 1 : 2;
  int                 count = 0;
  int                 pos[2][2];
  sc_array_t         *common_corners = args->common_corners;

  pos[0][0] = 0;
  pos[0][1] = dir ? 2 : 1;
  pos[1][0] = dir ? 1 : 2;
  pos[1][1] = 3;

  info->p4est = face_args->info.p4est;
  info->ghost_layer = face_args->info.ghost_layer;
  info->tree_boundary = face_args->info.tree_boundary;
  sc_array_init (&(info->sides), sizeof (p8est_iter_edge_side_t));
  args->num_sides = limit * ntc_str / 2;
  sc_array_resize (&(info->sides), (size_t) args->num_sides);
  sc_array_init (&(common_corners[0]), sizeof (int));
  sc_array_init (&(common_corners[1]), sizeof (int));
  sc_array_resize (&(common_corners[0]), (size_t) args->num_sides);
  sc_array_resize (&(common_corners[1]), (size_t) args->num_sides);
  e_start_idx2 = args->start_idx2 = P4EST_ALLOC (int, args->num_sides);
  args->loop_args = face_args->loop_args;

  for (j = 0; j < 2; j++) {
    for (k = 0; k < limit; k++) {
      cc = (int *) sc_array_index_int (&(common_corners[0]), count);
      c0 = *cc = num_to_child[k * ntc_str + pos[1 - j][0]];
      cc = (int *) sc_array_index_int (&(common_corners[1]), count);
      c1 = *cc = num_to_child[k * ntc_str + pos[1 - j][1]];
      fside = p4est_iter_fside_array_index_int (&face_args->info.sides, k);
      eside = p8est_iter_eside_array_index_int (&(info->sides), count);
      eside->orientation = (int8_t) ((c0 < c1) ? 0 : 1);
      eside->treeid = fside->treeid;
      eside->edge = (int8_t) p8est_child_corner_edges[c0][c1];
      e_start_idx2[count++] = num_to_child[k * ntc_str + pos[j][side]];
      if (p8est_edge_faces[eside->edge][0] == fside->face) {
        eside->faces[0] = j;
        eside->faces[1] = 2 + k;
      }
      else {
        eside->faces[0] = 2 + k;
        eside->faces[1] = j;
      }
    }
  }

  /* also initialize the corner that is in the middle of the edge */
  args->remote = face_args->remote;
  if (args->loop_args->loop_corner) {
    p8est_iter_init_corner_from_edge (&(args->corner_args), args);
  }
}
#endif

#if 0                           /* currently unused */

static int
fside_compare (const void *a, const void *b)
{
  p4est_locidx_t      diff;
  const p4est_iter_face_side_t *A = (const p4est_iter_face_side_t *) a;
  const p4est_iter_face_side_t *B = (const p4est_iter_face_side_t *) b;

  diff = (A->treeid - B->treeid);
  if (diff) {
    return (int) diff;
  }
  return (A->face - B->face);
}

#endif /* 0 */

/* initialize face args for a face between trees */
static void
p4est_iter_init_face (p4est_iter_face_args_t * args, p4est_t * p4est,
                      p4est_ghost_t * ghost_layer,
                      p4est_iter_loop_args_t * loop_args, p4est_topidx_t t,
                      int f)
{
  const int           ntc_str = P4EST_CHILDREN / 2;
  int                 i;
  int                 c;
  int                 count = 0;
  p4est_iter_face_info_t *info = &(args->info);
  p4est_iter_face_side_t *fside;
  int                *num_to_child = args->num_to_child;
  int                *start_idx2 = args->start_idx2;
#ifdef P4_TO_P8
  int                 ref, set;
#endif
  p4est_connectivity_t *conn = p4est->connectivity;
  p4est_topidx_t      nt = conn->tree_to_tree[t * P4EST_FACES + f];
  int                 nf = (int) conn->tree_to_face[t * P4EST_FACES + f];
  int                 o = nf / P4EST_FACES;

  nf %= P4EST_FACES;

  if (nt < t || (nt == t && nf < f)) {
    p4est_iter_init_face (args, p4est, ghost_layer, loop_args, nt, nf);
    return;
  }

  args->loop_args = loop_args;
  info->p4est = p4est;
  info->ghost_layer = ghost_layer;
  info->tree_boundary = (int8_t) P4EST_CONNECT_FACE;
  sc_array_init (&(info->sides), sizeof (p4est_iter_face_side_t));

#ifdef P4_TO_P8
  ref = p8est_face_permutation_refs[f][nf];
  set = p8est_face_permutation_sets[ref][o];
#endif

  if (t == nt && nf == f) {
    nt = -1;
  }

  args->outside_face = (nt == -1);

  fside = (p4est_iter_face_side_t *) sc_array_push (&(info->sides));
  fside->face = (int8_t) f;
  fside->treeid = t;
  start_idx2[count++] = 0;
  o = info->orientation = 0;

  if (nt != -1) {
    fside = (p4est_iter_face_side_t *) sc_array_push (&(info->sides));
    fside->treeid = nt;
    fside->face = (int8_t) nf;
    start_idx2[count++] = 0;
    o = info->orientation =
      conn->tree_to_face[t * P4EST_FACES + f] / P4EST_FACES;
  }

  /* for each corner, find the touching corner on the other tree */
  for (i = 0; i < ntc_str; i++) {
    c = p4est_face_corners[f][i];
    num_to_child[i] = c;
    if (nt != -1) {
#ifndef P4_TO_P8
      num_to_child[ntc_str + i] = p4est_face_corners[nf][o == 0 ? i : 1 - i];
#else
      num_to_child[ntc_str + i] = p8est_face_corners[nf]
        [p8est_face_permutations[set][i]];
#endif
    }
  }

#ifdef P4_TO_P8
  if (loop_args->loop_edge) {
    p8est_iter_init_edge_from_face (&(args->edge_args[0][0]), args, 0, 0);
    p8est_iter_init_edge_from_face (&(args->edge_args[0][1]), args, 0, 1);
    p8est_iter_init_edge_from_face (&(args->edge_args[1][0]), args, 1, 0);
    p8est_iter_init_edge_from_face (&(args->edge_args[1][1]), args, 1, 1);
  }
#endif
  if (loop_args->loop_corner) {
    p4est_iter_init_corner_from_face (&(args->corner_args), args);
  }

  if (nt != -1) {
    p4est_iter_init_loop_face (loop_args, t, nt, p4est, ghost_layer);
  }
  else {
    p4est_iter_init_loop_outside_face (loop_args, t, p4est, ghost_layer);
  }
}

static void
p4est_iter_reset_face (p4est_iter_face_args_t * args)
{
  if (args->loop_args->loop_corner) {
    p4est_iter_reset_corner (&(args->corner_args));
  }
#ifdef P4_TO_P8
  if (args->loop_args->loop_edge) {
    p8est_iter_reset_edge (&(args->edge_args[0][0]));
    p8est_iter_reset_edge (&(args->edge_args[0][1]));
    p8est_iter_reset_edge (&(args->edge_args[1][0]));
    p8est_iter_reset_edge (&(args->edge_args[1][1]));
  }
#endif
  sc_array_reset (&(args->info.sides));
}

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

  p4est_iter_loop_args_t *loop_args = args->loop_args;
  int                 start_level = loop_args->level;
  int                *start_idx2 = args->start_idx2;
  int                *level_num = loop_args->level_num;
  sc_array_t        **quadrants = loop_args->quadrants;
  size_t            **zindex = loop_args->index;
  size_t             *first_index = loop_args->first_index;
  int                *num_to_child = args->num_to_child;
  p4est_quadrant_t  **test = loop_args->test;
  size_t             *count = loop_args->count;
  int                *test_level = loop_args->test_level;
  int                *quad_idx2 = loop_args->quad_idx2;
  int8_t             *refine = loop_args->refine;
  int                 limit;

  int                 i;
  int                *Level = &(loop_args->level);
  int                 side;
  int                 type;
  int                 st;
  int                 level_idx2;
  p4est_iter_face_info_t *info = &(args->info);
  p4est_iter_face_side_t *fside;
  p4est_quadrant_t  **quads;
  p4est_locidx_t     *quadids;
  int8_t             *is_ghost;
  int                 child_corner;
  int8_t              has_local;
  sc_array_t          test_view;
  p4est_iter_corner_args_t *corner_args = &(args->corner_args);
  sc_array_t         *tier_rings = loop_args->tier_rings;
#ifdef P4_TO_P8
  int                 dir;
#endif
  int                 stop_refine;
  int                 all_empty;

  /* if we are at an outside face, then there is no right half to our search
   * that needs to be coordinated with the left half */
  limit = args->outside_face ? left : right;

  /* level_idx2 moves us to the correct set of bounds within the index arrays
   * for the level: it is a set of bounds because it includes all children at
   * this level */
  level_idx2 = start_level * P4EST_ITER_STRIDE;

  for (side = left; side <= limit; side++) {

    /* start_idx2 gives the ancestor id at level for the search area on this
     * side, so quad_idx2[side] now gives the correct location in
     * zindex[sidetype] of the bounds of the search area */
    quad_idx2[side] = level_idx2 + start_idx2[side];

    /* get the location in quadrants[sidetype] of the first quadrant in the
     * search area, and the count of quadrants in the search area */
    for (type = local; type <= ghost; type++) {
      st = side * 2 + type;
      first_index[st] = zindex[st][quad_idx2[side]];
      count[st] = (zindex[st][quad_idx2[side] + 1] - first_index[st]);
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
        st = side * 2 + type;
        if (count[st]) {
          test[st] = p4est_quadrant_array_index (quadrants[st],
                                                 first_index[st]);
          test_level[st] = (int) test[st]->level;
        }
        else {
          test[st] = NULL;
          test_level[st] = -1;
        }
      }
    }
    /* initially assume that each side needs to be refined */
    refine[left] = refine[right] = 1;
    stop_refine = 0;
    has_local = 0;

    /* get a candidate from each sidetype */
    for (side = left; side <= limit; side++) {
      for (type = local; type <= ghost; type++) {
        st = side * 2 + type;
        /* if the candidate from sidetype is the same size as the search area,
         * then we do not refine this side */
        if (test_level[st] == *Level) {
          if (!stop_refine) {
            stop_refine = 1;
            /* we are not going to recur on the next level, instead moving to
             * the next branch on this level */
            level_num[*Level]++;
            /* if there is no face callback (i.e., we are just running
             * face_iterate to find edges and corners), then we're done with
             * this branch */
            if (iter_face == NULL) {
              goto change_search_area;
            }
          }
          P4EST_ASSERT (count[st] == 1);
          P4EST_ASSERT (count[side * 2 + (type ^ 1)] == 0);
          refine[side] = 0;
          fside = p4est_iter_fside_array_index_int (&(info->sides), side);
          fside->is_hanging = 0;
          fside->is.full.quad = test[st];
          fside->is.full.quadid = (p4est_locidx_t) first_index[st];
          has_local = (has_local || (type == local));
          fside->is.full.is_ghost = (type == ghost);
        }
      }
    }
    for (side = left; side <= limit; side++) {
      if (refine[side]) {
        if (stop_refine && count[side * 2 + local] == 0 &&
            count[side * 2 + ghost] == 0) {
          fside = p4est_iter_fside_array_index_int (&info->sides, side);
          fside->is_hanging = 0;
          fside->is.full.quad = NULL;
          fside->is.full.is_ghost = 1;
          fside->is.full.quadid = -1;
          refine[side] = 0;
        }
      }
      if (refine[side]) {
        quad_idx2[side] = level_idx2 + P4EST_ITER_STRIDE;
        for (type = local; type <= ghost; type++) {
          st = side * 2 + type;
          sc_array_init_view (&test_view, quadrants[st],
                              first_index[st], count[st]);
          p4est_iter_tier_insert (&test_view, *Level, zindex[st] +
                                  quad_idx2[side], first_index[st],
                                  tier_rings, test[st]);
        }
        if (stop_refine) {
          fside = p4est_iter_fside_array_index_int (&(info->sides), side);
          fside->is_hanging = 1;
          quads = fside->is.hanging.quad;
          quadids = fside->is.hanging.quadid;
          is_ghost = fside->is.hanging.is_ghost;
          for (i = 0; i < P4EST_CHILDREN / 2; i++) {
            /* fside expects the hanging quadrants listed in z order, not search
             * order */
            child_corner = num_to_child[side * ntc_str + i];
            child_corner =
              p4est_corner_face_corners[child_corner][fside->face];
            quads[child_corner] = NULL;
            quadids[child_corner] = -1;
            is_ghost[child_corner] = 1;
            quad_idx2[side] = level_idx2 + P4EST_ITER_STRIDE +
              num_to_child[side * ntc_str + i];
            for (type = local; type <= ghost; type++) {
              st = side * 2 + type;
              first_index[st] = zindex[st][quad_idx2[side]];
              count[st] = zindex[st][quad_idx2[side] + 1] - first_index[st];
              /* if the search area is non-empty, by the two to one condition
               * it must contain exactly one quadrant: if one of the two types
               * is local, we run iter_face */
              if (count[st]) {
                quads[child_corner] = p4est_quadrant_array_index
                  (quadrants[st], first_index[st]);
                P4EST_ASSERT ((int) quads[child_corner]->level == *Level + 1);
                quadids[child_corner] = (p4est_locidx_t) first_index[st];
                is_ghost[child_corner] = (type == ghost);
                has_local = (has_local || (type == local));
              }
            }
          }
        }
      }
    }
    if (stop_refine) {
      if (has_local) {
        iter_face (info, user_data);
      }
    }
    else {
      /* if we refined both sides, we descend to the next level from this branch
       * and continue searching there */
      level_num[++(*Level)] = 0;
      level_idx2 += P4EST_ITER_STRIDE;
    }

  change_search_area:

    for (;;) {
      /* if we tried to advance the search area on start_level, we've completed
       * the search */
      if (level_num[start_level] > 0) {
        P4EST_ASSERT (*Level == start_level);
        return;
      }

      /* if we have tried to advance the search area past the number of
       * descendants, that means that we have completed all of the branches on
       * this level */
      if (level_num[*Level] == P4EST_CHILDREN / 2) {
#ifdef P4_TO_P8
        /* if we have an edge callback, we need to run it on all of the edges
         * between the face branches on this level */
        if (loop_args->loop_edge) {
          for (dir = 0; dir < 2; dir++) {
            for (side = 0; side < 2; side++) {
              P4EST_ASSERT (args->edge_args[dir][side].num_sides ==
                            2 * (limit + 1));
              p4est_iter_copy_indices (loop_args,
                                       args->edge_args[dir][side].start_idx2,
                                       limit + 1, 2);
              p8est_edge_iterate (&(args->edge_args[dir][side]), user_data,
                                  iter_edge, iter_corner);
            }
          }
        }
#endif
        /* if we have a corner callback, we need to run it on the corner between
         * the face branches on this level */
        if (iter_corner != NULL) {
          P4EST_ASSERT (corner_args->num_sides ==
                        (P4EST_CHILDREN / 2) * (limit + 1));
          p4est_iter_copy_indices (loop_args, corner_args->start_idx2,
                                   limit + 1, P4EST_HALF);
          p4est_corner_iterate (corner_args, user_data, iter_corner);
        }

        /* now that we're done on this level, go up a level and over a branch */
        level_num[--(*Level)]++;
        level_idx2 -= P4EST_ITER_STRIDE;
      }
      else {
        /* at this point, we need to initialize the bounds of the search areas
         * for this new branch */
        all_empty = 1;
        for (side = left; side <= limit; side++) {
          quad_idx2[side] =
            level_idx2 + num_to_child[side * ntc_str + level_num[*Level]];
        }
        for (side = left; side <= limit; side++) {
          for (type = local; type <= ghost; type++) {
            st = side * 2 + type;
            first_index[st] = zindex[st][quad_idx2[side]];
            count[st] = (zindex[st][quad_idx2[side] + 1] - first_index[st]);
            if (type == local && count[st]) {
              all_empty = 0;
            }
          }
        }
        if (all_empty) {
          /* if there are no local quadrants in either of the search areas, we're
           * done with this search area and proceed to the next branch on this
           * level */
          level_num[*Level]++;
        }
        else {
          /* otherwise we are done changing the search area */
          break;
        }
      }
    }
  }
}

/* volume iterate functions */
typedef struct p4est_iter_volume_args
{
  p4est_iter_loop_args_t *loop_args;
  int                 start_idx2;
  p4est_iter_face_args_t face_args[P4EST_DIM][P4EST_CHILDREN / 2];
#ifdef P4_TO_P8
  p8est_iter_edge_args_t edge_args[P4EST_DIM][2];
#endif
  p4est_iter_corner_args_t corner_args;
  p4est_iter_volume_info_t info;
  int                 remote;
}
p4est_iter_volume_args_t;

/* given valid volume arguments, setup face arguments for a face search that is
 * called for a face between two adjacent volumes: there are P4EST_DIM
 * directions the face can be oriented, and each direction can be run in a
 * different position, based on the child_ids of the two volumes on either side
 * of it.
 */
static void
p4est_iter_init_face_from_volume (p4est_iter_face_args_t * args,
                                  p4est_iter_volume_args_t * volume_args,
                                  int dir, int pos)
{
  const int           ntc_str = P4EST_CHILDREN / 2;
  int                 i, j;
  p4est_iter_face_info_t *info = &(args->info);
  p4est_iter_face_side_t *fside;

  info->p4est = volume_args->info.p4est;
  info->ghost_layer = volume_args->info.ghost_layer;
  info->orientation = 0;
  info->tree_boundary = 0;
  sc_array_init (&(info->sides), sizeof (p4est_iter_face_side_t));
  sc_array_resize (&(info->sides), 2);
  args->loop_args = volume_args->loop_args;
  args->outside_face = 0;

  args->start_idx2[0] = p4est_face_corners[dir * 2][pos];
  args->start_idx2[1] = p4est_face_corners[dir * 2 + 1][pos];

  for (i = 0; i < 2; i++) {
    for (j = 0; j < ntc_str; j++) {
      args->num_to_child[i * ntc_str + j] =
        p4est_face_corners[dir * 2 + (1 - i)][j];
    }
  }

  fside = p4est_iter_fside_array_index_int (&info->sides, 0);
  fside->treeid = volume_args->info.treeid;
  fside->face = (int8_t) (2 * dir + 1);
  fside = p4est_iter_fside_array_index_int (&info->sides, 1);
  fside->treeid = volume_args->info.treeid;
  fside->face = (int8_t) 2 *dir;

  args->remote = volume_args->remote;
#ifdef P4_TO_P8
  if (args->loop_args->loop_edge) {
    p8est_iter_init_edge_from_face (&(args->edge_args[0][0]), args, 0, 0);
    p8est_iter_init_edge_from_face (&(args->edge_args[0][1]), args, 0, 1);
    p8est_iter_init_edge_from_face (&(args->edge_args[1][0]), args, 1, 0);
    p8est_iter_init_edge_from_face (&(args->edge_args[1][1]), args, 1, 1);
  }
#endif
  if (args->loop_args->loop_corner) {
    p4est_iter_init_corner_from_face (&(args->corner_args), args);
  }
}

/* given valid volume arguments, setup edge arguments for an edge search that
 * is called for an edge between four adjacent volumes: there are P4EST_DIM
 * directions the edge can be oriented, and each direction can be run in oen
 * of two positioins, based on the child_ids of the four volumes surrounding
 * it.
 */
#ifdef P4_TO_P8
static void
p8est_iter_init_edge_from_volume (p8est_iter_edge_args_t * args,
                                  p4est_iter_volume_args_t * volume_args,
                                  int dir, int side)
{
  int                 i;
  int                *cc;
  p8est_iter_edge_info_t *info = &(args->info);
  p8est_iter_edge_side_t *eside;
  sc_array_t         *common_corners = args->common_corners;

  info->p4est = volume_args->info.p4est;
  info->ghost_layer = volume_args->info.ghost_layer;
  info->tree_boundary = 0;
  sc_array_init (&(info->sides), sizeof (p8est_iter_edge_side_t));
  sc_array_resize (&(info->sides), 4);
  sc_array_init (&(common_corners[0]), sizeof (int));
  sc_array_init (&(common_corners[1]), sizeof (int));
  sc_array_resize (&(common_corners[0]), 4);
  sc_array_resize (&(common_corners[1]), 4);
  args->start_idx2 = P4EST_ALLOC (int, 4);
  args->loop_args = volume_args->loop_args;
  args->num_sides = 4;

  for (i = 0; i < 4; i++) {
    args->start_idx2[i] = p8est_face_corners[dir * 2 + side][i];
    cc = (int *) sc_array_index_int (&(common_corners[0]), i);
    *cc = p8est_face_corners[dir * 2][3 - i];
    cc = (int *) sc_array_index_int (&(common_corners[1]), i);
    *cc = p8est_face_corners[dir * 2 + 1][3 - i];
    eside = p8est_iter_eside_array_index_int (&info->sides, i);
    eside->treeid = volume_args->info.treeid;
    eside->orientation = 0;
    eside->edge = (int8_t) (4 * dir + (3 - i));
    eside->faces[0] = (i >> 1);
    eside->faces[1] = 2 + (i & 1);
  }

  args->remote = volume_args->remote;
  if (args->loop_args->loop_corner) {
    p8est_iter_init_corner_from_edge (&(args->corner_args), args);
  }
}
#endif

/* given valid volume arguments, setup corner arguments for a corner search
 * that is called for a corner between P4EST_CHILDREN adjacent volumes.
 */
static void
p4est_iter_init_corner_from_volume (p4est_iter_corner_args_t * args,
                                    p4est_iter_volume_args_t * volume_args)
{
  int                 i;
  p4est_iter_corner_info_t *info = &(args->info);
  p4est_iter_corner_side_t *cside;

  info->p4est = volume_args->info.p4est;
  info->ghost_layer = volume_args->info.ghost_layer;
  info->tree_boundary = 0;
  sc_array_init (&(info->sides), sizeof (p4est_iter_corner_side_t));
  sc_array_resize (&(info->sides), P4EST_CHILDREN);
  args->start_idx2 = P4EST_ALLOC (int, P4EST_CHILDREN);
  args->num_sides = P4EST_CHILDREN;
  args->loop_args = volume_args->loop_args;

  for (i = 0; i < P4EST_CHILDREN; i++) {
    args->start_idx2[i] = i;
    cside = p4est_iter_cside_array_index_int (&(info->sides), i);
    cside->treeid = volume_args->info.treeid;
    cside->corner = (int8_t) (P4EST_CHILDREN - 1 - i);
#ifndef P4_TO_P8
    cside->faces[0] = (i >> 1);
    cside->faces[1] = 2 + (i & 1);
#else
    cside->faces[0] = (i >> 1);
    cside->faces[1] = 4 + ((i & 4) >> 1) + (i & 1);
    cside->faces[2] = 8 + (i & 3);
    cside->edges[0] = (i & 1);
    cside->edges[1] = 2 + ((i & 2) >> 1);
    cside->edges[2] = 4 + (i >> 2);
#endif
  }
  args->remote = volume_args->remote;
}

/* initialize volume arguments for a search in a tree */
static void
p4est_iter_init_volume (p4est_iter_volume_args_t * args, p4est_t * p4est,
                        p4est_ghost_t * ghost_layer,
                        p4est_iter_loop_args_t * loop_args, p4est_topidx_t t)
{
  int                 i, j;

  args->loop_args = loop_args;
  args->info.p4est = p4est;
  args->info.ghost_layer = ghost_layer;
  args->info.treeid = t;
  args->start_idx2 = 0;

  for (i = 0; i < P4EST_DIM; i++) {
    for (j = 0; j < P4EST_CHILDREN / 2; j++) {
      p4est_iter_init_face_from_volume (&(args->face_args[i][j]), args, i, j);
    }
#ifdef P4_TO_P8
    if (loop_args->loop_edge) {
      for (j = 0; j < 2; j++) {
        p8est_iter_init_edge_from_volume (&(args->edge_args[i][j]), args, i,
                                          j);
      }
    }
#endif
  }
  if (loop_args->loop_corner) {
    p4est_iter_init_corner_from_volume (&(args->corner_args), args);
  }

  p4est_iter_init_loop_volume (args->loop_args, t, p4est, ghost_layer);
}

static void
p4est_iter_reset_volume (p4est_iter_volume_args_t * args)
{
  int                 i, j;

  for (i = 0; i < P4EST_DIM; i++) {
    for (j = 0; j < P4EST_CHILDREN / 2; j++) {
      p4est_iter_reset_face (&(args->face_args[i][j]));
    }
#ifdef P4_TO_P8
    if (args->loop_args->loop_edge) {
      for (j = 0; j < 2; j++) {
        p8est_iter_reset_edge (&(args->edge_args[i][j]));
      }
    }
#endif
  }
  if (args->loop_args->loop_corner) {
    p4est_iter_reset_corner (&(args->corner_args));
  }
}

/* when there is only a volume callback, there is no coordination necessary, so
 * a simple loop is performed */
static void
p4est_volume_iterate_simple (p4est_t * p4est, p4est_ghost_t * ghost_layer,
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
    tree = p4est_tree_array_index (trees, t);
    quadrants = &(tree->quadrants);
    n_quads = quadrants->elem_count;
    for (si = 0; si < n_quads; si++) {
      info.quad = p4est_quadrant_array_index (quadrants, si);
      info.quadid = si;
      iter_volume (&info, user_data);
    }
  }
}

static void
p4est_volume_iterate (p4est_iter_volume_args_t * args, void *user_data,
                      p4est_iter_volume_t iter_volume,
                      p4est_iter_face_t iter_face,
#ifdef P4_TO_P8
                      p8est_iter_edge_t iter_edge,
#endif
                      p4est_iter_corner_t iter_corner)
{
  const int           local = 0;
  const int           ghost = 1;

  int                 dir, side, type;

  p4est_iter_loop_args_t *loop_args = args->loop_args;
  int                 start_level = loop_args->level;
  int                *Level = &(loop_args->level);
  int                 start_idx2 = args->start_idx2;
  int                *level_num = loop_args->level_num;
  sc_array_t        **quadrants = loop_args->quadrants;
  size_t            **zindex = loop_args->index;
  size_t             *first_index = loop_args->first_index;
  p4est_quadrant_t  **test = loop_args->test;
  size_t             *count = loop_args->count;
  int                *test_level = loop_args->test_level;
  sc_array_t         *tier_rings = loop_args->tier_rings;
  int                 quad_idx2;
  sc_array_t          test_view;
  p4est_iter_volume_info_t *info = &(args->info);
  int                 level_idx2;
  int                 refine;

  /* level_idx2 moves us to the correct set of bounds within the index arrays
   * for the level: it is a set of bounds because it includes all children at
   * this level */
  level_idx2 = start_level * P4EST_ITER_STRIDE;

  /* start_idx2 gives the ancestor id at level for the search area,
   * so quad_idx2 now gives the correct location in
   * index[type] of the bounds of the search area */
  quad_idx2 = level_idx2 + start_idx2;
  for (type = local; type <= ghost; type++) {
    first_index[type] = zindex[type][quad_idx2];
    count[type] = zindex[type][quad_idx2 + 1] - first_index[type];
  }

  /* if ther are no local quadrants, nothing to be done */
  if (!count[local]) {
    return;
  }

  /* we think of the search tree as being rooted at start_level, so we can
   * think the branch number at start_level as 0, even if it actually is not */
  level_num[start_level] = 0;

  for (;;) {

    refine = 1;
    /* for each type, get the first quadrant in the search area */
    for (type = local; type <= ghost; type++) {
      if (count[type]) {
        test[type] = p4est_quadrant_array_index (quadrants[type],
                                                 first_index[type]);
        test_level[type] = (int) test[type]->level;
        /* if the quadrant is the same size as the search area, we're done
         * search */
        if (test_level[type] == *Level) {
          refine = 0;
          P4EST_ASSERT (!count[type ^ 1]);
          /* if the quadrant is local, we run the callback */
          if (type == local) {
            info->quad = test[type];
            info->quadid = (p4est_locidx_t) first_index[type];
            if (iter_volume != NULL) {
              iter_volume (info, user_data);
            }
          }
          /* proceed to the next search area on this level */
          level_num[*Level]++;
        }
      }
      else {
        test[type] = NULL;
        test_level[type] = -1;
      }
    }

    if (refine) {
      /* we need to refine, we take the search area and split it up, taking the
       * indices for the refined search areas and placing them on the next tier in
       * index[type] */
      quad_idx2 = level_idx2 + P4EST_ITER_STRIDE;
      for (type = local; type <= ghost; type++) {
        sc_array_init_view (&test_view, quadrants[type],
                            first_index[type], count[type]);
        p4est_iter_tier_insert (&test_view, *Level, zindex[type] + quad_idx2,
                                first_index[type], tier_rings, test[type]);
      }

      /* we descend to the first descendant search area and search there */
      level_num[++(*Level)] = 0;
      level_idx2 += P4EST_ITER_STRIDE;
    }

    for (;;) {
      /* if we tried to advance the search area on start_level, we've completed
       * the search */
      if (level_num[start_level] > 0) {
        return;
      }

      /* if we have tried to advance the search area past the number of
       * descendants, that means that we have completed all of the branches on
       * this level. we can now run the face_iterate for all of the faces between
       * search areas on the level*/
      if (level_num[*Level] == P4EST_CHILDREN) {
        /* for each direction */
        for (dir = 0; dir < P4EST_DIM; dir++) {
          for (side = 0; side < P4EST_CHILDREN / 2; side++) {
            p4est_iter_copy_indices (loop_args,
                                     args->face_args[dir][side].start_idx2,
                                     1, 2);
            p4est_face_iterate (&(args->face_args[dir][side]), user_data,
                                iter_face,
#ifdef P4_TO_P8
                                iter_edge,
#endif
                                iter_corner);
          }
        }
#ifdef P4_TO_P8
        /* if there is an edge or a corner callback, we need to use
         * edge_iterate, so we set up the common corners and edge ids
         * for all of the edges between the search areas */
        if (loop_args->loop_edge) {
          for (dir = 0; dir < P4EST_DIM; dir++) {
            for (side = 0; side < 2; side++) {
              p4est_iter_copy_indices (loop_args,
                                       args->edge_args[dir][side].start_idx2,
                                       1, 4);
              p8est_edge_iterate (&(args->edge_args[dir][side]), user_data,
                                  iter_edge, iter_corner);
            }
          }
        }
#endif
        /* if there is a corner callback, we need to call corner_iterate on
         * the corner in the middle of the search areas */
        if (loop_args->loop_corner) {
          p4est_iter_copy_indices (loop_args, args->corner_args.start_idx2, 1,
                                   P4EST_CHILDREN);
          p4est_corner_iterate (&(args->corner_args), user_data, iter_corner);
        }
        /* we are done at the level, so we go up a level and over a branch */
        level_num[--(*Level)]++;
        level_idx2 -= P4EST_ITER_STRIDE;
      }
      else {
        /* quad_idx now gives the location in index[type] of the bounds
         * of the current search area, from which we get the first quad
         * and the count */
        quad_idx2 = level_idx2 + level_num[*Level];
        for (type = local; type <= ghost; type++) {
          first_index[type] = zindex[type][quad_idx2];
          count[type] = zindex[type][quad_idx2 + 1] - first_index[type];
        }
        if (!count[local]) {
          /* if there are no local quadrants, we are done with this search area,
           * and we advance to the next branch at this level */
          level_num[*Level]++;
        }
        else {
          /* otherwise we are done changing the search area */
          break;
        }
      }
    }
  }
  P4EST_ASSERT (*Level == start_level);
}

static int32_t     *
p4est_iter_get_boundaries (p4est_t * p4est, p4est_topidx_t * last_run_tree,
                           int remote)
{
  p4est_topidx_t      ti;
  int                 i;
  int                 rank = p4est->mpirank;
  p4est_connectivity_t *conn = p4est->connectivity;
  sc_array_t         *trees = p4est->trees;
  size_t              global_num_trees = trees->elem_count;
  int32_t            *init = P4EST_ALLOC_ZERO (int32_t, global_num_trees);
  int32_t            *owned = P4EST_ALLOC_ZERO (int32_t, global_num_trees);
  int32_t             touch;
  int32_t             mask;
  p4est_topidx_t      t, nt, ot;
  p4est_topidx_t      first_local_tree = p4est->first_local_tree;
  p4est_topidx_t      last_local_tree = p4est->last_local_tree;
  p4est_quadrant_t   *lq = &(p4est->global_first_position[rank]);
  p4est_quadrant_t    temp;
  p4est_quadrant_t   *uq = &(p4est->global_first_position[rank + 1]);
  uint64_t            uqid;
  p4est_quadrant_t   *tlq, *tuq;
  int                 f, nf, c, c2, nc, oc;
  p4est_topidx_t      corner;
  p4est_topidx_t     *ttc = conn->tree_to_corner;
  p4est_topidx_t     *ctt_offset = conn->ctt_offset;
  p4est_topidx_t     *ctt = conn->corner_to_tree;
  int8_t             *ctc = conn->corner_to_corner;
#ifdef P4_TO_P8
  int                 e, ne, oe;
  int                 nc2;
  p4est_topidx_t      edge;
  p4est_topidx_t     *tte = conn->tree_to_edge;
  p4est_topidx_t     *ett_offset = conn->ett_offset;
  p4est_topidx_t     *ett = conn->edge_to_tree;
  int8_t             *ete = conn->edge_to_edge;
  int                 ref, set;
  int                 this_o;
#endif
  p4est_topidx_t     *ttt = conn->tree_to_tree;
  int8_t             *ttf = conn->tree_to_face;
#ifndef P4_TO_P8
  int                 corner_offset = 4;
#else
  int                 edge_offset = 6;
  int                 corner_offset = 18;
#endif
  int                 o;

  *last_run_tree = -1;

  if (uq->p.which_tree > last_local_tree) {
    uq = NULL;
  }
  else {
    P4EST_ASSERT (uq->p.which_tree == last_local_tree);
    uqid = p4est_quadrant_linear_id (uq, P4EST_QMAXLEVEL);
    p4est_quadrant_set_morton (&temp, P4EST_QMAXLEVEL, uqid - 1);
    uq = &temp;
  }

  for (t = first_local_tree; t <= last_local_tree; t++) {
    if (!remote) {
      tlq = (t == first_local_tree) ? lq : NULL;
      tuq = (t == last_local_tree) ? uq : NULL;
    }
    else {
      tlq = NULL;
      tuq = NULL;
    }
    touch = p4est_find_range_boundaries (tlq, tuq, 0, NULL,
#ifdef P4_TO_P8
                                         NULL,
#endif
                                         NULL);
    if (!touch) {
      continue;
    }
    mask = 0x00000001;
    for (f = 0; f < P4EST_FACES; f++, mask <<= 1) {
      if ((touch & mask) && !(init[t] & mask)) {
        nt = ttt[t * P4EST_FACES + f];
        nf = (int) ttf[t * P4EST_FACES + f];
        nf %= P4EST_FACES;
        init[t] |= mask;
        init[nt] |= (((int32_t) 1) << nf);
        if (t > nt || ((t == nt) && (f >= nf))) {
          owned[t] |= mask;
          if (t > *last_run_tree) {
            *last_run_tree = t;
          }
        }
        else {
          owned[nt] |= (((int32_t) 1) << nf);
          if (nt > *last_run_tree) {
            *last_run_tree = nt;
          }
        }
      }
    }
#ifdef P4_TO_P8
    for (e = 0; e < 12; e++, mask <<= 1) {
      if ((touch & mask) && !(init[t] & mask)) {
        if (tte != NULL) {
          edge = tte[t * 12 + e];
        }
        else {
          edge = -1;
        }
        if (edge >= 0) {
          ot = -1;
          oe = -1;
          for (ti = ett_offset[edge]; ti < ett_offset[edge + 1]; ti++) {
            nt = ett[ti];
            ne = (int) ete[ti];
            ne %= 12;
            init[nt] |= (((int32_t) 1) << (ne + edge_offset));
            if (nt > ot || ((nt == ot) && (ne >= oe))) {
              ot = nt;
              oe = ne;
            }
          }
        }
        else {
          ot = t;
          oe = e;
          init[t] |= mask;
          for (i = 0; i < 2; i++) {
            c = p8est_edge_corners[e][0];
            c2 = p8est_edge_corners[e][1];
            f = p8est_edge_faces[e][i];
            c = p8est_corner_face_corners[c][f];
            c2 = p8est_corner_face_corners[c2][f];
            nt = ttt[t * P4EST_FACES + f];
            nf = (int) ttf[t * P4EST_FACES + f];
            o = nf / P4EST_FACES;
            nf %= P4EST_FACES;
            if (t == nt && f == nf) {
              continue;
            }
            ref = p8est_face_permutation_refs[f][nf];
            set = p8est_face_permutation_sets[ref][o];
            nc = p8est_face_permutations[set][c];
            nc2 = p8est_face_permutations[set][c2];
            nc = p8est_face_corners[nf][nc];
            nc2 = p8est_face_corners[nf][nc2];
            ne = p8est_child_corner_edges[nc][nc2];
            init[nt] |= (((int32_t) 1) << (ne + edge_offset));
            if (nt > ot || ((nt == ot) && (ne >= oe))) {
              ot = nt;
              oe = ne;
            }
          }
        }
        owned[ot] |= (((int32_t) 1) << (oe + edge_offset));
        if (ot > *last_run_tree) {
          *last_run_tree = ot;
        }
      }
    }
#endif
    for (c = 0; c < P4EST_CHILDREN; c++, mask <<= 1) {
      if ((touch & mask) && !(init[t] & mask)) {
        if (ttc != NULL) {
          corner = ttc[t * P4EST_CHILDREN + c];
        }
        else {
          corner = -1;
        }
        if (corner >= 0) {
          ot = -1;
          oc = -1;
          for (ti = ctt_offset[corner]; ti < ctt_offset[corner + 1]; ti++) {
            nt = ctt[ti];
            nc = (int) ctc[ti];
            init[nt] |= (((int32_t) 1) << (nc + corner_offset));
            if (nt > ot || ((nt == ot) && (nc >= oc))) {
              ot = nt;
              oc = nc;
            }
          }
        }
        else {
          ot = t;
          oc = c;
          init[t] |= mask;
          for (i = 0; i < P4EST_DIM; i++) {
            f = p4est_corner_faces[c][i];
            c2 = p4est_corner_face_corners[c][f];
            nt = ttt[t * P4EST_FACES + f];
            nf = (int) ttf[t * P4EST_FACES + f];
            o = nf / P4EST_FACES;
            nf %= P4EST_FACES;
            if (t == nt && f == nf) {
              continue;
            }
#ifndef P4_TO_P8
            nc = p4est_face_corners[nf][(o == 0) ? c2 : 1 - c2];
#else
            ref = p8est_face_permutation_refs[f][nf];
            set = p8est_face_permutation_sets[ref][o];
            nc2 = p8est_face_permutations[set][c2];
            nc = p8est_face_corners[nf][nc2];
#endif
            init[nt] |= (((int32_t) 1) << (nc + corner_offset));
            if (nt > ot || ((nt == ot) && (nc >= oc))) {
              ot = nt;
              oc = nc;
            }
          }
#ifdef P4_TO_P8
          for (i = 0; i < 3; i++) {
            e = p8est_corner_edges[c][i];
            c2 = (p8est_edge_corners[e][0] == c) ? 0 : 1;
            if (tte != NULL) {
              edge = tte[t * 12 + e];
            }
            else {
              edge = -1;
            }
            if (edge >= 0) {
              for (ti = ett_offset[edge]; ti < ett_offset[edge + 1]; ti++) {
                nt = ett[ti];
                ne = (int) ete[ti];
                this_o = ne / 12;
                ne %= 12;
                if (nt == t && ne == e) {
                  break;
                }
              }
              P4EST_ASSERT (ti < ett_offset[edge + 1]);
              for (ti = ett_offset[edge]; ti < ett_offset[edge + 1]; ti++) {
                nt = ett[ti];
                ne = (int) ete[ti];
                o = ne / 12;
                ne %= 12;
                nc = p8est_edge_corners[ne][(this_o == o) ? c2 : 1 - c2];
                init[nt] |= (((int32_t) 1) << (nc + corner_offset));
                if (nt > ot || ((nt == ot) && (nc >= oc))) {
                  ot = nt;
                  oc = nc;
                }
              }
            }
          }
#endif
        }
        owned[ot] |= (((int32_t) 1) << (oc + corner_offset));
        if (ot > *last_run_tree) {
          *last_run_tree = ot;
        }
      }
    }
  }

  P4EST_FREE (init);
  return owned;
}

void
p4est_iterate_ext (p4est_t * p4est, p4est_ghost_t * Ghost_layer,
                   void *user_data, p4est_iter_volume_t iter_volume,
                   p4est_iter_face_t iter_face,
#ifdef P4_TO_P8
                   p8est_iter_edge_t iter_edge,
#endif
                   p4est_iter_corner_t iter_corner, int remote)
{
  int                 f, c;
  p4est_topidx_t      t;
  p4est_ghost_t       empty_ghost_layer;
  p4est_ghost_t      *ghost_layer;
  sc_array_t         *trees = p4est->trees;
  p4est_connectivity_t *conn = p4est->connectivity;
  size_t              global_num_trees = trees->elem_count;
  p4est_iter_loop_args_t *loop_args;
  p4est_iter_face_args_t face_args;
#ifdef P4_TO_P8
  int                 e;
  p8est_iter_edge_args_t edge_args;
#endif
  p4est_iter_corner_args_t corner_args;
  p4est_iter_volume_args_t args;
  p4est_topidx_t      first_local_tree = p4est->first_local_tree;
  p4est_topidx_t      last_local_tree = p4est->last_local_tree;
  p4est_topidx_t      last_run_tree;
  int32_t            *owned;
  int32_t             mask, touch;

  P4EST_ASSERT (p4est_is_valid (p4est));

  if (p4est->first_local_tree < 0 ||
      (iter_face == NULL && iter_corner == NULL &&
#ifdef P4_TO_P8
       iter_edge == NULL &&
#endif
       iter_volume == NULL)) {
    return;
  }

  if (Ghost_layer == NULL) {
    sc_array_init (&(empty_ghost_layer.ghosts), sizeof (p4est_quadrant_t));
    empty_ghost_layer.tree_offsets = P4EST_ALLOC_ZERO (p4est_locidx_t,
                                                       global_num_trees + 1);
    empty_ghost_layer.proc_offsets = P4EST_ALLOC_ZERO (p4est_locidx_t,
                                                       p4est->mpisize + 1);
    ghost_layer = &empty_ghost_layer;
  }
  else {
    ghost_layer = Ghost_layer;
  }

  /* simple loop if there is only a volume callback */
  if (iter_face == NULL && iter_corner == NULL
#ifdef P4_TO_P8
      && iter_edge == NULL
#endif
    ) {
    p4est_volume_iterate_simple (p4est, ghost_layer, user_data, iter_volume);
    if (Ghost_layer == NULL) {
      P4EST_FREE (empty_ghost_layer.tree_offsets);
      P4EST_FREE (empty_ghost_layer.proc_offsets);
    }
    return;
  }

  /** initialize arrays that keep track of where we are in the search */
  loop_args = p4est_iter_loop_args_new (conn,
#ifdef P4_TO_P8
                                        iter_edge,
#endif
                                        iter_corner, ghost_layer,
                                        p4est->mpisize);

  owned = p4est_iter_get_boundaries (p4est, &last_run_tree, remote);
  last_run_tree = (last_run_tree < last_local_tree) ? last_local_tree :
    last_run_tree;

  /* start with the assumption that we only run on entities touches by the
   * local processor's domain */
  args.remote = remote;
  face_args.remote = remote;
#ifdef P4_TO_P8
  edge_args.remote = remote;
#endif
  corner_args.remote = remote;

  /** we have to loop over all trees and not just local trees because of the
   * ghost layer */
  for (t = first_local_tree; t <= last_run_tree; t++) {
    if (t >= first_local_tree && t <= last_local_tree) {
      p4est_iter_init_volume (&args, p4est, ghost_layer, loop_args, t);

      p4est_volume_iterate (&args, user_data, iter_volume, iter_face,
#ifdef P4_TO_P8
                            iter_edge,
#endif
                            iter_corner);

      p4est_iter_reset_volume (&args);
    }

    touch = owned[t];
    if (!touch) {
      continue;
    }
    mask = 0x00000001;
    /* Now we need to run face_iterate on the faces between trees */
    for (f = 0; f < 2 * P4EST_DIM; f++, mask <<= 1) {
      if ((touch & mask) == 0) {
        continue;
      }
      p4est_iter_init_face (&face_args, p4est, ghost_layer, loop_args, t, f);
      p4est_face_iterate (&face_args, user_data, iter_face,
#ifdef P4_TO_P8
                          iter_edge,
#endif
                          iter_corner);
      p4est_iter_reset_face (&face_args);
    }

    /* if there is an edge or a corner callback, we need to run
     * edge_iterate on the edges between trees */
#ifdef P4_TO_P8
    if (loop_args->loop_edge) {
      for (e = 0; e < 12; e++, mask <<= 1) {
        if ((touch & mask) == 0) {
          continue;
        }
        p8est_iter_init_edge (&edge_args, p4est, ghost_layer, loop_args, t,
                              e);
        p8est_edge_iterate (&edge_args, user_data, iter_edge, iter_corner);
        p8est_iter_reset_edge (&edge_args);
      }
    }
    else {
      mask <<= 12;
    }
#endif

    if (loop_args->loop_corner) {
      for (c = 0; c < P4EST_CHILDREN; c++, mask <<= 1) {
        if ((touch & mask) == 0) {
          continue;
        }
        p4est_iter_init_corner (&corner_args, p4est, ghost_layer, loop_args,
                                t, c);
        p4est_corner_iterate (&corner_args, user_data, iter_corner);
        p4est_iter_reset_corner (&corner_args);
      }
    }

  }

  if (Ghost_layer == NULL) {
    P4EST_FREE (empty_ghost_layer.tree_offsets);
    P4EST_FREE (empty_ghost_layer.proc_offsets);
  }

  P4EST_FREE (owned);
  p4est_iter_loop_args_destroy (loop_args);
}

void
p4est_iterate (p4est_t * p4est, p4est_ghost_t * Ghost_layer, void *user_data,
               p4est_iter_volume_t iter_volume, p4est_iter_face_t iter_face,
#ifdef P4_TO_P8
               p8est_iter_edge_t iter_edge,
#endif
               p4est_iter_corner_t iter_corner)
{
  p4est_iterate_ext (p4est, Ghost_layer, user_data, iter_volume, iter_face,
#ifdef P4_TO_P8
                     iter_edge,
#endif
                     iter_corner, 0);
}
