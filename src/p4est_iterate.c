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

/* tier ring functions */
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
    ring = sc_array_index_int (tier_rings, i);
    ring->next = 0;
    sc_array_init (&(ring->tiers), sizeof (p4est_iter_tier_t));
    sc_array_resize (&(ring->tiers), (size_t) tier_ring_max);
    for (j = 0; j < tier_ring_max; j++) {
      tier = sc_array_index_int (&(ring->tiers), j);
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
    ring = sc_array_index (tier_rings, zz);
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

/* loop arg functions */
typedef struct p4est_iter_loop_args
{
  int                 alloc_size;
#ifdef P4_TO_P8
  bool                loop_edge;
#endif
  bool                loop_corner;

  int                 level;
  int                *level_num;        /* an array that keeps track of which
                                           branch we take at each step in the
                                           heirarchical search areas */
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
  int                *test_level;       /* a testing variable used in
                                           edge_iterate: passed as an argument
                                           to avoid using alloc/free on each 
                                           call */
  bool               *refine;   /* a testing variable used in edge_iterate:
                                   passed as an argument to avoid using
                                   alloc/free on each call */
  sc_array_t         *tier_rings;
}
p4est_iter_loop_args_t;

static p4est_iter_loop_args_t *
p4est_iter_loop_args_new (p4est_connectivity_t * conn,
#ifdef P4_TO_P8
                          p8est_iter_edge_t iter_edge, int *Max_edge_size,
#endif
                          p4est_iter_corner_t iter_corner,
                          int *Max_corner_size, sc_array_t * ghost_layer,
                          int num_procs)
{
  int                 i;
  int                 c;
  int                 alloc_size;
#ifndef P4_TO_P8
  const p4est_topidx_t *ctt_offset = conn->vtt_offset;
  p4est_topidx_t      num_corners = conn->num_vertices;
#else
  p4est_topidx_t      num_corners = conn->num_corners;
  const p4est_topidx_t *ctt_offset = conn->ctt_offset;
  int                 e;
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
  loop_args->start_idx2 = P4EST_ALLOC (int, alloc_size / 2);
  loop_args->quad_idx2 = P4EST_ALLOC (int, alloc_size / 2);
  loop_args->quadrants = P4EST_ALLOC (sc_array_t *, alloc_size);
  loop_args->index = P4EST_ALLOC (size_t *, alloc_size);
  for (i = 0; i < alloc_size; i++) {
    loop_args->index[i] =
      P4EST_ALLOC (size_t, (P4EST_QMAXLEVEL + 1) * P4EST_ITER_STRIDE);
    if (i % 2) {
      loop_args->quadrants[i] = ghost_layer;
    }
  }
  loop_args->first_index = P4EST_ALLOC (size_t, alloc_size);
  loop_args->count = P4EST_ALLOC (size_t, alloc_size);
  loop_args->test = P4EST_ALLOC (p4est_quadrant_t *, alloc_size);
  loop_args->test_level = P4EST_ALLOC (int, alloc_size);
  loop_args->refine = P4EST_ALLOC (bool, alloc_size / 2);

  loop_args->tier_rings = p4est_iter_tier_rings_new (num_procs);

#ifdef P4_TO_P8
  loop_args->loop_edge = ((iter_corner != NULL) || (iter_edge != NULL));
#endif
  loop_args->loop_corner = (iter_corner != NULL);

#ifdef P4_TO_P8
  *Max_edge_size = max_edge_size;
#endif
  *Max_corner_size = max_corner_size;

  return loop_args;
}

static void
p4est_iter_loop_args_destroy (p4est_iter_loop_args_t * loop_args)
{
  int                 i;
  int                 alloc_size = loop_args->alloc_size;

  P4EST_FREE (loop_args->level_num);
  P4EST_FREE (loop_args->start_idx2);
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

static void
p4est_iter_init_loop_volume (p4est_iter_loop_args_t * loop_args,
                             p4est_topidx_t t, p4est_t * p4est,
                             sc_array_t * ghost_quads, size_t * ghost_offsets)
{
  const int           left = 0;
  const int           right = 1;
  const int           local = 0;
  const int           ghost = 1;
  int                 i;
  sc_array_t         *trees = p4est->trees;
  p4est_tree_t       *tree;
  sc_array_t         *local_quads;
  size_t              first_ghost_quad = ghost_offsets[t];
  size_t              stop_ghost_quad = ghost_offsets[t + 1];

  tree = p4est_array_index_topidx (trees, t);
  local_quads = &(tree->quadrants);

  loop_args->level = 0;
  loop_args->level_num[0] = 0;

  loop_args->start_idx2[0] = 0;
  for (i = left; i <= right; i++) {
    loop_args->index[i * 2 + local][0] = 0;
    loop_args->index[i * 2 + local][1] = local_quads->elem_count;
    loop_args->index[i * 2 + ghost][0] = first_ghost_quad;
    loop_args->index[i * 2 + ghost][1] = stop_ghost_quad;
  }
  for (i = 0; i < 4; i++) {
    loop_args->quadrants[i] = (i % 2) ? ghost_quads : local_quads;
  }
#ifdef P4_TO_P8
  if (loop_args->loop_edge) {
    for (; i < 8; i++) {
      loop_args->quadrants[i] = (i % 2) ? ghost_quads : local_quads;
    }
  }
#endif
  if (loop_args->loop_corner) {
    for (; i < 2 * P4EST_CHILDREN; i++) {
      loop_args->quadrants[i] = (i % 2) ? ghost_quads : local_quads;
    }
  }
}

static void
p4est_iter_init_loop_face (p4est_iter_loop_args_t * loop_args,
                           p4est_topidx_t t, p4est_topidx_t nt,
                           p4est_t * p4est, sc_array_t * ghost_quads,
                           size_t * ghost_offsets)
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
  size_t              left_first_ghost = ghost_offsets[t];
  size_t              left_stop_ghost = ghost_offsets[t + 1];
  size_t              right_first_ghost = ghost_offsets[nt];
  size_t              right_stop_ghost = ghost_offsets[nt + 1];

  tree = p4est_array_index_topidx (trees, t);
  left_local_quads = &(tree->quadrants);
  tree = p4est_array_index_topidx (trees, nt);
  right_local_quads = &(tree->quadrants);

  loop_args->level = 0;
  loop_args->level_num[0] = 0;

  loop_args->start_idx2[left] = 0;
  loop_args->start_idx2[right] = 0;

  loop_args->index[left * 2 + local][0] = 0;
  loop_args->index[left * 2 + local][1] = left_local_quads->elem_count;
  loop_args->index[left * 2 + ghost][0] = left_first_ghost;
  loop_args->index[left * 2 + ghost][1] = left_stop_ghost;

  loop_args->index[right * 2 + local][0] = 0;
  loop_args->index[right * 2 + local][1] = right_local_quads->elem_count;
  loop_args->index[right * 2 + ghost][0] = right_first_ghost;
  loop_args->index[right * 2 + ghost][1] = right_stop_ghost;

  loop_args->quadrants[left * 2 + local] = left_local_quads;
  loop_args->quadrants[left * 2 + ghost] = ghost_quads;
  loop_args->quadrants[right * 2 + local] = right_local_quads;
  loop_args->quadrants[right * 2 + ghost] = ghost_quads;
  i = 4;
#ifdef P4_TO_P8
  if (loop_args->loop_edge) {
    for (; i < 8; i++) {
      loop_args->quadrants[i] = (i % 2) ? ghost_quads :
        ((i / 2) % 2) ? right_local_quads : left_local_quads;
    }
  }
#endif
  if (loop_args->loop_corner) {
    for (; i < 2 * P4EST_CHILDREN; i++) {
      loop_args->quadrants[i] = (i % 2) ? ghost_quads :
        ((i / 2) % 2) ? right_local_quads : left_local_quads;
    }
  }
}

static void
p4est_iter_init_loop_outside_face (p4est_iter_loop_args_t * loop_args,
                                   p4est_topidx_t t, p4est_t * p4est,
                                   sc_array_t * ghost_quads,
                                   size_t * ghost_offsets)
{
  const int           local = 0;
  const int           ghost = 1;
  int                 i;
  sc_array_t         *trees = p4est->trees;
  p4est_tree_t       *tree;
  sc_array_t         *local_quads;
  size_t              first_ghost_quad = ghost_offsets[t];
  size_t              stop_ghost_quad = ghost_offsets[t + 1];

  tree = p4est_array_index_topidx (trees, t);
  local_quads = &(tree->quadrants);

  loop_args->level = 0;
  loop_args->level_num[0] = 0;

  loop_args->start_idx2[0] = 0;

  loop_args->index[local][0] = 0;
  loop_args->index[local][1] = local_quads->elem_count;
  loop_args->index[ghost][0] = first_ghost_quad;
  loop_args->index[ghost][1] = stop_ghost_quad;

  loop_args->quadrants[local] = local_quads;
  loop_args->quadrants[ghost] = ghost_quads;
  i = 2;
#ifdef P4_TO_P8
  if (loop_args->loop_edge) {
    for (; i < 4; i++) {
      loop_args->quadrants[i] = (i % 2) ? ghost_quads : local_quads;
    }
  }
#endif
  if (loop_args->loop_corner) {
    for (; i < P4EST_CHILDREN; i++) {
      loop_args->quadrants[i] = (i % 2) ? ghost_quads : local_quads;
    }
  }
}

#ifdef P4_TO_P8
static void
p8est_iter_init_loop_edge (p4est_iter_loop_args_t * loop_args,
                           p8est_t * p8est, sc_array_t * ghost_quads,
                           size_t * ghost_offsets,
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
    side = sc_array_index (&(info->sides), zz);
    t = side->treeid;
    tree = p4est_array_index_topidx (trees, t);
    local_quads = &(tree->quadrants);
    loop_args->start_idx2[zz] = 0;
    loop_args->index[zz * 2 + local][0] = 0;
    loop_args->index[zz * 2 + local][1] = local_quads->elem_count;
    loop_args->index[zz * 2 + ghost][0] = ghost_offsets[t];
    loop_args->index[zz * 2 + ghost][1] = ghost_offsets[t + 1];
    loop_args->quadrants[zz * 2 + local] = local_quads;
    loop_args->quadrants[zz * 2 + ghost] = ghost_quads;
    if (loop_args->loop_corner) {
      loop_args->quadrants[(limit + zz) * 2 + local] = local_quads;
      loop_args->quadrants[(limit + zz) * 2 + ghost] = ghost_quads;
    }
  }
}
#endif

static void
p4est_iter_init_loop_corner (p4est_iter_loop_args_t * loop_args,
                             p4est_t * p4est, sc_array_t * ghost_quads,
                             size_t * ghost_offsets,
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
    side = sc_array_index (&(info->sides), zz);
    t = side->treeid;
    tree = p4est_array_index_topidx (trees, t);
    local_quads = &(tree->quadrants);
    loop_args->start_idx2[zz] = 0;
    loop_args->index[zz * 2 + local][0] = 0;
    loop_args->index[zz * 2 + local][1] = local_quads->elem_count;
    loop_args->index[zz * 2 + ghost][0] = ghost_offsets[t];
    loop_args->index[zz * 2 + ghost][1] = ghost_offsets[t + 1];
    loop_args->quadrants[zz * 2 + local] = local_quads;
    loop_args->quadrants[zz * 2 + ghost] = ghost_quads;
  }
}

typedef struct p4est_iter_corner_args
{
  int                 num_sides;
  p4est_iter_loop_args_t *loop_args;
  p4est_iter_corner_info_t *info;
}
p4est_iter_corner_args_t;

static              bool
p4est_iter_init_corner (p4est_iter_corner_args_t * args, p4est_t * p4est,
                        p4est_topidx_t t, int c)
{
  int                 i;
  int                 f;
  int                 count = 0;
  sc_array_t          quads;
  sc_array_t          treeids;
  p4est_quadrant_t    orig, temp, *ptemp;
  p4est_topidx_t      nt;
  p4est_connectivity_t *conn = p4est->connectivity;
  p4est_iter_corner_info_t *info = args->info;
  p4est_iter_corner_side_t *cside;
  int                *start_idx2 = args->loop_args->start_idx2;
#ifdef P4_TO_P8
  int                 e;
#endif
  size_t              zz;

  temp.x = 0;
  temp.y = 0;
#ifdef P4_TO_P8
  temp.z = 0;
#endif
  temp.level = 1;

  p4est_quadrant_sibling (&temp, &orig, c);

  sc_array_resize (&(info->sides), count + 1);
  cside = sc_array_index_int (&(info->sides), count);
  cside->corner = c;
  cside->treeid = t;
  start_idx2[count++] = 0;

  for (i = 0; i < P4EST_DIM; i++) {
    f = p4est_corner_faces[c][i];
    nt = p4est_quadrant_face_neighbor_extra (&orig, t, f, &temp, conn);
    if (nt == -1) {
      continue;
    }
    sc_array_resize (&(info->sides), count + 1);
    cside = sc_array_index_int (&(info->sides), count);
    cside->corner = p4est_quadrant_child_id (&temp);
    cside->treeid = nt;
    start_idx2[count++] = 0;
  }
#ifdef P4_TO_P8
  for (i = 0; i < 3; i++) {
    e = p8est_corner_edges[c][i];

    sc_array_init (&quads, sizeof (p4est_quadrant_t));
    sc_array_init (&treeids, sizeof (p4est_topidx_t));
    p8est_quadrant_edge_neighbor_extra (&orig, t, e, &quads, &treeids, conn);
    for (zz = 0; zz < quads.elem_count; zz++) {
      nt = *((p4est_topidx_t *) sc_array_index (&treeids, zz));
      ptemp = sc_array_index (&quads, zz);
      sc_array_resize (&(info->sides), count + 1);
      cside = sc_array_index_int (&(info->sides), count);
      cside->corner = p4est_quadrant_child_id (ptemp);
      cside->treeid = nt;
      start_idx2[count++] = 0;
    }
    sc_array_reset (&quads);
    sc_array_reset (&treeids);
  }
#endif

  sc_array_init (&quads, sizeof (p4est_quadrant_t));
  sc_array_init (&treeids, sizeof (p4est_topidx_t));
  p4est_quadrant_corner_neighbor_extra (&orig, t, c, &quads, &treeids, conn);
  for (zz = 0; zz < quads.elem_count; zz++) {
    nt = *((p4est_topidx_t *) sc_array_index (&treeids, zz));
    ptemp = sc_array_index (&quads, zz);
    sc_array_resize (&(info->sides), count + 1);
    cside = sc_array_index_int (&(info->sides), count);
    cside->corner = p4est_quadrant_child_id (ptemp);
    cside->treeid = nt;
    start_idx2[count++] = 0;
  }
  sc_array_reset (&quads);
  sc_array_reset (&treeids);

  args->num_sides = count;

  for (i = 1; i < count; i++) {
    cside = sc_array_index_int (&(info->sides), i);
    if (cside->treeid > t) {
      return false;
    }
    if (cside->treeid == t && cside->corner > c) {
      return false;
    }
  }

  return true;
}

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

  p4est_iter_loop_args_t *loop_args = args->loop_args;
  int                 level = loop_args->level;
  int                 num_sides = args->num_sides;
  const int          *start_idx2 = loop_args->start_idx2;
  int                *quad_idx2 = loop_args->quad_idx2;
  int                 this_corner;
  sc_array_t        **quadrants = loop_args->quadrants;
  size_t            **index = loop_args->index;
  size_t             *first_index = loop_args->first_index;
  size_t             *count = loop_args->count;
  p4est_quadrant_t  **test = loop_args->test;
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
  p4est_iter_loop_args_t *loop_args;
  int                 num_sides;
  sc_array_t         *common_corners[2];        /* for each side of the edge,
                                                   there are two corners that
                                                   touch the edge */
  p4est_iter_corner_args_t *corner_args;
  p8est_iter_edge_info_t *info;
}
p8est_iter_edge_args_t;

static void
p8est_iter_init_corner_from_edge (p4est_iter_corner_args_t * args,
                                  p8est_iter_edge_args_t * edge_args)
{
  int                 j, k;
  p8est_iter_corner_info_t *info = args->info;
  p8est_iter_edge_side_t *eside;
  p8est_iter_corner_side_t *cside;
  sc_array_t        **common_corners = edge_args->common_corners;
  int                *c_start_idx2 = args->loop_args->start_idx2;

  args->num_sides = edge_args->num_sides * 2;
  sc_array_resize (&(info->sides), (size_t) args->num_sides);
  for (j = 0; j < args->num_sides; j++) {
    k = j % edge_args->num_sides;
    eside = sc_array_index_int (&(edge_args->info->sides), k);
    cside = sc_array_index_int (&(info->sides), j);
    cside->treeid = eside->treeid;
    if (j == k) {
      cside->corner = *((int *) sc_array_index (common_corners[1], k));
      c_start_idx2[j] = *((int *) sc_array_index (common_corners[0], k));
    }
    else {
      cside->corner = *((int *) sc_array_index (common_corners[0], k));
      c_start_idx2[j] = *((int *) sc_array_index (common_corners[1], k));
    }
  }
}

static              bool
p8est_iter_init_edge (p8est_iter_edge_args_t * args, p8est_t * p8est,
                      p4est_topidx_t t, int e)
{
  int                 i;
  int                 f;
  int                 c0, c1, *cc;
  int                 count = 0;
  int                 orig_orient;
  sc_array_t          quads;
  sc_array_t          treeids;
  p8est_quadrant_t    orig, tempq, tempr, *ptemp;
  p4est_topidx_t      nt;
  p8est_connectivity_t *conn = p8est->connectivity;
  p8est_iter_edge_info_t *info = args->info;
  p8est_iter_edge_side_t *eside;
  int                *start_idx2 = args->loop_args->start_idx2;
  size_t              zz;
  p4est_topidx_t      it;
  p4est_topidx_t      edge = conn->tree_to_edge[t * 12 + e];
  int8_t              ete;
  sc_array_t        **common_corners = args->common_corners;

  orig_orient = 0;
  for (it = conn->ett_offset[edge]; it < conn->ett_offset[edge + 1]; it++) {
    nt = conn->edge_to_tree[it];
    if (nt == t) {
      ete = conn->edge_to_edge[it];
      if ((int) (ete % 12) == e) {
        orig_orient = ete / 12;
      }
    }
  }

  tempq.x = 0;
  tempq.y = 0;
  tempq.z = 0;
  tempq.level = 1;

  c0 = p8est_edge_corners[e][0];
  c1 = p8est_edge_corners[e][1];

  if (orig_orient == 0) {
    p4est_quadrant_sibling (&tempq, &tempr, c1);
    tempr.level++;
    p4est_quadrant_sibling (&tempr, &orig, c0);
  }
  else {
    p4est_quadrant_sibling (&tempq, &tempr, c1);
    tempr.level++;
    p4est_quadrant_sibling (&tempr, &orig, c0);
  }
  P4EST_ASSERT (p8est_quadrant_touches_edge (&orig, e, true));

  sc_array_resize (common_corners[0], count + 1);
  cc = sc_array_index_int (common_corners[0], count);
  *cc = (orig_orient == 0) ? c0 : c1;
  sc_array_resize (common_corners[1], count + 1);
  cc = sc_array_index_int (common_corners[1], count);
  *cc = (orig_orient == 0) ? c1 : c0;
  sc_array_resize (&(info->sides), count + 1);
  eside = sc_array_index_int (&(info->sides), count);
  eside->edge = e;
  eside->treeid = t;
  eside->orientation = orig_orient;
  start_idx2[count++] = 0;

  for (i = 0; i < 2; i++) {
    f = p8est_edge_faces[e][i];
    nt = p8est_quadrant_face_neighbor_extra (&orig, t, f, &tempq, conn);
    if (nt == -1) {
      continue;
    }
    c0 = p8est_quadrant_child_id (&tempq);
    p4est_quadrant_parent (&tempq, &tempr);
    c1 = p8est_quadrant_child_id (&tempr);
    sc_array_resize (common_corners[0], count + 1);
    cc = sc_array_index_int (common_corners[0], count);
    *cc = (orig_orient == 0) ? c0 : c1;
    sc_array_resize (common_corners[1], count + 1);
    cc = sc_array_index_int (common_corners[1], count);
    *cc = (orig_orient == 0) ? c1 : c0;
    sc_array_resize (&(info->sides), count + 1);
    eside = sc_array_index_int (&(info->sides), count);
    eside->orientation = (c0 < c1) ? 0 : 1;
    eside->edge = p8est_child_corner_edges[c0][c1];
    P4EST_ASSERT (eside->edge >= 0);
    eside->treeid = nt;
    start_idx2[count++] = 0;
  }

  sc_array_init (&quads, sizeof (p4est_quadrant_t));
  sc_array_init (&treeids, sizeof (p4est_topidx_t));
  p8est_quadrant_edge_neighbor_extra (&orig, t, e, &quads, &treeids, conn);
  for (zz = 0; zz < quads.elem_count; zz++) {
    nt = *((p4est_topidx_t *) sc_array_index (&treeids, zz));
    ptemp = sc_array_index (&quads, zz);
    c0 = p8est_quadrant_child_id (ptemp);
    p4est_quadrant_parent (ptemp, &tempr);
    c1 = p8est_quadrant_child_id (&tempr);
    sc_array_resize (common_corners[0], count + 1);
    cc = sc_array_index_int (common_corners[0], count);
    *cc = (orig_orient == 0) ? c0 : c1;
    sc_array_resize (common_corners[1], count + 1);
    cc = sc_array_index_int (common_corners[1], count);
    *cc = (orig_orient == 0) ? c1 : c0;
    sc_array_resize (&(info->sides), count + 1);
    eside = sc_array_index_int (&(info->sides), count);
    eside->orientation = (c0 < c1) ? 0 : 1;
    eside->edge = p8est_child_corner_edges[c0][c1];
    P4EST_ASSERT (eside->edge >= 0);
    eside->treeid = nt;
    start_idx2[count++] = 0;
  }
  sc_array_reset (&quads);
  sc_array_reset (&treeids);

  args->num_sides = count;

  for (i = 1; i < count; i++) {
    eside = sc_array_index_int (&(info->sides), i);
    if (eside->treeid > t) {
      return false;
    }
    if (eside->treeid == t && eside->edge > e) {
      return false;
    }
  }

  return true;
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
  int                *start_idx2 = loop_args->start_idx2;
  int                *level_num = loop_args->level_num;
  sc_array_t        **quadrants = loop_args->quadrants;
  size_t            **index = loop_args->index;
  size_t             *first_index = loop_args->first_index;
  sc_array_t        **common_corners = args->common_corners;
  p8est_quadrant_t  **test = loop_args->test;
  size_t             *count = loop_args->count;
  int                *test_level = loop_args->test_level;
  int                *quad_idx2 = loop_args->quad_idx2;
  bool               *refine = loop_args->refine;
  int                *temp_int, *temp_int2;
  int                 i;
  int                *Level = &(loop_args->level);
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
        if (test_level[sidetype] == *Level) {
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
            level_num[*Level]++;
            goto change_search_area;
          }
        }
      }
    }
    if (stop_refine) {
      for (side = 0; side < num_sides; side++) {
        if (count[side * 2 + local] == 0 && count[side * 2 + ghost] == 0) {
          eside = sc_array_index_int (&(info->sides), side);
          eside->is_hanging = false;
          eside->is.full.quad = NULL;
          eside->is.full.is_local = false;
          eside->is.full.quadid = -1;
          refine[side] = false;
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
      level_num[*Level]++;
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
          p4est_iter_tier_insert (&test_view, *Level, index[sidetype] +
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
                quads[child_corner] = sc_array_index (quadrants[sidetype],
                                                      first_index[sidetype]);
                P4EST_ASSERT ((int) quads[child_corner]->level == *Level + 1);
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
      level_num[*Level]++;
      goto change_search_area;
    }
    /* if every side needed to be refined, then we descend along this branch to
     * this next level and search there */
    level_num[++(*Level)] = 0;
    level_idx2 += P4EST_ITER_STRIDE;

  change_search_area:
    /* if we tried to advance the search area on start_level, we've completed
     * the search */
    if (level_num[start_level] > 0) {
      break;
    }
    /* if we have tried to advance the search area past two branches, that
     * means that we have completed all of the branches on this level */
    if (level_num[*Level] == 2) {
      /* if we have a corner callback, we need to run it on the corner between
       * the edge branches on this level */
      if (iter_corner != NULL) {
        /* find the correct corner ids for the corner search areas, and copy
         * the index information necessary to run corner_iterate over to the
         * additional index arrays */
        for (i = 0; i < 2 * num_sides; i++) {
          temp_int = sc_array_index_int (common_corners[i / num_sides],
                                         i % num_sides);
          temp_idx2 = level_idx2 + *temp_int;
          index[i * 2 + local][temp_idx2] =
            index[(i % num_sides) * 2 + local][temp_idx2];
          index[i * 2 + local][temp_idx2 + 1] =
            index[(i % num_sides) * 2 + local][temp_idx2 + 1];
          index[i * 2 + ghost][temp_idx2] =
            index[(i % num_sides) * 2 + ghost][temp_idx2];
          index[i * 2 + ghost][temp_idx2 + 1] =
            index[(i % num_sides) * 2 + ghost][temp_idx2 + 1];
        }
        p8est_iter_init_corner_from_edge (corner_args, args);
        p4est_corner_iterate (corner_args, user_data, iter_corner);
      }
      /* now that we're done on this level, go up a level and over a branch */
      level_num[--(*Level)]++;
      level_idx2 -= P4EST_ITER_STRIDE;
      goto change_search_area;
    }

    /* at this point, we need to initialize the bounds of the search areas
     * for this new branch */
    all_empty = true;
    for (side = 0; side < num_sides; side++) {
      temp_int = sc_array_index_int (common_corners[level_num[*Level]], side);
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
      level_num[*Level]++;
      goto change_search_area;
    }
  }
  P4EST_ASSERT (*Level == start_level);
}
#endif

typedef struct p4est_iter_face_args
{
  p4est_iter_loop_args_t *loop_args;
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
  bool                outside_face;     /* indicates if we are at a tree
                                           boundary without a neighbor across
                                           the face */
#ifdef P4_TO_P8
  p8est_iter_edge_args_t *edge_args;
#endif
  p4est_iter_corner_args_t *corner_args;
  p4est_iter_face_info_t *info;
}
p4est_iter_face_args_t;

static void
p4est_iter_init_corner_from_face (p4est_iter_corner_args_t * args,
                                  p4est_iter_face_args_t * face_args)
{
  const int           ntc_str = P4EST_CHILDREN / 2;
  int                 j, k;
  p4est_iter_corner_info_t *info = args->info;
  p4est_iter_face_side_t *fside;
  p4est_iter_corner_side_t *cside;
  int                *num_to_child = face_args->num_to_child;
  int                *c_start_idx2 = args->loop_args->start_idx2;
  int                 limit = face_args->outside_face ? 1 : 2;
  int                 count = 0;

  args->num_sides = limit * ntc_str;
  sc_array_resize (&(info->sides), (size_t) args->num_sides);
  for (j = 0; j < ntc_str; j++) {
    for (k = 0; k < limit; k++) {
      fside = sc_array_index_int (&(face_args->info->sides), k);
      cside = sc_array_index_int (&(info->sides), count);
      cside->treeid = fside->treeid;
      cside->corner = num_to_child[k * ntc_str + (ntc_str - 1 - j)];
      c_start_idx2[count++] = num_to_child[k * ntc_str + j];
    }
  }
}

#ifdef P4_TO_P8
static void
p8est_iter_init_edge_from_face (p8est_iter_edge_args_t * args,
                                p4est_iter_face_args_t * face_args,
                                int dir, int side)
{
  const int           ntc_str = P4EST_CHILDREN / 2;
  int                 j, k;
  int                 c0, c1, *cc;
  p8est_iter_edge_info_t *info = args->info;
  p8est_iter_face_side_t *fside;
  p8est_iter_edge_side_t *eside;
  int                *num_to_child = face_args->num_to_child;
  int                *e_start_idx2 = args->loop_args->start_idx2;
  int                 limit = face_args->outside_face ? 1 : 2;
  int                 count = 0;
  int                 pos[2][2];
  sc_array_t        **common_corners = args->common_corners;

  pos[0][0] = 0;
  pos[0][1] = dir ? 2 : 1;
  pos[1][0] = dir ? 1 : 2;
  pos[1][1] = 3;

  args->num_sides = limit * ntc_str / 2;
  sc_array_resize (&(info->sides), (size_t) args->num_sides);
  sc_array_resize (args->common_corners[0], (size_t) args->num_sides);
  sc_array_resize (args->common_corners[1], (size_t) args->num_sides);
  for (j = 0; j < 2; j++) {
    for (k = 0; k < limit; k++) {
      cc = sc_array_index_int (common_corners[0], count);
      c0 = *cc = num_to_child[k * ntc_str + pos[1 - j][0]];
      cc = sc_array_index_int (common_corners[1], count);
      c1 = *cc = num_to_child[k * ntc_str + pos[1 - j][1]];
      fside = sc_array_index_int (&(face_args->info->sides), k);
      eside = sc_array_index_int (&(info->sides), count);
      eside->orientation = (c0 < c1) ? 0 : 1;
      eside->treeid = fside->treeid;
      eside->edge = p8est_child_corner_edges[c0][c1];
      e_start_idx2[count++] = num_to_child[k * ntc_str + pos[j][side]];
    }
  }
}
#endif

static              bool
p4est_iter_init_face (p4est_iter_face_args_t * args, p4est_t * p4est,
                      p4est_topidx_t t, int f, int *Nt)
{
  const int           ntc_str = P4EST_CHILDREN / 2;
  int                 i;
  int                 c;
#ifndef P4_TO_P8
  int                 rf = p4est_zface_to_rface[f];
  int                 nrf;
#endif
  int                 nf;
  int                 count = 0;
  p4est_quadrant_t    orig, tempq, tempr;
  p4est_topidx_t      nt;
  p4est_connectivity_t *conn = p4est->connectivity;
  p4est_iter_face_info_t *info = args->info;
  p4est_iter_face_side_t *fside;
  int                *num_to_child = args->num_to_child;
  int                *start_idx2 = args->loop_args->start_idx2;

  tempq.x = 0;
  tempq.y = 0;
#ifdef P4_TO_P8
  tempq.z = 0;
#endif
  tempq.level = 1;

  for (i = 0; i < ntc_str; i++) {
#ifndef P4_TO_P8
    c = p4est_face_corners[rf][i];
#else
    c = p4est_face_corners[f][i];
#endif
    num_to_child[i] = c;
    p4est_quadrant_sibling (&tempq, &orig, c);
#ifndef P4_TO_P8
    nt = p4est_quadrant_face_neighbor_extra (&orig, t, rf, &tempr, conn);
#else
    nt = p4est_quadrant_face_neighbor_extra (&orig, t, f, &tempr, conn);
#endif
    if (nt != -1) {
      num_to_child[ntc_str + i] = p4est_quadrant_child_id (&tempr);
    }
  }
  args->outside_face = (nt == -1);

  sc_array_resize (&(info->sides), count + 1);
  fside = sc_array_index_int (&(info->sides), count);
  fside->face = f;
  fside->treeid = t;
  start_idx2[count++] = 0;
  info->orientation = 0;

  if (nt != -1) {
    sc_array_resize (&(info->sides), count + 1);
    fside = sc_array_index_int (&(info->sides), count);
    fside->treeid = nt;
#ifndef P4_TO_P8
    nrf = conn->tree_to_face[4 * t + rf];
    info->orientation = nrf / (2 * P4EST_DIM);
    nrf = nrf % (2 * P4EST_DIM);
    nf = p4est_rface_to_zface[nrf];
#else
    nf = conn->tree_to_face[2 * P4EST_DIM * t + f];
    info->orientation = nf / (2 * P4EST_DIM);
    nf = nf % (2 * P4EST_DIM);
#endif
    fside->face = nf;
    start_idx2[count++] = 0;
  }

  if (nt != -1) {
    if (nt > t) {
      return false;
    }
    if (nt == t && nf > f) {
      return false;
    }
  }

  *Nt = nt;

  return true;
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
  int                *start_idx2 = loop_args->start_idx2;
  int                *level_num = loop_args->level_num;
  sc_array_t        **quadrants = loop_args->quadrants;
  size_t            **index = loop_args->index;
  size_t             *first_index = loop_args->first_index;
  int                *num_to_child = args->num_to_child;
  p4est_quadrant_t  **test = loop_args->test;
  size_t             *count = loop_args->count;
  int                *test_level = loop_args->test_level;
  int                *quad_idx2 = loop_args->quad_idx2;
  bool               *refine = loop_args->refine;
  int                 limit;

  int                 i, j;
  int                *Level = &(loop_args->level);
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
  sc_array_t         *tier_rings = loop_args->tier_rings;
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
  level_idx2 = start_level * P4EST_ITER_STRIDE;

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
        if (test_level[sidetype] == *Level) {
          P4EST_ASSERT (count[sidetype] == 1);
          P4EST_ASSERT (count[side * 2 + (type ^ 1)] == 0);
          if (iter_face != NULL) {
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
                    test_level[nsidentype] == *Level) {
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
                  level_num[*Level]++;
                  goto change_search_area;
                }
              }
              if (count[n_side * 2 + local] == 0 &&
                  count[n_side * 2 + ghost] == 0) {
                fside = sc_array_index_int (&(info->sides), n_side);
                fside->is_hanging = false;
                fside->is.full.quad = NULL;
                fside->is.full.is_local = false;
                fside->is.full.quadid = -1;
                iter_face (info, user_data);
                level_num[*Level]++;
                goto change_search_area;
              }
            }
            else {
              iter_face (info, user_data);
              level_num[*Level]++;
              goto change_search_area;
            }
          }
          /* if there is no face callback, i.e. we are running face_iterate only
           * to find the edges/corners that live on faces, then we are done once
           * we find a side that does not need to be refined */
          else {
            level_num[*Level]++;
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
          p4est_iter_tier_insert (&test_view, *Level, index[sidetype] +
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
              quads[child_corner] = sc_array_index (quadrants[nsidentype],
                                                    first_index[nsidentype]);
              P4EST_ASSERT ((int) quads[child_corner]->level == *Level + 1);
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
        level_num[*Level]++;
        goto change_search_area;
      }
    }
    /* if we refined both sides, we descend to the next level from this branch
     * and continue searching there */
    level_num[++(*Level)] = 0;
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
    if (level_num[*Level] == P4EST_CHILDREN / 2) {
#ifdef P4_TO_P8
      /* if we have an edge callback, we need to run it on all of the edges
       * between the face branches on this level */
      if (iter_edge != NULL || iter_corner != NULL) {
        if (!args->outside_face) {
          for (dir = 0; dir < 2; dir++) {

            /* we have to set up the common corners on either side of the
             * desired edge */
            for (side = 0; side < 2; side++) {
              p8est_iter_init_edge_from_face (edge_args, args, dir, side);
              for (j = 0; j < 4; j++) {
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
              p8est_iter_init_edge_from_face (edge_args, args, dir, side);
              for (j = 0; j < 2; j++) {
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
        /* find the correct corner ids for the corner search areas, and copy
         * the index information necessary to run corner_iterate over to the
         * additional index arrays */
        if (!args->outside_face) {
          for (j = 0; j < P4EST_CHILDREN; j++) {
            temp_idx2 =
              level_idx2 + num_to_child[(j % 2) * ntc_str + (j / 2)];
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
            temp_idx2 = level_idx2 + num_to_child[j];
            index[j * 2 + local][temp_idx2] = index[local][temp_idx2];
            index[j * 2 + local][temp_idx2 + 1] = index[local]
              [temp_idx2 + 1];
            index[j * 2 + ghost][temp_idx2] = index[ghost][temp_idx2];
            index[j * 2 + ghost][temp_idx2 + 1] = index[ghost]
              [temp_idx2 + 1];
          }
        }
        p4est_iter_init_corner_from_face (corner_args, args);
        p4est_corner_iterate (corner_args, user_data, iter_corner);
      }

      /* now that we're done on this level, go up a level and over a branch */
      level_num[--(*Level)]++;
      level_idx2 -= P4EST_ITER_STRIDE;
      goto change_search_area;
    }

    /* at this point, we need to initialize the bounds of the search areas
     * for this new branch */
    for (side = left; side <= limit; side++) {
      quad_idx2[side] =
        level_idx2 + num_to_child[side * ntc_str + level_num[*Level]];
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
        level_num[*Level]++;
        goto change_search_area;
      }
    }
    else {
      if (!count[left * 2 + local]) {
        level_num[*Level]++;
        goto change_search_area;
      }
    }
  }
  P4EST_ASSERT (*Level == start_level);
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
p4est_iterate (p4est_t * p4est, sc_array_t * Ghost_layer, void *user_data,
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
  int                 i, j, side, type;
  int                 left_side, right_side;
  p4est_topidx_t      t;
#ifdef P4_TO_P8
  int                 k;
  p4est_topidx_t      jt;
#endif

  /* counts */

  sc_array_t          empty_ghost_layer;
  sc_array_t         *ghost_layer;

  /* trees and quadrant arrays */
  sc_array_t         *trees = p4est->trees;
  size_t             *tree_first_ghost;

  /* topology */
  p4est_connectivity_t *conn = p4est->connectivity;
  const int8_t       *ttf = conn->tree_to_face;
  const p4est_topidx_t *ttt = conn->tree_to_tree;
  p4est_topidx_t      nt;
  size_t              global_num_trees;
#ifdef P4_TO_P8
  const p4est_topidx_t *tte = conn->tree_to_edge;
  p4est_topidx_t      num_edges = conn->num_edges;
  const p4est_topidx_t *ett_offset = conn->ett_offset;
  const p4est_topidx_t *edge_to_tree = conn->edge_to_tree;
  p4est_topidx_t      this_edge;
  const int8_t       *edge_to_edge = conn->edge_to_edge;
#endif

  /* allocation size variables */
  int                 max_corner_size;
#ifdef P4_TO_P8
  int                 max_edge_size;
#endif

  /* search loop variables */
  p4est_iter_loop_args_t *loop_args;

  int                *Level;
  int                *level_num;
  int                *start_idx2;
  int                *quad_idx2;
  sc_array_t        **quadrants;
  size_t            **index;
  size_t             *first_index;
  size_t             *count;
  p4est_quadrant_t  **test;
  int                *test_level;
  bool               *refine;
  sc_array_t         *tier_rings;

  size_t             *this_index;
  int                 level_idx2;
  sc_array_t          test_view;

  /* face variables */
  int                 face[2];
#ifndef P4_TO_P8
  int                 rface[2];
#endif
  int                 num_to_child[P4EST_CHILDREN];
  p4est_iter_face_args_t face_args;
  p4est_iter_face_info_t face_info;
  p4est_iter_face_side_t face_sides[2], *fside;

  /* edge variables */
#ifdef P4_TO_P8
  int                 orientation;
  p8est_iter_edge_args_t edge_args;
  sc_array_t          common_corners[2];
  int                *e_common_corners[2], *e_edges = NULL;
  int                 perm, corner0, corner1, ne, dir;
  p4est_topidx_t     *e_treeids = NULL;
  p8est_iter_edge_info_t edge_info;
  p8est_iter_edge_side_t *edge_sides = NULL, *eside;
#endif

  /* corner variables */
  int                *c_corners = NULL;
  p4est_topidx_t     *c_treeids = NULL;
  p4est_iter_corner_args_t corner_args;
  p4est_iter_corner_info_t corner_info;
  p4est_iter_corner_side_t *corner_sides = NULL, *cside;

  /* volume info */
  p4est_iter_volume_info_t info;

  /* multiuse variable */
  int                 modulus;

  P4EST_ASSERT (p4est_is_valid (p4est));
  if (Ghost_layer == NULL) {
    sc_array_init (&empty_ghost_layer, sizeof (p4est_quadrant_t));
    ghost_layer = &empty_ghost_layer;
  }
  else {
    ghost_layer = Ghost_layer;
  }
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
    if (iter_volume == NULL) {
      return;
    }
    p4est_volume_iterate (p4est, ghost_layer, user_data, iter_volume);
    return;
  }

  /** initialize arrays that keep track of where we are in the search */
  loop_args = p4est_iter_loop_args_new (conn,
#ifdef P4_TO_P8
                                        iter_edge, &max_edge_size,
#endif
                                        iter_corner, &max_corner_size,
                                        ghost_layer, p4est->mpisize);
  level_num = loop_args->level_num;
  start_idx2 = loop_args->start_idx2;
  quad_idx2 = loop_args->quad_idx2;
  quadrants = loop_args->quadrants;
  index = loop_args->index;
  first_index = loop_args->first_index;
  count = loop_args->count;
  test = loop_args->test;
  test_level = loop_args->test_level;
  refine = loop_args->refine;
  tier_rings = loop_args->tier_rings;
  Level = &(loop_args->level);

  /** Divide the ghost_layer by p.which_tree */
  global_num_trees = trees->elem_count;
  tree_first_ghost = P4EST_ALLOC (size_t, global_num_trees + 1);
  p4est_split_ghost_layer_by_tree (ghost_layer, tree_first_ghost,
                                   global_num_trees);

  /** set up the arguments passed to face_iterate that are invariant */
  face_args.loop_args = loop_args;
  face_args.num_to_child = num_to_child;
#ifdef P4_TO_P8
  face_args.edge_args = &edge_args;
#endif
  face_args.corner_args = &corner_args;
  face_info.p4est = p4est;
  face_info.ghost_layer = ghost_layer;
  sc_array_init_data (&face_info.sides, face_sides, sizeof (*face_sides), 2);
  face_args.info = &face_info;

  /** set up the arguments passed to corner_iterate that are invariant */
  if (iter_corner != NULL) {
    corner_args.loop_args = loop_args;
    c_corners = P4EST_ALLOC (int, max_corner_size);
    c_treeids = P4EST_ALLOC (p4est_topidx_t, max_corner_size);
    corner_sides = P4EST_ALLOC (p4est_iter_corner_side_t, max_corner_size);
    sc_array_init_data (&corner_info.sides, corner_sides,
                        sizeof (*corner_sides), (size_t) max_corner_size);
    corner_info.p4est = p4est;
    corner_info.ghost_layer = ghost_layer;
    corner_args.info = &corner_info;
  }

#ifdef P4_TO_P8
  /** set up the arguments passed to edge_iterate that are invariant */
  if (iter_edge != NULL || iter_corner != NULL) {
    edge_args.loop_args = loop_args;
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
    edge_args.common_corners[0] = &(common_corners[0]);
    edge_args.common_corners[1] = &(common_corners[1]);
    edge_args.corner_args = &corner_args;
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
    p4est_iter_init_loop_volume (loop_args, t, p4est, ghost_layer,
                                 tree_first_ghost);

    /* since we are intra-tree, both treeids are the current tree,
     * orientation is always 0 */
    face_args.outside_face = false;
    face_info.orientation = 0;
    sc_array_resize (&face_info.sides, 2);
    for (side = left; side <= right; side++) {
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
        cside = sc_array_index_int (&corner_info.sides, i);
        cside->treeid = t;
      }
    }

    info.treeid = t;

    level_idx2 = (*Level) * P4EST_ITER_STRIDE;

    /* if we try to advance the branch number at level 0, then we have
     * completed the search through the tree */
    while (level_num[0] == 0) {
      /* if we have advanced the branch number to P4EST_CHILDREN, we have
       * completed the search at this level, and we now set up the face_args
       * for the faces between the quad branches at this level */
      if (level_num[*Level] == P4EST_CHILDREN) {
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
        level_num[--(*Level)]++;
        level_idx2 -= P4EST_ITER_STRIDE;
        continue;
      }

      /* quad_idx[0] now gives the location in index[type] of the bounds
       * of the current search area, from which we get the first quad
       * and the count */
      quad_idx2[0] = level_idx2 + level_num[*Level];
      for (type = local; type <= ghost; type++) {
        this_index = index[type];
        first_index[type] = this_index[quad_idx2[0]];
        count[type] = (this_index[quad_idx2[0] + 1] - first_index[type]);
      }
      /* if there are no local quadrants, we are done with this search area,
       * and we advance to the next branch at this level */
      if (!count[local]) {
        level_num[*Level]++;
        continue;
      }

      for (type = local; type <= ghost; type++) {
        if (count[type]) {
          test[type] = sc_array_index (quadrants[type], first_index[type]);
          test_level[type] = (int) test[type]->level;
          /* if the test quad form this type is the size of the search area,
           * we run the quad callback */
          if (test_level[type] == *Level) {
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
        level_num[*Level]++;
        continue;
      }

      /* otherwise, we need to refine the search areas, and place the
       * bounds of the refined search areas on the next tier in index[type]*/
      quad_idx2[0] = level_idx2 + P4EST_ITER_STRIDE;
      for (type = local; type <= ghost; type++) {
        sc_array_init_view (&test_view, quadrants[type],
                            first_index[type], count[type]);
        this_index = index[type];
        p4est_iter_tier_insert (&test_view, *Level,
                                this_index + quad_idx2[0],
                                first_index[type], tier_rings, test[type]);
      }

      /* We descend to the next level and continue the search there */
      level_num[++(*Level)] = 0;
      level_idx2 += P4EST_ITER_STRIDE;
    }
    P4EST_ASSERT (*Level == 0);

    /* Now we need to run face_iterate on the faces between trees */
    for (i = 0; i < 2 * P4EST_DIM; i++) {
      if (p4est_iter_init_face (&face_args, p4est, t, i, &nt)) {
#ifdef P4_TO_P8
        /* we set the arrays in edge_args to the appropriate size and
         * make sure that the quadrants arrays are pointed at the appropriate
         * arrays, and that the tree ids reflect which side of the
         * face each edge-side is on*/
        if (iter_edge != NULL || iter_corner != NULL) {
          p8est_iter_init_edge_from_face (&edge_args, &face_args, 0, 0);
        }
#endif

        /* we set the arrays in edge_args the appropriate size and
         * make sure that the quadrants arrays are pointed at the appropriate
         * arrays, and that the tree ids reflect which side of the
         * face each corner-side is on*/
        if (iter_corner != NULL) {
          p4est_iter_init_corner_from_face (&corner_args, &face_args);
        }

        if (!face_args.outside_face) {
          p4est_iter_init_loop_face (loop_args, t, nt, p4est, ghost_layer,
                                     tree_first_ghost);
        }
        else {
          p4est_iter_init_loop_outside_face (loop_args, t, p4est, ghost_layer,
                                             tree_first_ghost);
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
        if (p8est_iter_init_edge (&edge_args, p4est, t, i)) {
          /* if there is a corner callback, initialize the quadrants arrays to
           * point to the correct arrays and make sure the arrays
           * are the correct sizes */
          if (iter_corner != NULL) {
            p8est_iter_init_corner_from_edge (&corner_args, &edge_args);
          }
          p8est_iter_init_loop_edge (loop_args, p4est, ghost_layer,
                                     tree_first_ghost, &edge_info);

          p8est_edge_iterate (&edge_args, user_data, iter_edge, iter_corner);
        }
      }
    }
#endif

    if (iter_corner != NULL) {
      for (i = 0; i < P4EST_CHILDREN; i++) {
        if (p4est_iter_init_corner (&corner_args, p4est, t, i)) {
          p4est_iter_init_loop_corner (loop_args, p4est, ghost_layer,
                                       tree_first_ghost, &corner_info);

          p4est_corner_iterate (&corner_args, user_data, iter_corner);
        }
      }
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

  p4est_iter_loop_args_destroy (loop_args);
}
