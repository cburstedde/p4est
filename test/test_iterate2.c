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
#include <p4est_algorithms.h>
#include <p4est_bits.h>
#include <p4est_extended.h>
#include <p4est_ghost.h>
#include <p4est_iterate.h>
#else
#include <p8est_algorithms.h>
#include <p8est_bits.h>
#include <p8est_extended.h>
#include <p8est_ghost.h>
#include <p8est_iterate.h>
#endif

#ifndef P4_TO_P8
static int          refine_level = 5;
#else
static int          refine_level = 3;
#endif

static int
refine_fn (p4est_t * p4est, p4est_topidx_t which_tree,
           p4est_quadrant_t * quadrant)
{
  int                 cid;

  if (which_tree == 2 || which_tree == 3) {
    return 0;
  }

  cid = p4est_quadrant_child_id (quadrant);

  if (cid == P4EST_CHILDREN - 1 ||
      (quadrant->x >= P4EST_LAST_OFFSET (P4EST_MAXLEVEL - 2) &&
       quadrant->y >= P4EST_LAST_OFFSET (P4EST_MAXLEVEL - 2)
#ifdef P4_TO_P8
       && quadrant->z >= P4EST_LAST_OFFSET (P4EST_MAXLEVEL - 2)
#endif
      )) {
    return 1;
  }
  if ((int) quadrant->level >= (refine_level - (int) (which_tree % 3))) {
    return 0;
  }
  if (quadrant->level == 1 && cid == 2) {
    return 1;
  }
  if (quadrant->x == P4EST_QUADRANT_LEN (2) &&
      quadrant->y == P4EST_LAST_OFFSET (2)) {
    return 1;
  }
  if (quadrant->y >= P4EST_QUADRANT_LEN (2)) {
    return 0;
  }

  return 1;
}

#ifndef P4_TO_P8
static int          face_offset = 1;
static int          corner_offset = 5;
static int          checks_per_quad = 9;
static int          check_to_type[9] = { 2, 1, 1, 1, 1, 0, 0, 0, 0 };
#else
static int          face_offset = 1;
static int          edge_offset = 7;
static int          corner_offset = 19;
static int          checks_per_quad = 27;
static int          check_to_type[27] = { 3,
  2, 2, 2, 2, 2, 2,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  0, 0, 0, 0, 0, 0, 0, 0
};
#endif

typedef struct iter_data
{
  int8_t              count_volume;
  int8_t              count_face;
  int8_t              ghost_face;
#ifdef P4_TO_P8
  int8_t              count_edge;
  int8_t              ghost_edge;
#endif
  int8_t              count_corner;
  int8_t              ghost_corner;
  int                *checks;
}
iter_data_t;

#if 0
/*@unused@*/
static void
volume_do_nothing (p4est_iter_volume_info_t * info, void *data)
{
};

/*@unused@*/
static void
face_do_nothing (p4est_iter_face_info_t * info, void *data)
{
};

#ifdef P4_TO_P8
/*@unused@*/
static void
edge_do_nothing (p8est_iter_edge_info_t * info, void *data)
{
};
#endif

/*@unused@*/
static void
corner_do_nothing (p4est_iter_corner_info_t * info, void *data)
{
};
#endif /* 0 */

static              int8_t
test_corner_side (p4est_t * p4est, p4est_iter_corner_side_t * side,
                  void *data)
{
  int                 corner = (int) side->corner;
  p4est_quadrant_t   *q = side->quad;
  iter_data_t        *iter_data = (iter_data_t *) data;
  int                *checks = iter_data->checks;
  p4est_locidx_t      qid = side->quadid;
  p4est_topidx_t      t = side->treeid;
  p4est_tree_t       *tree;
  p4est_locidx_t      offset;

  tree = p4est_tree_array_index (p4est->trees, t);
  offset = tree->quadrants_offset;

  q = side->quad;
  if (q == NULL) {
    SC_CHECK_ABORT (!(iter_data->ghost_corner), "Iterate: empty corner side");
    return 1;
  }
  else {
    if (!side->is_ghost && iter_data->count_corner) {
      qid += offset;
      checks[qid * checks_per_quad + corner + corner_offset]++;
    }
  }

  return side->is_ghost;
}

static              int8_t
quad_is_in_corner_sides (p4est_quadrant_t * q, p4est_topidx_t t,
                         sc_array_t * sides)
{
  size_t              zz;
  size_t              zcount = sides->elem_count;
  p4est_iter_corner_side_t *cside;
  p4est_quadrant_t   *r;

  for (zz = 0; zz < zcount; zz++) {
    cside = p4est_iter_cside_array_index (sides, zz);
    if (cside->treeid == t) {
      r = cside->quad;
      if (r == NULL) {
        continue;
      }
      if (p4est_quadrant_is_equal (q, r)) {
        return 1;
      }
      if (p4est_quadrant_is_ancestor (q, r)) {
        return 1;
      }
    }
  }

  return 0;
}

static void
test_corner_boundary (p4est_iter_corner_info_t * info)
{
  int                 c, f, dir, which;
  int                 k, count;
  p4est_qcoord_t      end, xyz[P4EST_DIM];
  p4est_iter_corner_side_t *cside;
  p4est_quadrant_t   *q;

  SC_CHECK_ABORT (info->sides.elem_count > 0, "Empty corner iteration");
  cside = (p4est_iter_corner_side_t *) sc_array_index (&info->sides, 0);
  c = cside->corner;

  if (!info->tree_boundary) {
    SC_CHECK_ABORT (c == P4EST_CHILDREN - 1, "Not the lowest corner");
    return;
  }

  /* grab information about this quadrant */
  q = cside->quad;
  if (q == NULL) {
    SC_CHECK_ABORT (cside->is_ghost && cside->quadid == -1,
                    "Not a corner ghost");
    return;
  }
  end = P4EST_LAST_OFFSET (q->level);
  xyz[0] = q->x;
  xyz[1] = q->y;
#ifdef P4_TO_P8
  xyz[2] = q->z;
#endif

  /* check how many tree faces it is touching */
  count = 0;
  for (k = 0; k < P4EST_DIM; ++k) {
    f = p4est_corner_faces[c][k];
    dir = f >> 1;
    which = f & 1;
    count += (which == 0 && xyz[dir] == 0) || (which == 1 && xyz[dir] == end);
  }
  if (info->tree_boundary == P4EST_CONNECT_CORNER) {
    /* we are at a true inter-tree corner */
    SC_CHECK_ABORT (count == P4EST_DIM, "Not a tree boundary corner");
  }
#ifdef P4_TO_P8
  else if (info->tree_boundary == P8EST_CONNECT_EDGE) {
    /* we are a corner inside an inter-tree edge */
    SC_CHECK_ABORT (count == 2, "Not a tree edge boundary corner");
  }
#endif
  else if (info->tree_boundary == P4EST_CONNECT_FACE) {
    /* we are a corner inside an inter-tree face */
    SC_CHECK_ABORT (count == 1, "Not a tree face boundary corner");
  }
  else {
    SC_ABORT_NOT_REACHED ();
  }
}

static void
test_corner_adjacency (p4est_iter_corner_info_t * info, void *data)
{
  int                 i, j;
  int8_t              has_local = 0;
  int                 limit = (int) info->sides.elem_count;
  p4est_iter_corner_side_t *cside;
  p4est_quadrant_t    tempq, tempr;
  int                 c, f;
  p4est_topidx_t      t, nt;
  sc_array_t          quads, treeids;
  p4est_connectivity_t *conn = info->p4est->connectivity;
  size_t              zz;
  p4est_quadrant_t   *ptemp;
#ifdef P4_TO_P8
  int                 e;
#endif
  int8_t              min_level = P4EST_QMAXLEVEL;
  iter_data_t        *iter_data = (iter_data_t *) data;

  test_corner_boundary (info);

  for (i = 0; i < limit; i++) {
    cside = p4est_iter_cside_array_index_int (&info->sides, i);
    has_local = (!test_corner_side (info->p4est, cside, data) || has_local);
    if (cside->quad != NULL) {
      min_level = (cside->quad->level < min_level) ? cside->quad->level :
        min_level;
    }
  }
  SC_CHECK_ABORT (has_local, "Iterate: non local corner");

  sc_array_init (&quads, sizeof (p4est_quadrant_t));
  sc_array_init (&treeids, sizeof (p4est_topidx_t));
  for (i = 0; i < limit; i++) {
    cside = p4est_iter_cside_array_index_int (&info->sides, i);
    if (cside->quad == NULL) {
      continue;
    }
    c = (int) cside->corner;
    t = cside->treeid;
    tempq = *(cside->quad);
    tempq.x &= ((p4est_qcoord_t) - 1) << (P4EST_MAXLEVEL - min_level);
    tempq.y &= ((p4est_qcoord_t) - 1) << (P4EST_MAXLEVEL - min_level);
#ifdef P4_TO_P8
    tempq.z &= ((p4est_qcoord_t) - 1) << (P4EST_MAXLEVEL - min_level);
#endif
    tempq.level = min_level;
    P4EST_ASSERT (p4est_quadrant_is_valid (&tempq));
    for (j = 0; j < P4EST_DIM; j++) {
      f = p4est_corner_faces[c][j];
      nt = p4est_quadrant_face_neighbor_extra (&tempq, t, f, &tempr, NULL,
                                               conn);
      if (nt == -1) {
        continue;
      }
      if (!quad_is_in_corner_sides (&tempr, nt, &(info->sides))) {
        SC_CHECK_ABORT (!(!cside->is_ghost && iter_data->ghost_face),
                        "Iterate: quad missing corner neighbor");
      }
    }
#ifdef P4_TO_P8
    for (j = 0; j < 3; j++) {
      e = p8est_corner_edges[c][j];
      sc_array_resize (&quads, 0);
      sc_array_resize (&treeids, 0);
      p8est_quadrant_edge_neighbor_extra (&tempq, t, e, &quads, &treeids,
                                          NULL, conn);
      for (zz = 0; zz < quads.elem_count; zz++) {
        ptemp = p4est_quadrant_array_index (&quads, zz);
        nt = *((p4est_topidx_t *) sc_array_index (&treeids, zz));
        if (!quad_is_in_corner_sides (ptemp, nt, &(info->sides))) {
          SC_CHECK_ABORT (!(!cside->is_ghost && iter_data->ghost_edge),
                          "Iterate: quad missing corner neighbor");
        }
      }
    }
#endif
    sc_array_resize (&quads, 0);
    sc_array_resize (&treeids, 0);
    p4est_quadrant_corner_neighbor_extra (&tempq, t, c, &quads, &treeids,
                                          NULL, conn);
    for (zz = 0; zz < quads.elem_count; zz++) {
      ptemp = p4est_quadrant_array_index (&quads, zz);
      nt = *((p4est_topidx_t *) sc_array_index (&treeids, zz));
      if (!quad_is_in_corner_sides (ptemp, nt, &(info->sides))) {
        SC_CHECK_ABORT (!(iter_data->ghost_corner),
                        "Iterate: quad missing corner neighbor");
      }
    }

  }

  sc_array_reset (&quads);
  sc_array_reset (&treeids);
}

#ifdef P4_TO_P8
static              int8_t
test_edge_side (p4est_t * p4est, p8est_iter_edge_side_t * side, void *data)
{
  int                 i;
  int                 edge = (int) side->edge;
  int                 child_id, opp_id;
  int8_t              has_local = 0;
  p4est_quadrant_t   *q;
  p4est_quadrant_t    tempq;
  iter_data_t        *iter_data = (iter_data_t *) data;
  int                *checks = iter_data->checks;
  p4est_locidx_t      qid;
  p4est_topidx_t      t = side->treeid;
  p4est_tree_t       *tree;
  p4est_locidx_t      offset;
  int                 quad_count = 0;

  tree = p4est_tree_array_index (p4est->trees, t);
  offset = tree->quadrants_offset;

  if (!side->is_hanging) {
    q = side->is.full.quad;
    if (q == NULL) {
      SC_CHECK_ABORT (!(iter_data->ghost_edge),
                      "Iterate: full edge side missing quadrant");
      return 1;
    }
    if (!side->is.full.is_ghost && iter_data->count_edge) {
      qid = side->is.full.quadid + offset;
      checks[qid * checks_per_quad + edge + edge_offset]++;
    }
    return side->is.full.is_ghost;
  }
  else {
    for (i = 0; i < 2; i++) {
      q = side->is.hanging.quad[i];
      if (q == NULL) {
        SC_CHECK_ABORT (!(iter_data->ghost_corner),
                        "Iterate: hanging edge missing quadrant");
        continue;
      }
      quad_count++;
      child_id = p4est_quadrant_child_id (q);
      SC_CHECK_ABORT (p8est_edge_corners[edge][i] == child_id,
                      "Iterate: edge side ordering");
      opp_id = p8est_edge_corners[edge][1 - i];
      if (!side->is.hanging.is_ghost[i]) {
        has_local = 1;
        qid = side->is.hanging.quadid[i] + offset;
        if (iter_data->count_edge) {
          checks[qid * checks_per_quad + edge + edge_offset]++;
        }
        if (iter_data->count_corner) {
          checks[qid * checks_per_quad + opp_id + corner_offset]++;
        }
      }
    }
    switch (quad_count) {
    case 0:
      SC_CHECK_ABORT (!(iter_data->ghost_edge),
                      "Iterate: hanging edge missing all quadrants");
      break;
    case 1:
      SC_CHECK_ABORT (!(has_local && iter_data->ghost_face),
                      "Iterate: hanging edge missing quadrant");
      break;
    default:
      q = side->is.hanging.quad[0];
      p4est_quadrant_parent (q, &tempq);
      q = side->is.hanging.quad[1];
      SC_CHECK_ABORT (p4est_quadrant_is_parent (&tempq, q),
                      "Iterate: non siblings share edge");
      break;
    }
    return !has_local;
  }
}

static              int8_t
quad_is_in_edge_sides (p8est_quadrant_t * q, p4est_topidx_t t,
                       sc_array_t * sides)
{
  size_t              zz;
  size_t              zcount = sides->elem_count;
  p8est_iter_edge_side_t *eside;

  for (zz = 0; zz < zcount; zz++) {
    eside = p8est_iter_eside_array_index (sides, zz);
    if (eside->treeid == t) {
      if (!eside->is_hanging) {
        if (eside->is.full.quad != NULL) {
          if (p4est_quadrant_is_equal (q, eside->is.full.quad)) {
            return 1;
          }
        }
      }
      else if (eside->is.hanging.quad[0] != NULL) {
        if (p4est_quadrant_is_parent (q, eside->is.hanging.quad[0])) {
          return 1;
        }
      }
      else if (eside->is.hanging.quad[1] != NULL) {
        if (p4est_quadrant_is_parent (q, eside->is.hanging.quad[1])) {
          return 1;
        }
      }
    }
  }
  return 0;
}

static void
test_edge_boundary (p8est_iter_edge_info_t * info)
{
  int                 e, f, dir, which;
  int                 k, count;
  p4est_qcoord_t      end, xyz[P4EST_DIM];
  p8est_iter_edge_side_t *eside;
  p4est_quadrant_t   *q;

  SC_CHECK_ABORT (info->sides.elem_count > 0, "Empty edge iteration");
  eside = (p8est_iter_edge_side_t *) sc_array_index (&info->sides, 0);
  e = eside->edge;

  if (!info->tree_boundary) {
    SC_CHECK_ABORT (e & 1, "Not the lowest edge");
    return;
  }

  /* grab information about this quadrant */
  q = eside->is_hanging ? eside->is.hanging.quad[0] : eside->is.full.quad;
  if (q == NULL) {
    return;
  }
  end = P4EST_LAST_OFFSET (q->level);
  xyz[0] = q->x;
  xyz[1] = q->y;
  xyz[2] = q->z;

  /* check how many tree faces it is touching */
  count = 0;
  for (k = 0; k < 2; ++k) {
    f = p8est_edge_faces[e][k];
    dir = f >> 1;
    which = f & 1;
    count += (which == 0 && xyz[dir] == 0) || (which == 1 && xyz[dir] == end);
  }
  if (info->tree_boundary == P8EST_CONNECT_EDGE) {
    /* we are at a true inter-tree edge */
    SC_CHECK_ABORT (count == 2, "Not a tree boundary edge");
  }
  else if (info->tree_boundary == P4EST_CONNECT_FACE) {
    /* we are an edge inside an inter-tree face */
    SC_CHECK_ABORT (count == 1, "Not a tree face boundary edge");
  }
  else {
    SC_ABORT_NOT_REACHED ();
  }
}

static void
test_edge_adjacency (p8est_iter_edge_info_t * info, void *data)
{
  int                 i, j;
  int8_t              has_local = 0;
  int                 limit = (int) info->sides.elem_count;
  p8est_iter_edge_side_t *eside;
  p8est_quadrant_t    tempq, tempr;
  int                 e, f;
  p4est_topidx_t      t, nt;
  sc_array_t          quads, treeids;
  p8est_connectivity_t *conn = info->p4est->connectivity;
  size_t              zz;
  p8est_quadrant_t   *ptemp;
  iter_data_t        *iter_data = (iter_data_t *) data;

  test_edge_boundary (info);

  sc_array_init (&quads, sizeof (p8est_quadrant_t));
  sc_array_init (&treeids, sizeof (p4est_topidx_t));
  for (i = 0; i < limit; i++) {
    eside = p8est_iter_eside_array_index_int (&info->sides, i);
    has_local = (!test_edge_side (info->p4est, eside, data) || has_local);
  }
  SC_CHECK_ABORT (has_local, "Iterate: non local edge");

  for (i = 0; i < limit; i++) {
    eside = p8est_iter_eside_array_index_int (&(info->sides), i);
    if (!eside->is_hanging) {
      if (eside->is.full.quad == NULL) {
        continue;
      }
      tempq = *(eside->is.full.quad);
      has_local = !eside->is.full.is_ghost;
    }
    else {
      if (eside->is.hanging.quad[0] != NULL) {
        p4est_quadrant_parent (eside->is.hanging.quad[0], &tempq);
        has_local = !eside->is.hanging.is_ghost[0];
      }
      else {
        if (eside->is.hanging.quad[1] == NULL) {
          continue;
        }
        p4est_quadrant_parent (eside->is.hanging.quad[1], &tempq);
        has_local = !eside->is.hanging.is_ghost[1];
      }
    }
    e = (int) eside->edge;
    t = eside->treeid;
    for (j = 0; j < 2; j++) {
      f = p8est_edge_faces[e][j];
      nt = p8est_quadrant_face_neighbor_extra (&tempq, t, f, &tempr, NULL,
                                               conn);
      if (nt == -1) {
        continue;
      }
      if (!quad_is_in_edge_sides (&tempr, nt, &(info->sides))) {
        SC_CHECK_ABORT (!(has_local && iter_data->ghost_face),
                        "Iterate: quad missing edge neighbor");
      }
    }
    sc_array_resize (&quads, 0);
    sc_array_resize (&treeids, 0);
    p8est_quadrant_edge_neighbor_extra (&tempq, t, e, &quads, &treeids, NULL,
                                        conn);
    for (zz = 0; zz < quads.elem_count; zz++) {
      ptemp = p4est_quadrant_array_index (&quads, zz);
      nt = *((p4est_topidx_t *) sc_array_index (&treeids, zz));
      if (!quad_is_in_edge_sides (ptemp, nt, &(info->sides))) {
        SC_CHECK_ABORT (!(iter_data->ghost_edge),
                        "Iterate: quad missing edge neighbor");
      }
    }
    sc_array_reset (&quads);
    sc_array_reset (&treeids);
  }
}
#endif

static              int8_t
test_face_side (p4est_t * p4est, p4est_iter_face_side_t * side, void *data)
{
  int                 i;
  int                 face = (int) side->face;
  int                 child_id, opp_id;
  int8_t              has_local = 0;
  p4est_quadrant_t   *q;
  p4est_quadrant_t    tempq;
  iter_data_t        *iter_data = (iter_data_t *) data;
  int                *checks = iter_data->checks;
  p4est_locidx_t      qid;
#ifdef P4_TO_P8
  int                 edge;
  int                 dir = face / 2;
  int                 j;
#endif
  p4est_topidx_t      t = side->treeid;
  p4est_tree_t       *tree;
  p4est_locidx_t      offset;
  int                 quad_count = 0;

  tree = p4est_tree_array_index (p4est->trees, t);
  offset = tree->quadrants_offset;

  if (!side->is_hanging) {
    q = side->is.full.quad;
    if (q == NULL) {
      SC_CHECK_ABORT (!(iter_data->ghost_face),
                      "Iterate: full face side missing quadrant");
      return 1;
    }
    if (!side->is.full.is_ghost && iter_data->count_face) {
      qid = side->is.full.quadid + offset;
      checks[qid * checks_per_quad + face + face_offset]++;
    }
    return side->is.full.is_ghost;
  }
  else {
    for (i = 0; i < P4EST_CHILDREN / 2; i++) {
      q = side->is.hanging.quad[i];
      if (q == NULL) {
#ifndef P4_TO_P8
        SC_CHECK_ABORT (!(iter_data->ghost_face),
                        "Iterate: hanging face side missing quadrant");
#else
        SC_CHECK_ABORT (!(iter_data->ghost_edge),
                        "Iterate: hanging face side missing quadrant");
#endif
        continue;
      }
      quad_count++;
      child_id = p4est_quadrant_child_id (q);
      SC_CHECK_ABORT (p4est_face_corners[face][i] == child_id,
                      "Iterate: face side ordering");
      opp_id = p4est_face_corners[face][P4EST_CHILDREN / 2 - 1 - i];
      if (!side->is.hanging.is_ghost[i]) {
        has_local = 1;
        qid = side->is.hanging.quadid[i] + offset;
        if (iter_data->count_face) {
          checks[qid * checks_per_quad + face + face_offset]++;
        }
        if (iter_data->count_corner) {
          checks[qid * checks_per_quad + opp_id + corner_offset]++;
        }
#ifdef P4_TO_P8
        if (iter_data->count_edge) {
          for (j = 1; j <= 2; j++) {
            edge = p8est_corner_edges[opp_id][(dir + j) % 3];
            checks[qid * checks_per_quad + edge + edge_offset]++;
          }
        }
#endif
      }
    }
    switch (quad_count) {
    case 0:
    case 1:
      SC_CHECK_ABORT (!(iter_data->ghost_face),
                      "Iterate: hanging face side missing all quadrants");
      break;
#ifdef P4_TO_P8
    case 2:
      SC_CHECK_ABORT (!(iter_data->ghost_face),
                      "Iterate: hanging face side missing all quadrants");
      break;
    case 3:
      SC_CHECK_ABORT (!(iter_data->ghost_edge),
                      "Iterate: hanging face missing quadrant");
      break;
#endif
    default:
      q = side->is.hanging.quad[0];
#ifdef P4_TO_P8
      if (q == NULL) {
        q = side->is.hanging.quad[1];
      }
#endif
      p4est_quadrant_parent (q, &tempq);
      for (i = 1; i < P4EST_CHILDREN / 2; i++) {
        q = side->is.hanging.quad[i];
        if (q != NULL) {
          SC_CHECK_ABORT (p4est_quadrant_is_parent (&tempq, q),
                          "Iterate: non siblings share face");
        }
      }
      break;
    }
    return !has_local;
  }
}

static void
test_face_boundary (p4est_iter_face_info_t * info)
{
  int                 f;
  int                 result;
  p4est_qcoord_t      end;
  p4est_iter_face_side_t *fside;
  p4est_quadrant_t   *q;

  SC_CHECK_ABORT (info->sides.elem_count > 0, "Empty face iteration");
  fside = (p4est_iter_face_side_t *) sc_array_index (&info->sides, 0);
  f = fside->face;

  if (!info->tree_boundary) {
    SC_CHECK_ABORT (f & 1, "Not the lowest face");
  }
  else {
    SC_CHECK_ABORT (info->tree_boundary == P4EST_CONNECT_FACE, "Not a face");
    q = fside->is_hanging ? fside->is.hanging.quad[0] : fside->is.full.quad;
    if (q == NULL) {
      return;
    }
    end = P4EST_LAST_OFFSET (q->level);

    result =
      (f == 0 && q->x == 0) || (f == 1 && q->x == end) ||
      (f == 2 && q->y == 0) || (f == 3 && q->y == end) ||
#ifdef P4_TO_P8
      (f == 4 && q->z == 0) || (f == 5 && q->z == end) ||
#endif
      0;
    SC_CHECK_ABORT (result, "Not a tree boundary face");
  }
}

static void
test_face_adjacency (p4est_iter_face_info_t * info, void *data)
{
  int                 i, j;
  int8_t              is_ghost[2];
  int                 limit = (int) info->sides.elem_count;
  p4est_iter_face_side_t *fside;
  p4est_quadrant_t    tempq[2], tempr[2];
  int                 face[2];
  p4est_topidx_t      treeid[2], nt[2];

  test_face_boundary (info);

  for (i = 0; i < 2; i++) {
    is_ghost[i] = 1;
  }
  for (i = 0; i < limit; i++) {
    fside = p4est_iter_fside_array_index_int (&info->sides, i);
    is_ghost[i] = test_face_side (info->p4est, fside, data);
  }
  SC_CHECK_ABORT (!(is_ghost[0] && is_ghost[1]), "Iterate: non local face");

  for (i = 0; i < limit; i++) {
    fside = p4est_iter_fside_array_index_int (&info->sides, i);
    face[i] = (int) fside->face;
    treeid[i] = fside->treeid;
    if (!fside->is_hanging) {
      if (fside->is.full.quad == NULL) {
        return;
      }
      tempq[i] = *(fside->is.full.quad);
    }
    else {
      for (j = 0; j < P4EST_CHILDREN / 2; j++) {
        if (fside->is.hanging.quad[0] != NULL) {
          p4est_quadrant_parent (fside->is.hanging.quad[0], &(tempq[i]));
          break;
        }
      }
      if (j == P4EST_CHILDREN / 2) {
        return;
      }
    }
    nt[1 - i] =
      p4est_quadrant_face_neighbor_extra (&(tempq[i]), treeid[i], face[i],
                                          &(tempr[1 - i]), NULL,
                                          info->p4est->connectivity);
  }
  if (limit == 2) {
    for (i = 0; i < 2; i++) {
      SC_CHECK_ABORT (nt[i] == treeid[i], "Iterate: face tree mismatch");
      SC_CHECK_ABORT (p4est_quadrant_is_equal (&(tempq[i]), &(tempr[i])),
                      "Iterate: face neighbor mismatch");
    }
  }
  else {
    SC_CHECK_ABORT (nt[1] == -1,
                    "Iterate: non boundary face without neighbor");
    fside = p4est_iter_fside_array_index_int (&info->sides, 0);
    SC_CHECK_ABORT (!fside->is_hanging,
                    "Iterate: hanging face without neighbor");
  }
}

static void
test_volume_adjacency (p4est_iter_volume_info_t * info, void *data)
{
  iter_data_t        *iter_data = (iter_data_t *) data;
  int                *checks = iter_data->checks;
  p4est_locidx_t      qid = info->quadid;
  p4est_topidx_t      t = info->treeid;
  p4est_tree_t       *tree = p4est_tree_array_index (info->p4est->trees, t);

  SC_CHECK_ABORT (info->quad != NULL, "Iterate: missing volume quad");
  if (iter_data->count_volume) {
    qid += tree->quadrants_offset;
    checks[qid * checks_per_quad]++;
  }
}

int
main (int argc, char **argv)
{
  sc_MPI_Comm         mpicomm;
  int                 mpiret;
  int                 mpisize, mpirank;
  p4est_t            *p4est;
  p4est_connectivity_t *connectivity;
  p4est_locidx_t      num_quads, li;
  p4est_locidx_t      num_checks;
  int                *checks;
  p4est_ghost_t      *ghost_layer;
  int                 ntests;
  int                 i, j, k;
  iter_data_t         iter_data;
  p4est_iter_volume_t iter_volume;
  p4est_iter_face_t   iter_face;
#ifdef P4_TO_P8
  p8est_iter_edge_t   iter_edge;
#endif
  p4est_iter_corner_t iter_corner;
  int                 volume_count;
  int                 face_count;
#ifdef P4_TO_P8
  int                 edge_count;
#endif
  int                 corner_count;

  ntests = 3;

  /* initialize MPI */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpicomm = sc_MPI_COMM_WORLD;
  mpiret = sc_MPI_Comm_size (mpicomm, &mpisize);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (mpicomm, &mpirank);
  SC_CHECK_MPI (mpiret);

  sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  for (i = 0; i < ntests; i++) {
    /* create connectivity and forest structures */
    switch (i) {
#ifndef P4_TO_P8
    case 0:
      connectivity = p4est_connectivity_new_moebius ();
      break;
    case 1:
      connectivity = p4est_connectivity_new_star ();
      break;
    default:
      connectivity = p4est_connectivity_new_periodic ();
      break;
#else
    case 0:
      connectivity = p8est_connectivity_new_periodic ();
      break;
    case 1:
      connectivity = p8est_connectivity_new_rotwrap ();
      break;
    default:
      connectivity = p8est_connectivity_new_rotcubes ();
      break;
#endif
    }
    p4est = p4est_new_ext (mpicomm, connectivity, 15, 0, 0, 0, NULL, NULL);

    /* refine to make the number of elements interesting */
    p4est_refine (p4est, 1, refine_fn, NULL);

    /* balance the forest */
    /* TODO: use BALANCE_FACE/EDGE when that is known to work */
    p4est_balance (p4est, P4EST_CONNECT_FULL, NULL);

    /* do a uniform partition */
    p4est_partition (p4est, 0, NULL);

    num_quads = p4est->local_num_quadrants;
    num_checks = checks_per_quad * num_quads;
    checks = P4EST_ALLOC_ZERO (int, num_checks);

    iter_data.checks = checks;

    volume_count = 0;
    face_count = 0;
#ifdef P4_TO_P8
    edge_count = 0;
#endif
    corner_count = 0;

    for (j = 0; j <= P4EST_DIM; j++) {

      iter_data.count_volume = 1;
      iter_data.count_face = (j > 0);
#ifdef P4_TO_P8
      iter_data.count_edge = (j > 1);
#endif
      iter_data.count_corner = (j == P4EST_DIM);

      for (k = 0; k <= P4EST_DIM; k++) {
        switch (k) {
        case 0:
          ghost_layer = NULL;
          iter_data.ghost_face = 0;
#ifdef P4_TO_P8
          iter_data.ghost_edge = 0;
#endif
          iter_data.ghost_corner = 0;
          break;
#ifndef P4_TO_P8
        case 1:
          ghost_layer = p4est_ghost_new (p4est, P4EST_CONNECT_FACE);
          iter_data.ghost_face = 1;
          iter_data.ghost_corner = 0;
          break;
        default:
          ghost_layer = p4est_ghost_new (p4est, P4EST_CONNECT_CORNER);
          iter_data.ghost_face = 1;
          iter_data.ghost_corner = 1;
          break;
#else
        case 1:
          ghost_layer = p4est_ghost_new (p4est, P8EST_CONNECT_FACE);
          iter_data.ghost_face = 1;
          iter_data.ghost_edge = 0;
          iter_data.ghost_corner = 0;
          break;
        case 2:
          ghost_layer = p4est_ghost_new (p4est, P8EST_CONNECT_EDGE);
          iter_data.ghost_face = 1;
          iter_data.ghost_edge = 1;
          iter_data.ghost_corner = 0;
          break;
        default:
          ghost_layer = p4est_ghost_new (p4est, P8EST_CONNECT_CORNER);
          iter_data.ghost_face = 1;
          iter_data.ghost_edge = 1;
          iter_data.ghost_corner = 1;
          break;
#endif
        }

        if (iter_data.count_volume) {
          iter_volume = test_volume_adjacency;
          volume_count++;
        }
        else {
          iter_volume = NULL;
        }
        if (iter_data.count_face) {
          iter_face = test_face_adjacency;
          face_count++;
        }
        else {
          iter_face = NULL;
        }
#ifdef P4_TO_P8
        if (iter_data.count_edge) {
          iter_edge = test_edge_adjacency;
          edge_count++;
        }
        else {
          iter_edge = NULL;
        }
#endif
        if (iter_data.count_corner) {
          iter_corner = test_corner_adjacency;
          corner_count++;
        }
        else {
          iter_corner = NULL;
        }

        P4EST_GLOBAL_PRODUCTIONF ("Begin adjacency test %d:%d:%d\n", i, j, k);

        p4est_iterate (p4est, ghost_layer, &iter_data, iter_volume, iter_face,
#ifdef P4_TO_P8
                       iter_edge,
#endif
                       iter_corner);

        for (li = 0; li < num_checks; li++) {
          switch (check_to_type[li % checks_per_quad]) {
          case P4EST_DIM:
            SC_CHECK_ABORT (checks[li] == volume_count,
                            "Iterate: completion check");
            break;
          case (P4EST_DIM - 1):
            SC_CHECK_ABORT (checks[li] == face_count,
                            "Iterate: completion check");
            break;
#ifdef P4_TO_P8
          case 1:
            SC_CHECK_ABORT (checks[li] == edge_count,
                            "Iterate: completion check");
            break;
#endif
          default:
            SC_CHECK_ABORT (checks[li] == corner_count,
                            "Iterate: completion check");
          }
        }
        /* clean up */
        if (k > 0) {
          p4est_ghost_destroy (ghost_layer);
        }
      }
    }
    P4EST_FREE (checks);

    p4est_destroy (p4est);
    p4est_connectivity_destroy (connectivity);
    P4EST_GLOBAL_PRODUCTIONF ("End adjacency test %d\n", i);
  }

  /* exit */
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}

/* EOF test_iterate2.c */
