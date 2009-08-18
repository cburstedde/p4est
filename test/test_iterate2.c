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

#ifndef P4_TO_P8
#include <p4est_bits.h>
#include <p4est_algorithms.h>
#include <p4est_ghost.h>
#include <p4est_iterate.h>
#else
#include <p8est_bits.h>
#include <p8est_algorithms.h>
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
  bool                count_volume;
  bool                count_face;
  bool                ghost_face;
#ifdef P4_TO_P8
  bool                count_edge;
  bool                ghost_edge;
#endif
  bool                count_corner;
  bool                ghost_corner;
  int                *checks;
}
iter_data_t;

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

static              bool
test_corner_side (p4est_t * p4est, p4est_iter_corner_side_t * side,
                  void *data)
{
  int                 corner = side->corner;
  p4est_quadrant_t   *q = side->quad;
  iter_data_t        *iter_data = (iter_data_t *) data;
  int                *checks = iter_data->checks;
  p4est_locidx_t      qid = side->quadid;
  p4est_topidx_t      t = side->treeid;
  p4est_tree_t       *tree;
  p4est_locidx_t      offset;

  tree = p4est_array_index_topidx (p4est->trees, t);
  offset = tree->quadrants_offset;

  q = side->quad;
  if (q == NULL) {
    SC_CHECK_ABORT (!(iter_data->ghost_corner), "Iterate: empty corner side");
    return false;
  }
  else {
    if (side->is_local && iter_data->count_corner) {
      qid += offset;
      checks[qid * checks_per_quad + corner + corner_offset]++;
    }
  }

  return side->is_local;
}

static              bool
quad_is_in_corner_sides (p4est_quadrant_t * q, p4est_topidx_t t,
                         sc_array_t * sides)
{
  size_t              zz;
  size_t              zcount = sides->elem_count;
  p4est_iter_corner_side_t *cside;
  p4est_quadrant_t   *r;

  for (zz = 0; zz < zcount; zz++) {
    cside = sc_array_index (sides, zz);
    if (cside->treeid == t) {
      r = cside->quad;
      if (r == NULL) {
        continue;
      }
      if (p4est_quadrant_is_equal (q, r)) {
        return true;
      }
      if (p4est_quadrant_is_ancestor (q, r)) {
        return true;
      }
    }
  }

  return false;
}

static void
test_corner_adjacency (p4est_iter_corner_info_t * info, void *data)
{
  int                 i, j;
  bool                has_local = false;
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

  for (i = 0; i < limit; i++) {
    cside = sc_array_index_int (&(info->sides), i);
    has_local = (test_corner_side (info->p4est, cside, data) || has_local);
    if (cside->quad != NULL) {
      min_level = (cside->quad->level < min_level) ? cside->quad->level :
        min_level;
    }
  }
  SC_CHECK_ABORT (has_local, "Iterate: non local corner");

  sc_array_init (&quads, sizeof (p4est_quadrant_t));
  sc_array_init (&treeids, sizeof (p4est_topidx_t));
  for (i = 0; i < limit; i++) {
    cside = sc_array_index_int (&(info->sides), i);
    if (cside->quad == NULL) {
      continue;
    }
    c = cside->corner;
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
      nt = p4est_quadrant_face_neighbor_extra (&tempq, t, f, &tempr, conn);
      if (nt == -1) {
        continue;
      }
      if (!quad_is_in_corner_sides (&tempr, nt, &(info->sides))) {
        SC_CHECK_ABORT (!(cside->is_local && iter_data->ghost_face),
                        "Iterate: quad missing corner neighbor");
      }
    }
#ifdef P4_TO_P8
    for (j = 0; j < 3; j++) {
      e = p8est_corner_edges[c][j];
      sc_array_resize (&quads, 0);
      sc_array_resize (&treeids, 0);
      p8est_quadrant_edge_neighbor_extra (&tempq, t, e, &quads, &treeids,
                                          conn);
      for (zz = 0; zz < quads.elem_count; zz++) {
        ptemp = sc_array_index (&quads, zz);
        nt = *((p4est_topidx_t *) sc_array_index (&treeids, zz));
        if (!quad_is_in_corner_sides (ptemp, nt, &(info->sides))) {
          SC_CHECK_ABORT (!(cside->is_local && iter_data->ghost_edge),
                          "Iterate: quad missing corner neighbor");
        }
      }
    }
#endif
    sc_array_resize (&quads, 0);
    sc_array_resize (&treeids, 0);
    p4est_quadrant_corner_neighbor_extra (&tempq, t, c, &quads, &treeids,
                                          conn);
    for (zz = 0; zz < quads.elem_count; zz++) {
      ptemp = sc_array_index (&quads, zz);
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
static              bool
test_edge_side (p4est_t * p4est, p8est_iter_edge_side_t * side, void *data)
{
  int                 i;
  int                 edge = side->edge;
  int                 child_id, opp_id;
  bool                has_local = false;
  p4est_quadrant_t   *q;
  p4est_quadrant_t    tempq;
  iter_data_t        *iter_data = (iter_data_t *) data;
  int                *checks = iter_data->checks;
  p4est_locidx_t      qid;
  p4est_topidx_t      t = side->treeid;
  p4est_tree_t       *tree;
  p4est_locidx_t      offset;
  int                 quad_count = 0;

  tree = p4est_array_index_topidx (p4est->trees, t);
  offset = tree->quadrants_offset;

  if (!side->is_hanging) {
    q = side->is.full.quad;
    if (q == NULL) {
      SC_CHECK_ABORT (!(iter_data->ghost_edge),
                      "Iterate: full edge side missing quadrant");
      return false;
    }
    if (side->is.full.is_local && iter_data->count_edge) {
      qid = side->is.full.quadid + offset;
      checks[qid * checks_per_quad + edge + edge_offset]++;
    }
    return side->is.full.is_local;
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
      if (side->is.hanging.is_local[i]) {
        has_local = true;
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
    return has_local;
  }
}

static              bool
quad_is_in_edge_sides (p8est_quadrant_t * q, p4est_topidx_t t,
                       sc_array_t * sides)
{
  size_t              zz;
  size_t              zcount = sides->elem_count;
  p8est_iter_edge_side_t *eside;

  for (zz = 0; zz < zcount; zz++) {
    eside = sc_array_index (sides, zz);
    if (eside->treeid == t) {
      if (!eside->is_hanging) {
        if (eside->is.full.quad != NULL) {
          if (p4est_quadrant_is_equal (q, eside->is.full.quad)) {
            return true;
          }
        }
      }
      else if (eside->is.hanging.quad[0] != NULL) {
        if (p4est_quadrant_is_parent (q, eside->is.hanging.quad[0])) {
          return true;
        }
      }
      else if (eside->is.hanging.quad[1] != NULL) {
        if (p4est_quadrant_is_parent (q, eside->is.hanging.quad[1])) {
          return true;
        }
      }
    }
  }
  return false;
}

static void
test_edge_adjacency (p8est_iter_edge_info_t * info, void *data)
{
  int                 i, j;
  bool                has_local = false;
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

  sc_array_init (&quads, sizeof (p8est_quadrant_t));
  sc_array_init (&treeids, sizeof (p4est_topidx_t));
  for (i = 0; i < limit; i++) {
    eside = sc_array_index_int (&(info->sides), i);
    has_local = (test_edge_side (info->p4est, eside, data) || has_local);
  }
  SC_CHECK_ABORT (has_local, "Iterate: non local edge");

  for (i = 0; i < limit; i++) {
    eside = sc_array_index_int (&(info->sides), i);
    if (!eside->is_hanging) {
      if (eside->is.full.quad == NULL) {
        continue;
      }
      tempq = *(eside->is.full.quad);
      has_local = eside->is.full.is_local;
    }
    else {
      if (eside->is.hanging.quad[0] != NULL) {
        p4est_quadrant_parent (eside->is.hanging.quad[0], &tempq);
        has_local = eside->is.hanging.is_local[0];
      }
      else {
        if (eside->is.hanging.quad[1] == NULL) {
          continue;
        }
        p4est_quadrant_parent (eside->is.hanging.quad[1], &tempq);
        has_local = eside->is.hanging.is_local[1];
      }
    }
    e = eside->edge;
    t = eside->treeid;
    for (j = 0; j < 2; j++) {
      f = p8est_edge_faces[e][j];
      nt = p8est_quadrant_face_neighbor_extra (&tempq, t, f, &tempr, conn);
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
    p8est_quadrant_edge_neighbor_extra (&tempq, t, e, &quads, &treeids, conn);
    for (zz = 0; zz < quads.elem_count; zz++) {
      ptemp = sc_array_index (&quads, zz);
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

static              bool
test_face_side (p4est_t * p4est, p4est_iter_face_side_t * side, void *data)
{
  int                 i;
  int                 face = side->face;
  int                 child_id, opp_id;
  bool                has_local = false;
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

  tree = p4est_array_index_topidx (p4est->trees, t);
  offset = tree->quadrants_offset;

#ifndef P4_TO_P8
  face = p4est_zface_to_rface[face];
#endif

  if (!side->is_hanging) {
    q = side->is.full.quad;
    if (q == NULL) {
      SC_CHECK_ABORT (!(iter_data->ghost_face),
                      "Iterate: full face side missing quadrant");
      return false;
    }
    if (side->is.full.is_local && iter_data->count_face) {
      qid = side->is.full.quadid + offset;
      checks[qid * checks_per_quad + face + face_offset]++;
    }
    return side->is.full.is_local;
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
      if (side->is.hanging.is_local[i]) {
        has_local = true;
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
    return has_local;
  }
}

static void
test_face_adjacency (p4est_iter_face_info_t * info, void *data)
{
  int                 i, j;
  bool                is_local[2];
  int                 limit = (int) info->sides.elem_count;
  p4est_iter_face_side_t *fside;
  p4est_quadrant_t    tempq[2], tempr[2];
  int                 face[2];
  p4est_topidx_t      treeid[2], nt[2];

  for (i = 0; i < 2; i++) {
    is_local[i] = false;
  }
  for (i = 0; i < limit; i++) {
    fside = sc_array_index_int (&(info->sides), i);
    is_local[i] = test_face_side (info->p4est, fside, data);
  }
  SC_CHECK_ABORT (is_local[0] || is_local[1], "Iterate: non local face");

  for (i = 0; i < limit; i++) {
    fside = sc_array_index_int (&(info->sides), i);
    face[i] = fside->face;
#ifndef P4_TO_P8
    face[i] = p4est_zface_to_rface[face[i]];
#endif
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
    nt[1 - i] = p4est_quadrant_face_neighbor_extra (&(tempq[i]), treeid[i],
                                                    face[i], &(tempr[1 - i]),
                                                    info->p4est->
                                                    connectivity);
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
    fside = sc_array_index_int (&(info->sides), 0);
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
  p4est_tree_t       *tree = p4est_array_index_topidx (info->p4est->trees, t);

  SC_CHECK_ABORT (info->quad != NULL, "Iterate: missing volume quad");
  if (iter_data->count_volume) {
    qid += tree->quadrants_offset;
    checks[qid * checks_per_quad]++;
  }
}

int
main (int argc, char **argv)
{
  MPI_Comm            mpicomm;
  int                 mpiret;
  int                 mpisize, mpirank;
  p4est_t            *p4est;
  p4est_connectivity_t *connectivity;
  p4est_locidx_t      num_quads, li;
  p4est_locidx_t      num_checks;
  int                *checks;
  sc_array_t          Ghost_layer;
  sc_array_t         *ghost_layer;
  int                *ghost_owner;
  bool                success;
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

#ifndef P4_TO_P8
  ntests = 3;
#else
  ntests = 4;
#endif

  /* initialize MPI */
  mpiret = MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpicomm = MPI_COMM_WORLD;
  mpiret = MPI_Comm_size (mpicomm, &mpisize);
  SC_CHECK_MPI (mpiret);
  mpiret = MPI_Comm_rank (mpicomm, &mpirank);
  SC_CHECK_MPI (mpiret);

  sc_init (mpicomm, true, true, NULL, SC_LP_DEFAULT);
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
    case 2:
      connectivity = p8est_connectivity_new_rotcubes ();
      break;
    default:
      connectivity = p8est_connectivity_new_shell ();
      break;
#endif
    }
    p4est = p4est_new (mpicomm, connectivity, 15, 0, NULL, NULL);

    /* refine to make the number of elements interesting */
    p4est_refine (p4est, true, refine_fn, NULL);

    /* balance the forest */
    p4est_balance (p4est, P4EST_BALANCE_DEFAULT, NULL);

    /* do a uniform partition */
    p4est_partition (p4est, NULL);

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

      iter_data.count_volume = true;
      iter_data.count_face = (j > 0);
#ifdef P4_TO_P8
      iter_data.count_edge = (j > 1);
#endif
      iter_data.count_corner = (j == P4EST_DIM);

      for (k = 0; k <= P4EST_DIM; k++) {
        sc_array_init (&Ghost_layer, sizeof (p4est_quadrant_t));
        switch (k) {
        case 0:
          ghost_layer = NULL;
          iter_data.ghost_face = false;
#ifdef P4_TO_P8
          iter_data.ghost_edge = false;
#endif
          iter_data.ghost_corner = false;
          break;
#ifndef P4_TO_P8
        case 1:
          ghost_layer = &Ghost_layer;
          success = p4est_build_ghost_layer (p4est, P4EST_BALANCE_FACE,
                                             ghost_layer, &ghost_owner);
          P4EST_ASSERT (success);
          iter_data.ghost_face = true;
          iter_data.ghost_corner = false;
          break;
        default:
          ghost_layer = &Ghost_layer;
          success = p4est_build_ghost_layer (p4est, P4EST_BALANCE_CORNER,
                                             ghost_layer, &ghost_owner);
          P4EST_ASSERT (success);
          iter_data.ghost_face = true;
          iter_data.ghost_corner = true;
          break;
#else
        case 1:
          ghost_layer = &Ghost_layer;
          success = p4est_build_ghost_layer (p4est, P8EST_BALANCE_FACE,
                                             ghost_layer, &ghost_owner);
          P4EST_ASSERT (success);
          iter_data.ghost_face = true;
          iter_data.ghost_edge = false;
          iter_data.ghost_corner = false;
          break;
        case 2:
          ghost_layer = &Ghost_layer;
          success = p4est_build_ghost_layer (p4est, P8EST_BALANCE_EDGE,
                                             ghost_layer, &ghost_owner);
          P4EST_ASSERT (success);
          iter_data.ghost_face = true;
          iter_data.ghost_edge = true;
          iter_data.ghost_corner = false;
          break;
        default:
          ghost_layer = &Ghost_layer;
          success = p4est_build_ghost_layer (p4est, P8EST_BALANCE_CORNER,
                                             ghost_layer, &ghost_owner);
          P4EST_ASSERT (success);
          iter_data.ghost_face = true;
          iter_data.ghost_edge = true;
          iter_data.ghost_corner = true;
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
        sc_array_reset (&Ghost_layer);
        if (k > 0) {
          P4EST_FREE (ghost_owner);
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

  mpiret = MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}

/* EOF test_iterate2.c */
