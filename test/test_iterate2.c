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

static void
face_do_nothing (p4est_iter_face_info_t * info, void *data)
{
};

static void
corner_do_nothing (p4est_iter_corner_info_t * info, void *data)
{
};

static int
quad_is_in_array (p4est_t * p4est, p4est_quadrant_t * q, p4est_topidx_t t,
                  sc_array_t * quads, sc_array_t * treeids,
                  int *level_diff, int *child_id)
{
  int                 i, n = (int) quads->elem_count;
  p4est_topidx_t      nt;
  p4est_quadrant_t   *p;
  P4EST_ASSERT (treeids->elem_count == (size_t) n);
  for (i = 0; i < n; i++) {
    nt = *((p4est_topidx_t *) sc_array_index_int (treeids, i));
    if (nt != t) {
      continue;
    }
    p = *((p4est_quadrant_t **) sc_array_index_int (quads, i));
    if (p4est_quadrant_is_equal (q, p)) {
      *level_diff = 0;
      *child_id = -1;
      return i;
    }
    if (p4est_quadrant_is_ancestor (q, p)) {
      *level_diff = (int) (q->level - p->level);
      *child_id = p4est_quadrant_child_id (p);
      return i;
    }
    if (p4est_quadrant_is_ancestor (p, q)) {
      *level_diff = (int) (q->level - p->level);
      *child_id = p4est_quadrant_child_id (q);
      return i;
    }
  }
  return -1;
}

static              p4est_topidx_t
quad_face_neighbor (p4est_t * p4est, p4est_quadrant_t * q, p4est_topidx_t t,
                    int f, p4est_quadrant_t * p)
{
  p4est_connectivity_t *conn = p4est->connectivity;
#ifndef P4_TO_P8
  int                 transform;
#else
  int                 transform[9];
  p4est_topidx_t      flag;
#endif
  p4est_quadrant_t    temp;

  p4est_quadrant_face_neighbor (q, f, p);
  if (p4est_quadrant_is_inside_root (p)) {
    return t;
  }
#ifndef P4_TO_P8
  transform = p4est_find_face_transform (conn, t, f);
  if (transform == -1) {
    *p = *q;
    return t;
  }
  p4est_quadrant_translate_face (p, f);
  p4est_quadrant_transform_face (p, &temp, transform);
#else
  flag = p8est_find_face_transform (conn, t, f, transform);
  if (flag == -1) {
    *p = *q;
    return t;
  }
  p8est_quadrant_transform_face (p, &temp, transform);
#endif
  *p = temp;

  return conn->tree_to_tree[t * 2 * P4EST_DIM + f];
}

#ifdef P4_TO_P8
static              p4est_topidx_t
quad_edge_neighbor (p4est_t * p4est, p4est_quadrant_t * q, p4est_topidx_t t,
                    int e, p4est_quadrant_t * p)
{
  int                 f;
  p4est_quadrant_t    temp;
  p8est_quadrant_edge_neighbor (q, e, p);
  if (p4est_quadrant_is_inside_root (p)) {
    return t;
  }
  P4EST_ASSERT (p4est_quadrant_is_outside_face (p));
  f = p8est_edge_faces[e][0];
  p4est_quadrant_face_neighbor (q, f, &temp);
  if (p4est_quadrant_is_inside_root (&temp)) {
    return quad_face_neighbor (p4est, &temp, t, p8est_edge_faces[e][1], p);
  }
  f = p8est_edge_faces[e][1];
  p4est_quadrant_face_neighbor (q, f, &temp);
  P4EST_ASSERT (p4est_quadrant_is_inside_root (&temp));

  return quad_face_neighbor (p4est, &temp, t, p8est_edge_faces[e][0], p);
}
#endif

static              p4est_topidx_t
quad_corner_neighbor (p4est_t * p4est, p4est_quadrant_t * q, p4est_topidx_t t,
                      int c, p4est_quadrant_t * p)
{
  int                 f;
  p4est_quadrant_t    temp;
  p4est_quadrant_corner_neighbor (q, c, p);
  if (p4est_quadrant_is_inside_root (p)) {
    return t;
  }
  P4EST_ASSERT (p4est_quadrant_is_outside_face (p));
#ifndef P4_TO_P8
  f = p4est_corner_faces[c][0];
  p4est_quadrant_face_neighbor (q, f, &temp);
  if (p4est_quadrant_is_inside_root (&temp)) {
    return quad_face_neighbor (p4est, &temp, t, p4est_corner_faces[c][1], p);
  }
  f = p4est_corner_faces[c][1];
  p4est_quadrant_face_neighbor (q, f, &temp);
  P4EST_ASSERT (p4est_quadrant_is_inside_root (&temp));

  return quad_face_neighbor (p4est, &temp, t, p4est_corner_faces[c][0], p);
#else
  f = p8est_corner_faces[c][0];
  p4est_quadrant_face_neighbor (q, f, &temp);
  if (p4est_quadrant_is_inside_root (&temp)) {
    return quad_edge_neighbor (p4est, &temp, t, p8est_corner_edges[c][0], p);
  }
  f = p8est_corner_faces[c][1];
  p4est_quadrant_face_neighbor (q, f, &temp);
  if (p4est_quadrant_is_inside_root (&temp)) {
    return quad_edge_neighbor (p4est, &temp, t, p8est_corner_edges[c][1], p);
  }
  f = p8est_corner_faces[c][2];
  p4est_quadrant_face_neighbor (q, f, &temp);
  P4EST_ASSERT (p4est_quadrant_is_inside_root (&temp));

  return quad_edge_neighbor (p4est, &temp, t, p8est_corner_edges[c][2], p);
#endif
}

static void
corner_test_adjacency (p4est_iter_corner_info_t * info, void *data)
{
  p4est_t            *p4est = info->p4est;
  int                 i, j, k;
  int                 n = (int) info->quads->elem_count;
  int                 c;
  int                 child_id, level_diff;
  p4est_quadrant_t   *q, temp;
  p4est_topidx_t      t, nt;
  sc_array_t         *corners = info->corners;
  p4est_connectivity_t *conn = p4est->connectivity;

  size_t              ctree;
#ifndef P4_TO_P8
  p4est_corner_transform_t *ct;
  sc_array_t          ctransforms, *cta;
  cta = &ctransforms;
#else
  p4est_quadrant_t    temp2;
  p8est_corner_info_t ci;
  p8est_corner_transform_t *ct;
  sc_array_t         *cta = &ci.corner_transforms;

  p8est_edge_info_t   ei;
  p8est_edge_transform_t *et;
  sc_array_t         *eta = &ei.edge_transforms;
  size_t              etree;
  int                 e;
#endif

  if (n == 1) {
    c = *((int *) sc_array_index (info->corners, 0));
    q = *((p4est_quadrant_t **) sc_array_index (info->quads, 0));
    SC_CHECK_ABORT (p4est_quadrant_touches_corner (q, c, true),
                    "Iterate: mismatched corner 1");
    return;
  }
#ifndef P4_TO_P8
  sc_array_init (cta, sizeof (p4est_corner_transform_t));
#else
  sc_array_init (eta, sizeof (p8est_edge_transform_t));
  sc_array_init (cta, sizeof (p8est_corner_transform_t));
#endif /* !P4_TO_P8 */
  for (i = 0; i < n; i++) {
    t = *((p4est_topidx_t *) sc_array_index_int (info->treeids, i));
    q = *((p4est_quadrant_t **) sc_array_index_int (info->quads, i));
    c = *((int *) sc_array_index_int (corners, i));
    for (j = 0; j < P4EST_DIM; j++) {
      nt = quad_face_neighbor (p4est, q, t, p4est_corner_faces[c][j], &temp);
      k = quad_is_in_array (p4est, &temp, nt, info->quads, info->treeids,
                            &level_diff, &child_id);
      SC_CHECK_ABORT (k >= 0, "Iterate: mismatched corner 2");
      SC_CHECK_ABORT (-1 <= level_diff && level_diff <= 1,
                      "Iterate: mismatched corner 3");
      if (level_diff != 0) {
        SC_CHECK_ABORT
          (child_id == *((int *) sc_array_index_int (corners, k)),
           "Iterate: mismatched corner 4");
      }

#ifdef P4_TO_P8
      e = p8est_corner_edges[c][j];
      if (!p8est_quadrant_touches_edge (q, e, true)) {
        nt = quad_edge_neighbor (p4est, q, t, e, &temp);
        k = quad_is_in_array (p4est, &temp, nt, info->quads, info->treeids,
                              &level_diff, &child_id);
        SC_CHECK_ABORT (k >= 0, "Iterate: mismatched corner 5");
        SC_CHECK_ABORT (-1 <= level_diff && level_diff <= 1,
                        "Iterate: mismatched corner 6");
        if (level_diff != 0) {
          SC_CHECK_ABORT
            (child_id == *((int *) sc_array_index_int (corners, k)),
             "Iterate: mismatched corner 7");
        }
      }
      else {
        p8est_quadrant_edge_neighbor (q, e, &temp2);
        p8est_find_edge_transform (conn, t, e, &ei);
        for (etree = 0; etree < eta->elem_count; etree++) {
          et = sc_array_index (eta, etree);
          p8est_quadrant_transform_edge (&temp2, &temp, &ei, et, true);
          nt = et->ntree;
          k = quad_is_in_array (p4est, &temp, nt, info->quads, info->treeids,
                                &level_diff, &child_id);
          SC_CHECK_ABORT (k >= 0, "Iterate: mismatched corner 8");
          SC_CHECK_ABORT (-1 <= level_diff && level_diff <= 1,
                          "Iterate: mismatched corner 9");
          if (level_diff != 0) {
            SC_CHECK_ABORT
              (child_id == *((int *) sc_array_index_int (corners, k)),
               "Iterate: mismatched corner 10");
          }
        }
      }
#endif
    }

    if (!p4est_quadrant_touches_corner (q, c, true)) {
#ifndef P4_TO_P8
      nt = quad_corner_neighbor (p4est, q, t, c, &temp);
      k = quad_is_in_array (p4est, &temp, nt, info->quads, info->treeids,
                            &level_diff, &child_id);
      SC_CHECK_ABORT (k >= 0, "Iterate: mismatched corner 11");
      if (level_diff != 0) {
        SC_CHECK_ABORT
          (child_id == *((int *) sc_array_index_int (corners, k)),
           "Iterate: mismatched corner 12");
      }
#else
      for (j = 0; j < P4EST_DIM; j++) {
        e = p8est_corner_edges[c][j];
        if (p8est_quadrant_touches_edge (q, e, true)) {
          break;
        }
      }
      if (j == P4EST_DIM) {
        nt = quad_corner_neighbor (p4est, q, t, c, &temp);
        k = quad_is_in_array (p4est, &temp, nt, info->quads, info->treeids,
                              &level_diff, &child_id);
        SC_CHECK_ABORT (k >= 0, "Iterate: mismatched corner 13");
        if (level_diff != 0) {
          SC_CHECK_ABORT
            (child_id == *((int *) sc_array_index_int (corners, k)),
             "Iterate: mismatched corner 14");
        }
      }
      else {
        p4est_quadrant_face_neighbor (q, p8est_corner_faces[c][j], &temp);
        p8est_quadrant_edge_neighbor (&temp, e, &temp2);
        p8est_find_edge_transform (conn, t, e, &ei);
        for (etree = 0; etree < eta->elem_count; etree++) {
          et = sc_array_index (eta, etree);
          p8est_quadrant_transform_edge (&temp2, &temp, &ei, et, true);
          nt = et->ntree;
          k = quad_is_in_array (p4est, &temp, nt, info->quads, info->treeids,
                                &level_diff, &child_id);
          SC_CHECK_ABORT (k >= 0, "Iterate: mismatched corner 15");
          if (level_diff != 0) {
            SC_CHECK_ABORT
              (child_id == *((int *) sc_array_index_int (corners, k)),
               "Iterate: mismatched edge 16");
          }
        }
      }
#endif
    }
    else {
      p4est_quadrant_corner_neighbor (q, c, &temp);
#ifndef P4_TO_P8
      p4est_find_corner_transform (conn, t, p4est_corner_to_zorder[c], cta);
#else
      p8est_find_corner_transform (conn, t, c, &ci);
#endif
      for (ctree = 0; ctree < cta->elem_count; ++ctree) {
        ct = sc_array_index (cta, ctree);
        p4est_quadrant_transform_corner (&temp, (int) ct->ncorner, true);
        nt = ct->ntree;
        k = quad_is_in_array (p4est, &temp, nt, info->quads, info->treeids,
                              &level_diff, &child_id);
        SC_CHECK_ABORT (k >= 0, "Iterate: mismatched corner 16");
        if (level_diff != 0) {
          SC_CHECK_ABORT
            (child_id == *((int *) sc_array_index_int (corners, k)),
             "Iterate: mismatched edge 17");
        }
      }
    }
  }
  sc_array_reset (cta);
#ifdef P4_TO_P8
  sc_array_reset (eta);
#endif
}

#ifdef P4_TO_P8
static void
edge_do_nothing (p8est_iter_edge_info_t * info, void *data)
{
};

static void
edge_test_adjacency (p8est_iter_edge_info_t * info, void *data)
{
  p4est_t            *p4est = info->p4est;
  p8est_quadrant_t    temp, temp2;
  p8est_quadrant_t   *q;
  int                 n = (int) info->quads->elem_count;
  int                 i, j, k, e;
  int                 child_id, level_diff;
  p4est_topidx_t      t, nt;
  sc_array_t         *common_corners = info->common_corners;
  p8est_connectivity_t *conn = p4est->connectivity;

  p8est_edge_info_t   ei;
  p8est_edge_transform_t *et;
  sc_array_t         *eta;
  size_t              etree;
  eta = &ei.edge_transforms;

  if (n == 1) {
    e = *((int *) sc_array_index (info->edges, 0));
    q = *((p8est_quadrant_t **) sc_array_index (info->quads, 0));
    SC_CHECK_ABORT (p8est_quadrant_touches_edge (q, e, true),
                    "Iterate: mismatched edge 1");
    return;
  }
  sc_array_init (eta, sizeof (p8est_edge_transform_t));
  for (i = 0; i < n; i++) {
    q = *((p8est_quadrant_t **) sc_array_index_int (info->quads, i));
    t = *((p4est_topidx_t *) sc_array_index_int (info->treeids, i));
    e = *((int *) sc_array_index_int (info->edges, i));
    for (j = 0; j < 2; j++) {
      nt = quad_face_neighbor (p4est, q, t, p8est_edge_faces[e][j], &temp);
      k = quad_is_in_array (p4est, &temp, nt, info->quads, info->treeids,
                            &level_diff, &child_id);
      SC_CHECK_ABORT (k >= 0, "Iterate: mismatched edge 2");
      SC_CHECK_ABORT (-1 <= level_diff && level_diff <= 1,
                      "Iterate: mismatched edge 3");
      if (level_diff != 0) {
        SC_CHECK_ABORT
          (child_id == *((int *) sc_array_index_int (common_corners, k)),
           "Iterate: mismatched edge 4");
      }
    }
    if (!p8est_quadrant_touches_edge (q, e, true)) {
      nt = quad_edge_neighbor (p4est, q, t, e, &temp);
      k = quad_is_in_array (p4est, &temp, nt, info->quads, info->treeids,
                            &level_diff, &child_id);
      SC_CHECK_ABORT (k >= 0, "Iterate: mismatched edge 5");
      SC_CHECK_ABORT (-1 <= level_diff && level_diff <= 1,
                      "Iterate: mismatched edge 6");
      if (level_diff != 0) {
        SC_CHECK_ABORT
          (child_id == *((int *) sc_array_index_int (common_corners, k)),
           "Iterate: mismatched edge 7");
      }
    }
    else {
      p8est_quadrant_edge_neighbor (q, e, &temp2);
      p8est_find_edge_transform (conn, t, e, &ei);
      for (etree = 0; etree < eta->elem_count; etree++) {
        et = sc_array_index (eta, etree);
        p8est_quadrant_transform_edge (&temp2, &temp, &ei, et, true);
        nt = et->ntree;
        k = quad_is_in_array (p4est, &temp, nt, info->quads, info->treeids,
                              &level_diff, &child_id);
        SC_CHECK_ABORT (k >= 0, "Iterate: mismatched edge 8");
        SC_CHECK_ABORT (-1 <= level_diff && level_diff <= 1,
                        "Iterate: mismatched edge 9");
        if (level_diff != 0) {
          SC_CHECK_ABORT
            (child_id == *((int *) sc_array_index_int (common_corners, k)),
             "Iterate: mismatched edge 10");
        }
      }
    }
  }
  sc_array_reset (eta);
}
#endif

static void
face_test_adjacency (p4est_iter_face_info_t * info, void *data)
{
  p4est_t            *p4est = info->p4est;
  p4est_quadrant_t    temp;
  p4est_quadrant_t   *left = info->left_quad;
  p4est_quadrant_t   *right = info->right_quad;
  int                 left_outgoing_face = info->left_outgoing_face;
  int                 right_outgoing_face = info->right_outgoing_face;
  int                 level_diff, child_id, j;
  p4est_topidx_t      left_treeid = info->left_treeid;
  p4est_topidx_t      right_treeid = info->right_treeid;
  p4est_topidx_t      nt;
  sc_array_t          view;
  sc_array_t          view_tree;
#ifndef P4_TO_P8
  left_outgoing_face = p4est_zface_to_rface[left_outgoing_face];
  right_outgoing_face = p4est_zface_to_rface[right_outgoing_face];
#endif

  nt =
    quad_face_neighbor (p4est, left, left_treeid, left_outgoing_face, &temp);
  view.elem_count = 1;
  view.byte_alloc = -1;
  view.elem_size = sizeof (p4est_quadrant_t *);
  view.array = (void *) (&right);
  view_tree.elem_count = 1;
  view_tree.byte_alloc = -1;
  view_tree.elem_size = sizeof (p4est_topidx_t);
  view_tree.array = (void *) (&right_treeid);
  j = quad_is_in_array (p4est, &temp, nt, &view, &view_tree, &level_diff,
                        &child_id);
  SC_CHECK_ABORT (j == 0, "Iterate: mismatched face 1");
  if (!info->is_hanging) {
    SC_CHECK_ABORT (level_diff == 0, "Iterate: mismatched face 2");
  }
  else {
    SC_CHECK_ABORT (level_diff == -1, "Iterate: mismatched face 3");
    SC_CHECK_ABORT (child_id == info->right_corner,
                    "Iterate: mismatched face 4");
  }

  nt = quad_face_neighbor (p4est, right, right_treeid,
                           right_outgoing_face, &temp);
  view.elem_count = 1;
  view.byte_alloc = -1;
  view.elem_size = sizeof (p4est_quadrant_t *);
  view.array = (void *) (&left);
  view_tree.elem_count = 1;
  view_tree.byte_alloc = -1;
  view_tree.elem_size = sizeof (p4est_topidx_t);
  view_tree.array = (void *) (&left_treeid);
  j = quad_is_in_array (p4est, &temp, nt, &view, &view_tree, &level_diff,
                        &child_id);
  SC_CHECK_ABORT (j == 0, "Iterate: mismatched face 1");
  if (!info->is_hanging) {
    SC_CHECK_ABORT (level_diff == 0, "Iterate: mismatched face 2");
  }
  else {
    SC_CHECK_ABORT (level_diff == 1, "Iterate: mismatched face 3");
    SC_CHECK_ABORT (child_id == info->left_corner,
                    "Iterate: mismatched face 4");
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
  p4est_locidx_t      num_quads, array_size, num_ghosts, ghost_size;
  sc_array_t         *QtoQ, *QtoFO, *QtoH;
  sc_array_t         *gQtoQ, *gQtoFO, *gQtoH;
  sc_array_t          ghost_layer;
  int                *ghost_owner;
  bool                success;
  int                 ntests;
  int                 i;

#ifndef P4_TO_P8
  ntests = 5;
# else
  ntests = 6;
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
      connectivity = p4est_connectivity_new_unitsquare ();
      break;
    case 1:
      connectivity = p4est_connectivity_new_corner ();
      break;
    case 2:
      connectivity = p4est_connectivity_new_moebius ();
      break;
    case 3:
      connectivity = p4est_connectivity_new_star ();
      break;
    default:
      connectivity = p4est_connectivity_new_periodic ();
      break;
#else
    case 0:
      connectivity = p8est_connectivity_new_unitcube ();
      break;
    case 1:
      connectivity = p8est_connectivity_new_periodic ();
      break;
    case 2:
      connectivity = p8est_connectivity_new_rotwrap ();
      break;
    case 3:
      connectivity = p8est_connectivity_new_twocubes ();
      break;
    case 4:
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

    /* create ghost layer: n.b. BALANCE_FULL required if iter_corner != NULL */
    sc_array_init (&ghost_layer, sizeof (p4est_quadrant_t));
    success = p4est_build_ghost_layer (p4est, P4EST_BALANCE_FULL,
                                       &ghost_layer, &ghost_owner);
    P4EST_ASSERT (success);

    /* create and fill adjacency arrays */
    num_quads = p4est->local_num_quadrants;
    array_size = num_quads * P4EST_DIM * 2;
    QtoQ = sc_array_new (sizeof (p4est_gloidx_t));
    sc_array_resize (QtoQ, (size_t) array_size);
    QtoFO = sc_array_new (sizeof (uint8_t));
    sc_array_resize (QtoFO, (size_t) array_size);
    QtoH = sc_array_new (sizeof (int8_t));
    sc_array_resize (QtoH, (size_t) array_size);
    num_ghosts = ghost_layer.elem_count;
    ghost_size = num_ghosts * P4EST_DIM * 2;
    gQtoQ = sc_array_new (sizeof (p4est_gloidx_t));
    sc_array_resize (gQtoQ, (size_t) ghost_size);
    gQtoFO = sc_array_new (sizeof (uint8_t));
    sc_array_resize (gQtoFO, (size_t) ghost_size);
    gQtoH = sc_array_new (sizeof (int8_t));
    sc_array_resize (gQtoH, (size_t) ghost_size);

    P4EST_GLOBAL_PRODUCTIONF ("Begin adjacency test %d\n", i);
#ifndef P4_TO_P8
    p4est_iterate (p4est, &ghost_layer, NULL, NULL, face_test_adjacency,
                   corner_test_adjacency);
#else
    p8est_iterate (p4est, &ghost_layer, NULL,
                   NULL, face_test_adjacency, edge_test_adjacency,
                   corner_test_adjacency);
#endif

    /* clean up */
    sc_array_destroy (QtoQ);
    sc_array_destroy (QtoFO);
    sc_array_destroy (QtoH);
    sc_array_destroy (gQtoQ);
    sc_array_destroy (gQtoFO);
    sc_array_destroy (gQtoH);
    sc_array_reset (&ghost_layer);
    P4EST_FREE (ghost_owner);

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
