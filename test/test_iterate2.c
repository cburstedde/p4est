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
#else
static int          face_offset = 1;
static int          edge_offset = 7;
static int          corner_offset = 19;
static int          checks_per_quad = 27;
#endif

static void
volume_do_nothing (p4est_iter_volume_info_t * info, void *data)
{
};

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

static void
corner_test_adjacency (p4est_iter_corner_info_t * info, void *data)
{
  p4est_t            *p4est = info->p4est;
  int                 i, j, k, f;
  int                 n = (int) info->quads->elem_count;
  int                 c;
  int                 child_id, level_diff;
  p4est_quadrant_t   *q, temp, *ptemp;
  p4est_topidx_t      t, nt;
  sc_array_t         *corners = info->corners;
  p4est_connectivity_t *conn = p4est->connectivity;
  int                *checks = (int *) data;
  ssize_t             qid;
  sc_array_t         *trees = p4est->trees;
  p4est_tree_t       *tree;
#ifdef P4_TO_P8
  int                 e, l;
#endif
  sc_array_t          qs, ts;
  sc_array_t         *quads = &qs;
  sc_array_t         *treeids = &ts;

  sc_array_init (quads, sizeof (p4est_quadrant_t));
  sc_array_init (treeids, sizeof (p4est_topidx_t));

  if (n == 1) {
    c = *((int *) sc_array_index (info->corners, 0));
    q = *((p4est_quadrant_t **) sc_array_index (info->quads, 0));
    SC_CHECK_ABORT (p4est_quadrant_touches_corner (q, c, true),
                    "Iterate: mismatched corner 1");

    qid = *((ssize_t *) sc_array_index (info->quadids, 0));
    t = *((p4est_topidx_t *) sc_array_index_int (info->treeids, 0));
    if (qid >= 0) {
      tree = p4est_array_index_topidx (trees, t);
      qid += tree->quadrants_offset;
      checks[qid * checks_per_quad + corner_offset + c]++;
    }
    return;
  }

  for (i = 0; i < n; i++) {
    t = *((p4est_topidx_t *) sc_array_index_int (info->treeids, i));
    q = *((p4est_quadrant_t **) sc_array_index_int (info->quads, i));
    c = *((int *) sc_array_index_int (corners, i));
    for (j = 0; j < P4EST_DIM; j++) {
      f = p4est_corner_faces[c][j];
      nt = p4est_quadrant_face_neighbor_extra (q, t, f, &temp, conn);
      if (nt == -1) {
        continue;
      }
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
      p8est_quadrant_edge_neighbor_extra (q, t, e, quads, treeids, conn);
      for (l = 0; (size_t) l < quads->elem_count; l++) {
        ptemp = sc_array_index_int (quads, l);
        nt = *((p4est_topidx_t *) sc_array_index_int (treeids, l));
        k = quad_is_in_array (p4est, ptemp, nt, info->quads, info->treeids,
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
      sc_array_reset (quads);
      sc_array_reset (treeids);
#endif
    }

    p4est_quadrant_corner_neighbor_extra (q, t, c, quads, treeids, conn);
    for (j = 0; (size_t) j < quads->elem_count; j++) {
      ptemp = sc_array_index_int (quads, j);
      nt = *((p4est_topidx_t *) sc_array_index_int (treeids, j));
      k = quad_is_in_array (p4est, ptemp, nt, info->quads, info->treeids,
                            &level_diff, &child_id);
      SC_CHECK_ABORT (k >= 0, "Iterate: mismatched corner 8");
      if (level_diff != 0) {
        SC_CHECK_ABORT
          (child_id == *((int *) sc_array_index_int (corners, k)),
           "Iterate: mismatched corner 9");
      }
    }
    sc_array_reset (quads);
    sc_array_reset (treeids);

    qid = *((ssize_t *) sc_array_index_int (info->quadids, i));
    if (qid >= 0) {
      tree = p4est_array_index_topidx (trees, t);
      qid += tree->quadrants_offset;
      checks[qid * checks_per_quad + corner_offset + c]++;
    }
  }
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
  p8est_quadrant_t    temp, *ptemp;
  p8est_quadrant_t   *q;
  int                 n = (int) info->quads->elem_count;
  int                 i, j, k, e, c, c2, f;
  int                 child_id, level_diff;
  p4est_topidx_t      t, nt;
  sc_array_t         *common_corners = info->common_corners;
  p8est_connectivity_t *conn = p4est->connectivity;

  int                *checks = (int *) data;
  ssize_t             qid;
  bool                hanging;
  p4est_tree_t       *tree;
  sc_array_t         *trees = p4est->trees;
  sc_array_t          qs, ts;
  sc_array_t         *quads = &qs;
  sc_array_t         *treeids = &ts;

  sc_array_init (quads, sizeof (p4est_quadrant_t));
  sc_array_init (treeids, sizeof (p4est_topidx_t));

  if (n == 1) {
    e = *((int *) sc_array_index (info->edges, 0));
    q = *((p8est_quadrant_t **) sc_array_index (info->quads, 0));
    SC_CHECK_ABORT (p8est_quadrant_touches_edge (q, e, true),
                    "Iterate: mismatched edge 1");

    qid = *((ssize_t *) sc_array_index (info->quadids, 0));
    t = *((p4est_topidx_t *) sc_array_index_int (info->treeids, 0));
    if (qid >= 0) {
      tree = p4est_array_index_topidx (trees, t);
      qid += tree->quadrants_offset;
      checks[qid * checks_per_quad + edge_offset + e]++;
    }
    return;
  }
  for (i = 0; i < n; i++) {
    hanging = false;
    q = *((p8est_quadrant_t **) sc_array_index_int (info->quads, i));
    t = *((p4est_topidx_t *) sc_array_index_int (info->treeids, i));
    e = *((int *) sc_array_index_int (info->edges, i));
    for (j = 0; j < 2; j++) {
      f = p8est_edge_faces[e][j];
      nt = p8est_quadrant_face_neighbor_extra (q, t, f, &temp, conn);
      if (nt == -1) {
        continue;
      }
      k = quad_is_in_array (p4est, &temp, nt, info->quads, info->treeids,
                            &level_diff, &child_id);
      SC_CHECK_ABORT (k >= 0, "Iterate: mismatched edge 2");
      SC_CHECK_ABORT (-1 <= level_diff && level_diff <= 1,
                      "Iterate: mismatched edge 3");
      if (level_diff == 1) {
        hanging = true;
      }
      if (level_diff != 0) {
        SC_CHECK_ABORT
          (child_id == *((int *) sc_array_index_int (common_corners, k)),
           "Iterate: mismatched edge 4");
      }
    }
    p8est_quadrant_edge_neighbor_extra (q, t, e, quads, treeids, conn);
    for (j = 0; (size_t) j < quads->elem_count; j++) {
      ptemp = sc_array_index_int (quads, j);
      nt = *((p4est_topidx_t *) sc_array_index_int (treeids, j));
      k = quad_is_in_array (p4est, ptemp, nt, info->quads, info->treeids,
                            &level_diff, &child_id);
      SC_CHECK_ABORT (k >= 0, "Iterate: mismatched edge 5");
      SC_CHECK_ABORT (-1 <= level_diff && level_diff <= 1,
                      "Iterate: mismatched edge 6");
      if (level_diff == 1) {
        hanging = true;
      }
      if (level_diff != 0) {
        SC_CHECK_ABORT
          (child_id == *((int *) sc_array_index_int (common_corners, k)),
           "Iterate: mismatched edge 7");
      }
    }
    sc_array_reset (quads);
    sc_array_reset (treeids);

    qid = *((ssize_t *) sc_array_index_int (info->quadids, i));
    if (qid >= 0) {
      tree = p4est_array_index_topidx (trees, t);
      qid += tree->quadrants_offset;
      checks[qid * checks_per_quad + edge_offset + e]++;
      if (hanging) {
        c = *((int *) sc_array_index_int (common_corners, i));
        if (p8est_edge_corners[e][0] == c) {
          c2 = p8est_edge_corners[e][1];
        }
        else {
          c2 = p8est_edge_corners[e][0];
        }
        checks[qid * checks_per_quad + corner_offset + c2]++;
      }
      else if (checks[qid * checks_per_quad + edge_offset + e] == 2) {
        checks[qid * checks_per_quad + edge_offset + e] = 1;
      }
    }
  }
}
#endif

static void
face_test_adjacency (p4est_iter_face_info_t * info, void *data)
{
  p4est_t            *p4est = info->p4est;
  p4est_quadrant_t    temp;
  p4est_quadrant_t   *left = info->left_quad;
  p4est_quadrant_t   *right = info->right_quad;
  int                 left_face = info->left_outgoing_face;
  int                 right_face = info->right_outgoing_face;
  int                 level_diff, child_id, j;
  p4est_topidx_t      left_treeid = info->left_treeid;
  p4est_topidx_t      right_treeid = info->right_treeid;
  p4est_topidx_t      nt;
  sc_array_t          view;
  sc_array_t          view_tree;
  int                *checks = (int *) data;
  sc_array_t         *trees = p4est->trees;
  p4est_tree_t       *tree;
  ssize_t             qid;
  int                 c;
#ifdef P4_TO_P8
  int                 e, dir;
#endif
  p4est_connectivity_t *conn = p4est->connectivity;

#ifndef P4_TO_P8
  nt = p4est_quadrant_face_neighbor_extra (left, left_treeid,
                                           p4est_zface_to_rface[left_face],
                                           &temp, conn);
#else
  nt =
    p4est_quadrant_face_neighbor_extra (left, left_treeid, left_face, &temp,
                                        conn);
#endif
  if (nt != -1) {
    sc_array_init_data (&view, &right, sizeof (p4est_quadrant_t *), 1);
    sc_array_init_data (&view_tree, &right_treeid, sizeof (p4est_topidx_t),
                        1);
    j =
      quad_is_in_array (p4est, &temp, nt, &view, &view_tree, &level_diff,
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

#ifndef P4_TO_P8
    nt = p4est_quadrant_face_neighbor_extra (right, right_treeid,
                                             p4est_zface_to_rface[right_face],
                                             &temp, conn);
#else
    nt = p4est_quadrant_face_neighbor_extra (right, right_treeid, right_face,
                                             &temp, conn);
#endif
    sc_array_init_data (&view, &left, sizeof (p4est_quadrant_t *), 1);
    sc_array_init_data (&view_tree, &left_treeid, sizeof (p4est_topidx_t), 1);
    j = quad_is_in_array (p4est, &temp, nt, &view, &view_tree, &level_diff,
                          &child_id);
    SC_CHECK_ABORT (j == 0, "Iterate: mismatched face 5");
    if (!info->is_hanging) {
      SC_CHECK_ABORT (level_diff == 0, "Iterate: mismatched face 6");
    }
    else {
      SC_CHECK_ABORT (level_diff == 1, "Iterate: mismatched face 7");
      SC_CHECK_ABORT (child_id == info->left_corner,
                      "Iterate: mismatched face 8");
    }
  }

  qid = info->left_quadid;
  if (qid >= 0) {
    tree = p4est_array_index_topidx (trees, left_treeid);
    qid += tree->quadrants_offset;
    checks[qid * checks_per_quad + face_offset + left_face]++;
    if (checks[qid * checks_per_quad + face_offset + left_face] ==
        P4EST_CHILDREN / 2) {
      checks[qid * checks_per_quad + face_offset + left_face] = 1;
    }
  }
  if (left_face == right_face && left_treeid == right_treeid) {
    return;
  }
  qid = info->right_quadid;
  if (qid >= 0) {
    tree = p4est_array_index_topidx (trees, right_treeid);
    qid += tree->quadrants_offset;
    checks[qid * checks_per_quad + face_offset + right_face]++;
    if (!info->is_hanging) {
      return;
    }
#ifndef P4_TO_P8
    if (p4est_face_corners[p4est_zface_to_rface[right_face]][0]
        == info->right_corner) {
      c = p4est_face_corners[p4est_zface_to_rface[right_face]][1];
    }
    else {
      c = p4est_face_corners[p4est_zface_to_rface[right_face]][0];
    }
    checks[qid * checks_per_quad + corner_offset + c]++;
#else
    c = p8est_corner_face_corners[info->right_corner][right_face];
    P4EST_ASSERT (c >= 0);
    c = p8est_face_corners[right_face][3 - c];
    checks[qid * checks_per_quad + corner_offset + c]++;
    dir = right_face / 2;
    e = p8est_corner_edges[c][(dir + 1) % 3];
    checks[qid * checks_per_quad + edge_offset + e]++;
    e = p8est_corner_edges[c][(dir + 2) % 3];
    checks[qid * checks_per_quad + edge_offset + e]++;
#endif
  }
}

static void
volume_test_adjacency (p4est_iter_volume_info_t * info, void *data)
{
  int                *checks = (int *) data;
  size_t              qid = info->quadid;
  p4est_topidx_t      t = info->treeid;
  p4est_tree_t       *tree = p4est_array_index_topidx (info->p4est->trees, t);

  qid += tree->quadrants_offset;
  checks[qid * checks_per_quad]++;
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

    num_quads = p4est->local_num_quadrants;
    num_checks = checks_per_quad * num_quads;
    checks = P4EST_ALLOC_ZERO (int, num_checks);

    P4EST_GLOBAL_PRODUCTIONF ("Begin adjacency test %d\n", i);
    p4est_iterate (p4est, &ghost_layer, checks, volume_test_adjacency,
                   face_test_adjacency,
#ifdef P4_TO_P8
                   edge_test_adjacency,
#endif
                   corner_test_adjacency);

    for (li = 0; li < num_checks; li++) {
      if (checks[li] != 1) {
        P4EST_VERBOSEF ("qid = %d\n", li / checks_per_quad);
        P4EST_VERBOSEF ("check = %d\n", li % checks_per_quad);
        P4EST_VERBOSEF ("count = %d\n", checks[li]);
      }
      SC_CHECK_ABORT (checks[li] == 1, "Iterate: completion check");
    }

    /* clean up */
    P4EST_FREE (checks);
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
