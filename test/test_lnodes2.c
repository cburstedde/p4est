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
#include <p4est_lnodes.h>
#else
#include <p8est_algorithms.h>
#include <p8est_bits.h>
#include <p8est_extended.h>
#include <p8est_ghost.h>
#include <p8est_lnodes.h>
#endif

#ifndef P4_TO_P8
static int          refine_level = 5;
#else
static int          refine_level = 3;
#endif

#ifndef P4_TO_P8
static p4est_connectivity_t *
p4est_connectivity_new_lnodes_test (void)
{
  const p4est_topidx_t num_vertices = 6;
  const p4est_topidx_t num_trees = 2;
  const p4est_topidx_t num_ctt = 0;
  const double        vertices[6 * 3] = {
    0, 0, 0,
    1, 0, 0,
    0, 1, 0,
    1, 1, 0,
    0, 2, 0,
    1, 2, 0,
  };
  const p4est_topidx_t tree_to_vertex[2 * 4] = {
    2, 3, 4, 5, 0, 1, 2, 3,
  };
  const p4est_topidx_t tree_to_tree[2 * 4] = {
    0, 0, 1, 0, 1, 1, 1, 0,
  };
  const int8_t        tree_to_face[2 * 4] = {
    0, 1, 3, 3, 0, 1, 2, 2,
  };

  return p4est_connectivity_new_copy (num_vertices, num_trees, 0,
                                      vertices, tree_to_vertex,
                                      tree_to_tree, tree_to_face,
                                      NULL, &num_ctt, NULL, NULL);
}
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
static int
refine_fn_lnodes_test (p4est_t * p4est, p4est_topidx_t which_tree,
                       p4est_quadrant_t * quadrant)
{
  int                 cid;

  cid = p4est_quadrant_child_id (quadrant);

  if (!which_tree && cid == 1) {
    return 1;
  }
  if (which_tree == 1 && cid == 2) {
    return 1;
  }
  return 0;
}
#endif

typedef struct tpoint
{
  p4est_topidx_t      tree;
  double              point[P4EST_DIM];
}
tpoint_t;

static void
get_point (double point[P4EST_DIM], p4est_quadrant_t * q, int i, int j,
#ifdef P4_TO_P8
           int k,
#endif
           int degree)
{
  p4est_qcoord_t      len = P4EST_QUADRANT_LEN (q->level);
  double              rlen = (double) P4EST_ROOT_LEN;
  double              deg = (double) degree;
  double              qlen = ((double) len) / rlen;

  P4EST_ASSERT (0 <= i && i < degree + 1);
  P4EST_ASSERT (0 <= j && j < degree + 1);
#ifdef P4_TO_P8
  P4EST_ASSERT (0 <= k && k < degree + 1);
#endif

  point[0] = ((double) q->x) / rlen + (((double) i) / deg) * qlen;
  point[1] = ((double) q->y) / rlen + (((double) j) / deg) * qlen;
#ifdef P4_TO_P8
  point[2] = ((double) q->z) / rlen + (((double) k) / deg) * qlen;
#endif
}

static int
same_point (tpoint_t * a, tpoint_t * b, p4est_connectivity_t * conn)
{
  double              tol = 1e-10;
  int                 a_count = 0;
  int                 b_count = 0;
  int                 a_x_pos = -1;
  int                 b_x_pos = -1;
  int                 a_y_pos = -1;
  int                 b_y_pos = -1;
#ifdef P4_TO_P8
  int                 a_z_pos = -1;
  int                 b_z_pos = -1;
#endif
  p4est_topidx_t      a_t = a->tree;
  p4est_topidx_t      b_t = b->tree;
  int                 ftrans[P4EST_FTRANSFORM];
  double              a_trans[P4EST_DIM];
  int                 a_f, b_f;
  int                 a_nf, b_nf;
  int                 a_o;
  int                 a_c, b_c;
  int                 i;
  int                 a_c_pos, b_c_pos;
  p4est_topidx_t      corner;
  p4est_topidx_t      tz;
  p4est_topidx_t      a_nt;
#ifdef P4_TO_P8
  int                 a_ref, a_set;
  int                 a_e, b_e;
  p4est_topidx_t      edge;
  int                 a_ne;
  double              a_edge_pos, b_edge_pos;
  int                 b_o;
  int                 a_e_c0, a_e_c1;
  int                 b_e_c0, b_e_c1;
#endif

  if (a_t == b_t
      && fabs (a->point[0] - b->point[0]) < tol
      && fabs (a->point[1] - b->point[1]) < tol
#ifdef P4_TO_P8
      && fabs (a->point[2] - b->point[2]) < tol
#endif
    ) {
    return 1;
  }

  if (fabs (a->point[0]) < tol) {
    a_x_pos = 0;
    a_count++;
  }
  else if (fabs (a->point[0] - 1.) < tol) {
    a_x_pos = 1;
    a_count++;
  }
  if (fabs (b->point[0]) < tol) {
    b_x_pos = 0;
    b_count++;
  }
  else if (fabs (b->point[0] - 1.) < tol) {
    b_x_pos = 1;
    b_count++;
  }

  if (fabs (a->point[1]) < tol) {
    a_y_pos = 0;
    a_count++;
  }
  else if (fabs (a->point[1] - 1.) < tol) {
    a_y_pos = 1;
    a_count++;
  }
  if (fabs (b->point[1]) < tol) {
    b_y_pos = 0;
    b_count++;
  }
  else if (fabs (b->point[1] - 1.) < tol) {
    b_y_pos = 1;
    b_count++;
  }

#ifdef P4_TO_P8
  if (fabs (a->point[2]) < tol) {
    a_z_pos = 0;
    a_count++;
  }
  else if (fabs (a->point[2] - 1.) < tol) {
    a_z_pos = 1;
    a_count++;
  }
  if (fabs (b->point[2]) < tol) {
    b_z_pos = 0;
    b_count++;
  }
  else if (fabs (b->point[2] - 1.) < tol) {
    b_z_pos = 1;
    b_count++;
  }
#endif

  if (a_count != b_count) {
    return 0;
  }

  switch (a_count) {
  case 0:
    return 0;
    break;
  case 1:
    a_f = -1;
    a_f = a_x_pos >= 0 ? a_x_pos : a_y_pos >= 0 ? a_y_pos + 2 :
#ifdef P4_TO_P8
      a_z_pos >= 0 ? a_z_pos + 4 :
#endif
      -1;
    P4EST_ASSERT (a_f >= 0 && a_f < P4EST_FACES);

    b_f = -1;
    b_f = b_x_pos >= 0 ? b_x_pos : b_y_pos >= 0 ? b_y_pos + 2 :
#ifdef P4_TO_P8
      b_z_pos >= 0 ? b_z_pos + 4 :
#endif
      -1;
    P4EST_ASSERT (b_f >= 0 && b_f < P4EST_FACES);

    if (conn->tree_to_tree[a_t * P4EST_FACES + a_f] != b_t ||
        conn->tree_to_tree[b_t * P4EST_FACES + b_f] != a_t) {
      return 0;
    }

    a_nf = (int) conn->tree_to_face[a_t * P4EST_FACES + a_f];
    a_o = a_nf / P4EST_FACES;
    a_nf %= P4EST_FACES;

    b_nf = (int) conn->tree_to_face[b_t * P4EST_FACES + b_f];
    b_nf %= P4EST_FACES;

    if (a_nf != b_f || b_nf != a_f) {
      return 0;
    }

    (void) p4est_find_face_transform (conn, a_t, a_f, ftrans);

    a_trans[ftrans[3]] = !ftrans[6] ? a->point[ftrans[0]] :
      1. - a->point[ftrans[0]];
#ifdef P4_TO_P8
    a_trans[ftrans[4]] = !ftrans[7] ? a->point[ftrans[1]] :
      1. - a->point[ftrans[1]];
#endif
    switch (ftrans[8]) {
    case 0:
      a_trans[ftrans[5]] = -a->point[ftrans[2]];
      break;
    case 1:
      a_trans[ftrans[5]] = 1. + a->point[ftrans[2]];
      break;
    case 2:
      a_trans[ftrans[5]] = a->point[ftrans[2]] - 1.;
      break;
    case 3:
      a_trans[ftrans[5]] = 2. - a->point[ftrans[2]];
      break;
    default:
      SC_ABORT_NOT_REACHED ();
    }

    if (fabs (a_trans[0] - b->point[0]) < tol
        && fabs (a_trans[1] - b->point[1]) < tol
#ifdef P4_TO_P8
        && fabs (a_trans[2] - b->point[2]) < tol
#endif
      ) {
      return 1;
    }
    else {
      return 0;
    }
  case P4EST_DIM:
    a_c = a_x_pos + 2 * a_y_pos;
    b_c = b_x_pos + 2 * b_y_pos;
#ifdef P4_TO_P8
    a_c += 4 * a_z_pos;
    b_c += 4 * b_z_pos;
#endif
    for (i = 0; i < P4EST_DIM; i++) {
      a_f = p4est_corner_faces[a_c][i];
      if (conn->tree_to_tree[a_t * P4EST_FACES + a_f] == b_t) {
        a_nf = (int) conn->tree_to_face[a_t * P4EST_FACES + a_f];
        a_o = a_nf / P4EST_FACES;
        a_nf %= P4EST_FACES;
        a_c_pos = p4est_corner_face_corners[a_c][a_f];
#ifndef P4_TO_P8
        b_c_pos = !a_o ? a_c_pos : 1 - a_c_pos;
#else
        a_ref = p8est_face_permutation_refs[a_f][a_nf];
        a_set = p8est_face_permutation_sets[a_ref][a_o];
        b_c_pos = p8est_face_permutations[a_set][a_c_pos];
#endif
        if (p4est_face_corners[a_nf][b_c_pos] == b_c) {
          return 1;
        }
      }
    }
#ifdef P4_TO_P8
    if (conn->tree_to_edge != NULL) {
      for (i = 0; i < 3; i++) {
        a_e = p8est_corner_edges[a_c][i];
        edge = conn->tree_to_edge[a_t * 12 + a_e];
        if (edge != -1) {
          a_o = -1;
          a_c_pos = (a_c >> i) & 1;
          for (tz = conn->ett_offset[edge]; tz < conn->ett_offset[edge + 1];
               tz++) {
            a_nt = conn->edge_to_tree[tz];
            a_ne = (int) conn->edge_to_edge[tz];
            if (a_nt == a_t && a_ne % 12 == a_e) {
              a_o = a_ne / 12;
            }
          }
          P4EST_ASSERT (a_o >= 0);
          for (tz = conn->ett_offset[edge]; tz < conn->ett_offset[edge + 1];
               tz++) {
            a_nt = conn->edge_to_tree[tz];
            if (a_nt != b_t) {
              continue;
            }
            b_e = (int) conn->edge_to_edge[tz];
            b_o = b_e / 12;
            b_e %= 12;
            if (b_o == a_o) {
              if (p8est_edge_corners[b_e][a_c_pos] == b_c) {
                return 1;
              }
            }
            else {
              if (p8est_edge_corners[b_e][1 - a_c_pos] == b_c) {
                return 1;
              }
            }
          }
        }
      }
    }
#endif
    if (conn->tree_to_corner == NULL) {
      return 0;
    }
    corner = conn->tree_to_corner[a_t * P4EST_CHILDREN + a_c];
    if (corner == -1) {
      return 0;
    }
    for (tz = conn->ctt_offset[corner]; tz < conn->ctt_offset[corner + 1];
         tz++) {
      a_nt = conn->corner_to_tree[tz];
      if (a_nt == b_t && (int) conn->corner_to_corner[tz] == b_c) {
        return 1;
      }
    }
    return 0;
#ifdef P4_TO_P8
  case 2:
    a_e = 0;
    if (a_z_pos >= 0) {
      a_e += a_z_pos;
    }
    if (a_y_pos >= 0) {
      a_e <<= 1;
      a_e += a_y_pos;
    }
    if (a_x_pos >= 0) {
      a_e <<= 1;
      a_e += a_x_pos;
    }
    P4EST_ASSERT (0 <= a_e && a_e < 4);
    if (a_z_pos < 0) {
      a_e += 8;
      a_edge_pos = a->point[2];
    }
    else if (a_y_pos < 0) {
      a_e += 4;
      a_edge_pos = a->point[1];
    }
    else {
      P4EST_ASSERT (a_x_pos < 0);
      a_edge_pos = a->point[0];
    }
    b_e = 0;
    if (b_z_pos >= 0) {
      b_e += b_z_pos;
    }
    if (b_y_pos >= 0) {
      b_e <<= 1;
      b_e += b_y_pos;
    }
    if (b_x_pos >= 0) {
      b_e <<= 1;
      b_e += b_x_pos;
    }
    P4EST_ASSERT (0 <= b_e && b_e < 4);
    if (b_z_pos < 0) {
      b_e += 8;
      b_edge_pos = b->point[2];
    }
    else if (b_y_pos < 0) {
      b_e += 4;
      b_edge_pos = b->point[1];
    }
    else {
      P4EST_ASSERT (b_x_pos < 0);
      b_edge_pos = b->point[0];
    }
    for (i = 0; i < 2; i++) {
      a_f = p8est_edge_faces[a_e][i];
      a_nt = conn->tree_to_tree[a_t * P4EST_FACES + a_f];
      if (a_nt != b_t) {
        continue;
      }
      a_nf = (int) conn->tree_to_face[a_t * P4EST_FACES + a_f];
      a_o = a_nf / P4EST_FACES;
      a_nf %= P4EST_FACES;
      if (p8est_edge_faces[b_e][0] != a_nf &&
          p8est_edge_faces[b_e][1] != a_nf) {
        continue;
      }
      a_e_c0 = p8est_edge_corners[a_e][0];
      a_e_c1 = p8est_edge_corners[a_e][1];

      a_ref = p8est_face_permutation_refs[a_f][a_nf];
      a_set = p8est_face_permutation_sets[a_ref][a_o];

      a_c_pos = p8est_corner_face_corners[a_e_c0][a_f];
      b_c_pos = p8est_face_permutations[a_set][a_c_pos];
      b_e_c0 = p8est_face_corners[a_nf][b_c_pos];

      if (p8est_edge_corners[b_e][0] != b_e_c0 &&
          p8est_edge_corners[b_e][1] != b_e_c0) {
        continue;
      }

      a_c_pos = p8est_corner_face_corners[a_e_c1][a_f];
      b_c_pos = p8est_face_permutations[a_set][a_c_pos];
      b_e_c1 = p8est_face_corners[a_nf][b_c_pos];

      if (p8est_edge_corners[b_e][0] != b_e_c1 &&
          p8est_edge_corners[b_e][1] != b_e_c1) {
        continue;
      }

      if (p8est_edge_corners[b_e][0] == b_e_c0) {
        if (fabs (a_edge_pos - b_edge_pos) < tol) {
          return 1;
        }
        else {
          return 0;
        }
      }
      else {
        if (fabs (a_edge_pos - (1. - b_edge_pos)) < tol) {
          return 1;
        }
        else {
          return 0;
        }
      }
    }
    if (conn->tree_to_edge == NULL) {
      return 0;
    }
    edge = conn->tree_to_edge[a_t * 12 + a_e];
    if (edge == -1) {
      return 0;
    }
    a_o = -1;
    for (tz = conn->ett_offset[edge]; tz < conn->ett_offset[edge + 1]; tz++) {
      a_nt = conn->edge_to_tree[tz];
      a_ne = (int) conn->edge_to_edge[tz];
      if (a_nt == a_t && a_ne % 12 == a_e) {
        a_o = a_ne / 12;
      }
    }
    P4EST_ASSERT (a_o != -1);
    for (tz = conn->ett_offset[edge]; tz < conn->ett_offset[edge + 1]; tz++) {
      a_nt = conn->edge_to_tree[tz];
      a_ne = (int) conn->edge_to_edge[tz];
      b_o = a_ne / 12;
      a_ne %= 12;
      if (a_ne == b_e) {
        if (a_o == b_o) {
          if (fabs (a_edge_pos - b_edge_pos) < tol) {
            return 1;
          }
          else {
            return 0;
          }
        }
        else {
          if (fabs (a_edge_pos - (1. - b_edge_pos)) < tol) {
            return 1;
          }
          else {
            return 0;
          }
        }
      }
    }
    return 0;
#endif
  default:
    SC_ABORT_NOT_REACHED ();
    return -1;
  }

}

int
main (int argc, char **argv)
{
  sc_MPI_Comm         mpicomm;
  int                 mpiret;
  int                 mpisize, mpirank;
  p4est_t            *p4est;
  p4est_connectivity_t *conn;
  p4est_ghost_t      *ghost_layer;
  p4est_ghost_t      *face_ghost_layer;
#ifdef P4_TO_P8
  p4est_ghost_t      *edge_ghost_layer;
#endif
  int                 ntests;
  int                 i, j, k;
  p4est_lnodes_t     *lnodes;
  p4est_locidx_t      nin;
  tpoint_t           *tpoints, tpoint, *tpoint_p;
  p4est_locidx_t      elid;
  p4est_locidx_t      elnid;
  p4est_locidx_t      nid;
  p4est_topidx_t      t, flt, llt;
  p4est_tree_t       *tree;
  size_t              zz, zy, count;
  p4est_lnodes_code_t fcode;
  int                 hface[P4EST_FACES];
  p4est_quadrant_t   *q, p, *q_ptr;
  int                 bcount;
  int                 iind, jind;
  int                 ib, jb;
  int                 is_hanging;
  int                 f;
  int                 c;
  sc_array_t          tpoint_array;
  p4est_lnodes_buffer_t *buffer;
  sc_array_t         *peer_buffer;
  p4est_lnodes_rank_t *lrank;
  sc_array_t         *shared_nodes;
#ifdef P4_TO_P8
  int                 hedge[12];
  int                 kind;
  int                 kb;
  int                 e;
#endif
  sc_array_t         *global_nodes;
  p4est_gloidx_t      gn;

#ifndef P4_TO_P8
  ntests = 4;
#else
  ntests = 4;
#endif

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
      conn = p4est_connectivity_new_moebius ();
      break;
    case 1:
      conn = p4est_connectivity_new_star ();
      break;
    case 2:
      conn = p4est_connectivity_new_periodic ();
      break;
    case 3:
      conn = p4est_connectivity_new_lnodes_test ();
      break;
#else
    case 0:
      conn = p8est_connectivity_new_periodic ();
      break;
    case 1:
      conn = p8est_connectivity_new_rotwrap ();
      break;
    case 2:
      conn = p8est_connectivity_new_rotcubes ();
      break;
    case 3:
      conn = p8est_connectivity_new_shell ();
      break;
#endif
    default:
      SC_ABORT_NOT_REACHED ();
      break;
    }
#ifndef P4_TO_P8
    if (i == 3) {
      p4est = p4est_new_ext (mpicomm, conn, 0, 1, 1, 0, NULL, NULL);
    }
    else {
      p4est = p4est_new_ext (mpicomm, conn, 15, 0, 0, 0, NULL, NULL);
    }
#else
    p4est = p4est_new_ext (mpicomm, conn, 15, 0, 0, 0, NULL, NULL);
#endif

    /* refine to make the number of elements interesting */
#ifndef P4_TO_P8
    if (i == 3) {
      p4est_refine (p4est, 0, refine_fn_lnodes_test, NULL);
      P4EST_ASSERT (p4est_is_balanced (p4est, P4EST_CONNECT_FULL));
    }
    else {
      p4est_refine (p4est, 1, refine_fn, NULL);
    }
#else
    p4est_refine (p4est, 1, refine_fn, NULL);
#endif

    /* balance the forest */
#ifndef P4_TO_P8
    p4est_balance (p4est, P4EST_CONNECT_FULL, NULL);
#else
    p4est_balance (p4est, P8EST_CONNECT_FULL, NULL);
#endif

    /* do a uniform partition */
#ifndef P4_TO_P8
    if (i == 3 && mpisize == 3) {
      p4est_locidx_t      num_quads[3] = { 3, 8, 3 };
      p4est_partition_given (p4est, num_quads);
    }
    else {
      p4est_partition (p4est, 0, NULL);
    }
#else
    p4est_partition (p4est, 0, NULL);
#endif

    ghost_layer = p4est_ghost_new (p4est, P4EST_CONNECT_FULL);
    face_ghost_layer = p4est_ghost_new (p4est, P4EST_CONNECT_FACE);
#ifdef P4_TO_P8
    edge_ghost_layer = p4est_ghost_new (p4est, P8EST_CONNECT_EDGE);
#endif

    flt = p4est->first_local_tree;
    llt = p4est->last_local_tree;

    for (j = -P4EST_DIM; j <= 4; j++) {
      if (!j) {
        continue;
      }
      P4EST_GLOBAL_PRODUCTIONF ("Begin lnodes test %d:%d\n", i, j);
      p4est_log_indent_push ();
      switch (j) {
#ifdef P4_TO_P8
      case -2:
        lnodes = p4est_lnodes_new (p4est, edge_ghost_layer, j);
        break;
#endif
      case -1:
        lnodes = p4est_lnodes_new (p4est, face_ghost_layer, j);
        break;
      default:
        lnodes = p4est_lnodes_new (p4est, ghost_layer, j);
        break;
      }

      if (j < 0) {
        p4est_lnodes_destroy (lnodes);
        p4est_log_indent_pop ();
        continue;
      }
      nin = lnodes->num_local_nodes;
      tpoints = P4EST_ALLOC (tpoint_t, nin);
      memset (tpoints, -1, nin * sizeof (tpoint_t));
      for (elid = 0, elnid = 0, t = flt; t <= llt; t++) {
        tree = p4est_tree_array_index (p4est->trees, t);
        count = tree->quadrants.elem_count;
        tpoint.tree = t;
        /* for every node of every element,
         * determine what kind of node it is */
        for (zz = 0; zz < count; zz++, elid++) {
          fcode = lnodes->face_code[elid];
          p4est_lnodes_decode (fcode, hface
#ifdef P4_TO_P8
                               , hedge
#endif
            );
          q = p4est_quadrant_array_index (&tree->quadrants, zz);
          p4est_quadrant_parent (q, &p);
#ifdef P4_TO_P8
          for (kind = 0; kind < j + 1; kind++) {
#endif
            for (jind = 0; jind < j + 1; jind++) {
              for (iind = 0; iind < j + 1; iind++) {
                bcount = 0;
                if (iind == 0) {
                  ib = 0;
                  bcount++;
                }
                else if (iind == j) {
                  ib = 1;
                  bcount++;
                }
                else {
                  ib = -1;
                }
                if (jind == 0) {
                  jb = 0;
                  bcount++;
                }
                else if (jind == j) {
                  jb = 1;
                  bcount++;
                }
                else {
                  jb = -1;
                }
#ifdef P4_TO_P8
                if (kind == 0) {
                  kb = 0;
                  bcount++;
                }
                else if (kind == j) {
                  kb = 1;
                  bcount++;
                }
                else {
                  kb = -1;
                }
#endif
                is_hanging = 0;
                /* if it touches a hanging part of the boundary, then the
                 * location of the element node is really the location of a
                 * node with the same index in the element's parent */
                if (fcode) {
                  switch (bcount) {
                  case 0:
                    break;
                  case 1:
                    f = (ib >= 0) ? ib : (jb >= 0) ? 2 + jb
#ifdef P4_TO_P8
                      : (kb >= 0) ? 4 + kb
#endif
                      : -1;
                    P4EST_ASSERT (f >= 0);
                    if (fcode && hface[f] >= 0) {
                      is_hanging = 1;
                    }
                    break;
                  case P4EST_DIM:
                    c = ib + 2 * jb
#ifdef P4_TO_P8
                      + 4 * kb
#endif
                      ;
                    for (k = 0; k < P4EST_DIM; k++) {
                      f = p4est_corner_faces[c][k];
                      if (fcode && hface[f] >= 0) {
                        is_hanging = 1;
                      }
                    }
#ifdef P4_TO_P8
                    for (k = 0; k < 3; k++) {
                      e = p8est_corner_edges[c][k];
                      if (fcode && hedge[e] >= 0) {
                        is_hanging = 1;
                      }
                    }
#endif
                    break;
#ifdef P4_TO_P8
                  case 2:
                    e = 0;
                    if (kb >= 0) {
                      e += kb;
                    }
                    if (jb >= 0) {
                      e <<= 1;
                      e += jb;
                    }
                    if (ib >= 0) {
                      e <<= 1;
                      e += ib;
                    }
                    P4EST_ASSERT (0 <= e && e < 4);
                    if (kb < 0) {
                      e = e + 8;
                    }
                    else if (jb < 0) {
                      e = e + 4;
                    }
#ifdef P4EST_ENABLE_DEBUG
                    else {
                      P4EST_ASSERT (ib < 0);
                    }
#endif
                    if (fcode && hedge[e] >= 0) {
                      is_hanging = 1;
                    }
                    break;
#endif
                  default:
                    SC_ABORT_NOT_REACHED ();
                  }
                }
                q_ptr = (!is_hanging) ? q : &p;
                get_point (tpoint.point, q_ptr, iind, jind,
#ifdef P4_TO_P8
                           kind,
#endif
                           j);
                nid = lnodes->element_nodes[elnid];
                if (tpoints[nid].tree == -1) {
                  tpoints[nid].tree = t;
                  tpoints[nid].point[0] = tpoint.point[0];
                  tpoints[nid].point[1] = tpoint.point[1];
#ifdef P4_TO_P8
                  tpoints[nid].point[2] = tpoint.point[2];
#endif
                }
                else {
                  SC_CHECK_ABORT (same_point (&tpoint, tpoints + nid, conn),
                                  "Lnodes: bad element-to-global node map");
                }
                elnid++;
              }
            }
#ifdef P4_TO_P8
          }
#endif
        }
      }

      sc_array_init_data (&tpoint_array, tpoints, sizeof (tpoint_t), nin);
      buffer = p4est_lnodes_share_all (&tpoint_array, lnodes);

      for (zz = 0; zz < lnodes->sharers->elem_count; zz++) {
        lrank = p4est_lnodes_rank_array_index (lnodes->sharers, zz);
        if (lrank->rank == mpirank) {
          continue;
        }
        peer_buffer =
          (sc_array_t *) sc_array_index (buffer->recv_buffers, zz);
        shared_nodes = &(lrank->shared_nodes);
        P4EST_ASSERT (shared_nodes->elem_count == peer_buffer->elem_count);
        for (zy = 0; zy < shared_nodes->elem_count; zy++) {
          nid = *((p4est_locidx_t *) sc_array_index (shared_nodes, zy));
          tpoint_p = (tpoint_t *) sc_array_index (peer_buffer, zy);
          SC_CHECK_ABORT (same_point (tpoint_p, tpoints + nid, conn),
                          "Lnodes: bad element-to-global node map across processors");
        }
      }

      p4est_lnodes_buffer_destroy (buffer);

      global_nodes = sc_array_new (sizeof (p4est_gloidx_t));
      sc_array_resize (global_nodes, lnodes->num_local_nodes);
      for (zz = 0; zz < global_nodes->elem_count; zz++) {
        *((p4est_gloidx_t *) sc_array_index (global_nodes, zz)) =
          p4est_lnodes_global_index (lnodes, zz);
      }

      p4est_lnodes_share_owned (global_nodes, lnodes);

      for (zz = 0; zz < global_nodes->elem_count; zz++) {
        gn = *((p4est_gloidx_t *) sc_array_index (global_nodes, zz));
        SC_CHECK_ABORT (gn == p4est_lnodes_global_index (lnodes, zz),
                        "Lnodes: bad global index across procesors");
      }

      sc_array_destroy (global_nodes);

      p4est_lnodes_destroy (lnodes);
      P4EST_FREE (tpoints);
      p4est_log_indent_pop ();
      P4EST_GLOBAL_PRODUCTIONF ("End lnodes test %d:%d\n", i, j);
    }

    /* clean up */
    p4est_ghost_destroy (ghost_layer);
    p4est_ghost_destroy (face_ghost_layer);
#ifdef P4_TO_P8
    p4est_ghost_destroy (edge_ghost_layer);
#endif

    p4est_destroy (p4est);
    p4est_connectivity_destroy (conn);
  }

  /* exit */
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}

/* EOF test_lnodes2.c */
