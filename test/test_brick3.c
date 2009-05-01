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

/* This file is contributed by Toby Isaac. */

#include <p8est.h>

static void
check_brick (p8est_connectivity_t * conn, p4est_topidx_t m, p4est_topidx_t n,
             p4est_topidx_t p, bool periodic_a, bool periodic_b,
             bool periodic_c)
{
  int                 i;
  p4est_topidx_t      ti, tj, tk;
  p4est_topidx_t     *tree_to_vertex = conn->tree_to_vertex;
  p4est_topidx_t     *tree_to_corner = conn->tree_to_corner;
  p4est_topidx_t     *tree_to_edge = conn->tree_to_edge;
  p4est_topidx_t     *tree_to_tree = conn->tree_to_tree;
  int8_t             *tree_to_face = conn->tree_to_face;
  p4est_topidx_t     *ett_offset = conn->ett_offset;
  p4est_topidx_t     *edge_to_tree = conn->edge_to_tree;
  int8_t             *edge_to_edge = conn->edge_to_edge;
  p4est_topidx_t     *ctt_offset = conn->ctt_offset;
  p4est_topidx_t     *corner_to_tree = conn->corner_to_tree;
  int8_t             *corner_to_corner = conn->corner_to_corner;
  p4est_topidx_t      num_trees = conn->num_trees;
  p4est_topidx_t      num_vertices = conn->num_vertices;
  p4est_topidx_t      num_corners = conn->num_corners;
  p4est_topidx_t      num_edges = conn->num_edges;
  double             *vertices = conn->vertices;
  double             *vertex[8];
  int8_t             *vert_counter, *corn_counter, *edge_counter;
  p4est_topidx_t   ***quad_counter;
  int8_t              total, face1, face2, face3, edge1, edge2, corn1;
  p4est_topidx_t      tx, ty, tz, ttree1, ttree2, ttree3, tedge1, tedge2,
    tcorn1;
  p4est_topidx_t      diffx, diffy, diffz;

  SC_CHECK_ABORT (num_trees > 0, "no trees");
  SC_CHECK_ABORT (num_trees == m * n * p, "bad dimensions");
  SC_CHECK_ABORT (num_vertices == (m + 1) * (n + 1) * (p + 1),
                  "wrong number of vertices");
  quad_counter = P4EST_ALLOC (p4est_topidx_t **, m);
  vert_counter = P4EST_ALLOC_ZERO (int8_t, num_vertices);
  if (num_corners > 0) {
    corn_counter = P4EST_ALLOC_ZERO (int8_t, num_corners);
  }
  if (num_edges > 0) {
    edge_counter = P4EST_ALLOC_ZERO (int8_t, num_edges);
  }

  for (ti = 0; ti < m; ti++) {
    quad_counter[ti] = P4EST_ALLOC (p4est_topidx_t *, n);
    for (tj = 0; tj < n; tj++) {
      quad_counter[ti][tj] = P4EST_ALLOC (p4est_topidx_t, p);
      for (tk = 0; tk < p; tk++) {
        quad_counter[ti][tj][tk] = -1;
      }
    }
  }

  for (ti = 0; ti < num_trees; ti++) {
    for (i = 0; i < 8; i++) {
      vertex[i] = vertices + 3 * tree_to_vertex[ti * 8 + i];
      vert_counter[tree_to_vertex[ti * 8 + i]]++;
      if (num_corners > 0 && tree_to_corner[ti * 8 + i] != -1) {
        corn_counter[tree_to_corner[ti * 8 + i]]++;
      }
    }
    tx = (p4est_topidx_t) vertex[0][0];
    ty = (p4est_topidx_t) vertex[0][1];
    tz = (p4est_topidx_t) vertex[0][2];
    SC_CHECK_ABORT (tx < m, "vertex coordinates out of range");
    SC_CHECK_ABORT (ty < n, "vertex coordinates out of range");
    SC_CHECK_ABORT (tz < p, "vertex coordinates out of range");
    quad_counter[tx][ty][tz] = ti;
    for (i = 1; i < 8; i++) {
      tx = (p4est_locidx_t) (vertex[i][0] - vertex[0][0]);
      ty = (p4est_locidx_t) (vertex[i][1] - vertex[0][1]);
      tz = (p4est_locidx_t) (vertex[i][2] - vertex[0][2]);
      if ((i % 2) == 1) {
        SC_CHECK_ABORT (tx == 1, "non-unit vertex difference");
      }
      else {
        SC_CHECK_ABORT (tx == 0, "non-unit vertex difference");
      }
      if ((i / 2) % 2 == 1) {
        SC_CHECK_ABORT (ty == 1, "non-unit vertex difference");
      }
      else {
        SC_CHECK_ABORT (ty == 0, "non-unit vertex difference");
      }
      if (i / 4 == 1) {
        SC_CHECK_ABORT (tz == 1, "non-unit vertex difference");
      }
      else {
        SC_CHECK_ABORT (tz == 0, "non-unit vertex difference");
      }
    }
    if (num_edges > 0) {
      for (i = 0; i < 12; i++) {
        if (tree_to_edge[ti * 12 + i] != -1) {
          edge_counter[tree_to_edge[ti * 12 + i]]++;
        }
      }
    }
  }

  for (ti = 0; ti < m; ti++) {
    for (tj = 0; tj < n; tj++) {
      for (tk = 0; tk < p; tk++) {
        SC_CHECK_ABORT (quad_counter[ti][tj][tk] != -1,
                        "grid points has no tree");
      }
    }
  }

  for (ti = 0; ti < num_vertices; ti++) {
    tx = (p4est_topidx_t) vertices[ti * 3];
    ty = (p4est_topidx_t) vertices[ti * 3 + 1];
    tz = (p4est_topidx_t) vertices[ti * 3 + 2];
    total = 8;
    if (tx == m || tx == 0) {
      total /= 2;
    }
    if (ty == n || ty == 0) {
      total /= 2;
    }
    if (tz == p || tz == 0) {
      total /= 2;
    }
    SC_CHECK_ABORT (vert_counter[ti] == total,
                    "vertex has too many or too few trees");
  }

  if (num_corners > 0) {
    for (ti = 0; ti < num_corners; ti++) {
      SC_CHECK_ABORT (corn_counter[ti] == 8,
                      "corner has too many or too few trees");
      SC_CHECK_ABORT (ctt_offset[ti] == 8 * ti, "corner offset incorrect");
    }
    SC_CHECK_ABORT (ctt_offset[ti] == 8 * ti, "corner offset incorrect");
  }

  if (num_edges > 0) {
    for (ti = 0; ti < num_edges; ti++) {
      SC_CHECK_ABORT (edge_counter[ti] == 4,
                      "edge has too many or too few trees");
      SC_CHECK_ABORT (ett_offset[ti] == 4 * ti, "edge offset incorrect");
    }
    SC_CHECK_ABORT (ett_offset[ti] == 4 * ti, "edge offset incorrect");
  }

  for (ti = 0; ti < m; ti++) {
    for (tj = 0; tj < n; tj++) {
      for (tk = 0; tk < p; tk++) {
        ttree1 = quad_counter[ti][tj][tk];
        for (face1 = 0; face1 < 6; face1++) {
          ttree2 = tree_to_tree[ttree1 * 6 + face1];
          face2 = tree_to_face[ttree1 * 6 + face1];
          if (periodic_a == false &&
              ((face1 == 0 && ti == 0) || (face1 == 1 && ti == m - 1))) {
            SC_CHECK_ABORT (ttree2 == ttree1 && face2 == face1,
                            "boundary tree without boundary face");
          }
          else if (periodic_b == false &&
                   ((face1 == 2 && tj == 0) || (face1 == 3 && tj == n - 1))) {
            SC_CHECK_ABORT (ttree2 == ttree1 && face2 == face1,
                            "boundary tree without boundary face");
          }
          else if (periodic_c == false &&
                   ((face1 == 4 && tk == 0) || (face1 == 5 && tk == p - 1))) {
            SC_CHECK_ABORT (ttree2 == ttree1 && face2 == face1,
                            "boundary tree without boundary face");
          }
          else {
            switch (face1) {
            case 0:
              ttree3 = quad_counter[(ti + m - 1) % m][tj][tk];
              break;
            case 1:
              ttree3 = quad_counter[(ti + 1) % m][tj][tk];
              break;
            case 2:
              ttree3 = quad_counter[ti][(tj + n - 1) % n][tk];
              break;
            case 3:
              ttree3 = quad_counter[ti][(tj + 1) % n][tk];
              break;
            case 4:
              ttree3 = quad_counter[ti][tj][(tk + p - 1) % p];
              break;
            default:
              ttree3 = quad_counter[ti][tj][(tk + 1) % p];
              break;
            }
            face3 = face1 ^ 1;
            SC_CHECK_ABORT (ttree3 == ttree2 && face2 == face3,
                            "tree has incorrect neighbor");
            ttree3 = tree_to_tree[ttree2 * 6 + face2];
            SC_CHECK_ABORT (ttree1 == ttree3, "tree mismatch");
            face3 = tree_to_face[ttree2 * 6 + face2];
            SC_CHECK_ABORT (face1 == face3, "face mismatch");
          }
        }
        if (num_edges > 0) {
          for (edge1 = 0; edge1 < 12; edge1++) {
            if ((periodic_b == false &&
                 (((edge1 == 0 || edge1 == 2) && (tj == 0)) ||
                  ((edge1 == 1 || edge1 == 3) && (tj == n - 1)))) ||
                (periodic_c == false &&
                 (((edge1 == 0 || edge1 == 1) && (tk == 0)) ||
                  ((edge1 == 2 || edge1 == 3) && (tk == p - 1))))) {
              SC_CHECK_ABORT (tree_to_edge[ttree1 * 12 + edge1] == -1,
                              "boundary tree without boundary edge");
            }
            else if ((periodic_a == false &&
                      (((edge1 == 4 || edge1 == 6) && (ti == 0)) ||
                       ((edge1 == 5 || edge1 == 7) && (ti == m - 1)))) ||
                     (periodic_c == false &&
                      (((edge1 == 4 || edge1 == 5) && (tk == 0)) ||
                       ((edge1 == 6 || edge1 == 7) && (tk == p - 1))))) {
              SC_CHECK_ABORT (tree_to_edge[ttree1 * 12 + edge1] == -1,
                              "boundary tree without boundary edge");
            }
            else if ((periodic_a == false &&
                      (((edge1 == 8 || edge1 == 10) && (ti == 0)) ||
                       ((edge1 == 9 || edge1 == 11) && (ti == m - 1)))) ||
                     (periodic_b == false &&
                      (((edge1 == 8 || edge1 == 9) && (tj == 0)) ||
                       ((edge1 == 10 || edge1 == 11) && (tj == n - 1))))) {
              SC_CHECK_ABORT (tree_to_edge[ttree1 * 12 + edge1] == -1,
                              "boundary tree without boundary edge");
            }
            else {
              tedge1 = tree_to_edge[ttree1 * 12 + edge1];
              SC_CHECK_ABORT (edge_to_tree[4 * tedge1 + (3 - (edge1 % 4))] ==
                              ttree1, "edge_to_tree mismatch");
              SC_CHECK_ABORT (edge_to_edge[4 * tedge1 + (3 - (edge1 % 4))] ==
                              edge1, "edge_to_edge mismatch");
              ttree2 = tree_to_tree[ttree1 * 6 + p8est_edge_faces[edge1][0]];
              edge2 = edge1 ^ 1;
              tedge2 = tree_to_edge[ttree2 * 12 + edge2];
              SC_CHECK_ABORT (tedge1 == tedge2,
                              "face neighbor trees do not share edge");
              SC_CHECK_ABORT (edge_to_tree[4 * tedge1 + (3 - (edge2 % 4))] ==
                              ttree2,
                              "edge does not recognize face neighbors");
              SC_CHECK_ABORT (edge_to_edge[4 * tedge1 + (3 - (edge2 % 4))] ==
                              edge2,
                              "edge does not recognize face neighbors' edges");
              ttree2 = tree_to_tree[ttree1 * 6 + p8est_edge_faces[edge1][1]];
              edge2 = edge1 ^ 2;
              tedge2 = tree_to_edge[ttree2 * 12 + edge2];
              SC_CHECK_ABORT (tedge1 == tedge2,
                              "face neighbor trees do not share edge");
              SC_CHECK_ABORT (edge_to_tree[4 * tedge1 + (3 - (edge2 % 4))] ==
                              ttree2,
                              "edge does not recognize face neighbors");
              SC_CHECK_ABORT (edge_to_edge[4 * tedge1 + (3 - (edge2 % 4))] ==
                              edge2,
                              "edge does not recognize face neighbors' edges");
              ttree2 =
                tree_to_tree[ttree2 * 6 + p8est_edge_faces[edge1 ^ 2][0]];
              edge2 = edge1 ^ 3;
              tedge2 = tree_to_edge[ttree2 * 12 + edge2];
              SC_CHECK_ABORT (tedge1 == tedge2,
                              "diagonal trees do not share edge");
              SC_CHECK_ABORT (edge_to_tree[4 * tedge1 + (3 - (edge2 % 4))] ==
                              ttree2,
                              "edge does not recognize diagonal trees");
              SC_CHECK_ABORT (edge_to_edge[4 * tedge1 + (3 - (edge2 % 4))] ==
                              edge2,
                              "edge does not recognize diagonal trees' edges");
            }
          }
        }
        if (num_corners > 0) {
          for (corn1 = 0; corn1 < 8; corn1++) {
            if ((periodic_a == false &&
                 ((corn1 % 2 == 0 && ti == 0) ||
                  (corn1 % 2 == 1 && ti == m - 1))) ||
                (periodic_b == false &&
                 (((corn1 / 2) % 2 == 0 && tj == 0) ||
                  ((corn1 / 2) % 2 == 1 && tj == n - 1))) ||
                (periodic_c == false &&
                 ((corn1 / 4 == 0 && tk == 0) ||
                  (corn1 / 4 == 1 && tk == p - 1)))) {
              SC_CHECK_ABORT (tree_to_corner[ttree1 * 8 + corn1] == -1,
                              "boundary tree without boundary corner");
            }
            else {
              tcorn1 = tree_to_corner[ttree1 * 8 + corn1];
              SC_CHECK_ABORT (corner_to_tree[tcorn1 * 8 + (7 - corn1)] ==
                              ttree1, "corner_to_tree mismatch");
              SC_CHECK_ABORT (corner_to_corner[tcorn1 * 8 + 7 - corn1] ==
                              corn1, "corner_to_corner mismatch");
              for (i = 0; i < 8; i++) {
                ttree2 = corner_to_tree[tcorn1 * 8 + i];
                tx =
                  (p4est_topidx_t) vertices[3 * tree_to_vertex[ttree2 * 8]];
                ty = (p4est_topidx_t)
                  vertices[3 * tree_to_vertex[ttree2 * 8] + 1];
                tz = (p4est_topidx_t)
                  vertices[3 * tree_to_vertex[ttree2 * 8] + 2];
                diffx = (i % 2) - ((7 - corn1) % 2);
                diffy = ((i / 2) % 2) - (((7 - corn1) / 2) % 2);
                diffz = (i / 4) - ((7 - corn1) / 4);
                SC_CHECK_ABORT ((ti + diffx + m) % m == tx,
                                "unexpected trees around corner");
                SC_CHECK_ABORT ((tj + diffy + n) % n == ty,
                                "unexpected trees around corner");
                SC_CHECK_ABORT ((tk + diffz + p) % p == tz,
                                "unexpected trees around corner");
              }
            }
          }
        }
      }
    }
  }

  if (num_edges > 0) {
    P4EST_FREE (edge_counter);
  }
  P4EST_FREE (vert_counter);
  if (num_corners > 0) {
    P4EST_FREE (corn_counter);
  }
  for (ti = 0; ti < m; ti++) {
    for (tj = 0; tj < n; tj++) {
      P4EST_FREE (quad_counter[ti][tj]);
    }
    P4EST_FREE (quad_counter[ti]);
  }
  P4EST_FREE (quad_counter);

}

int
main (int argc, char **argv)
{
  p4est_topidx_t      i, j, k;
  int                 l, m, n;
  MPI_Comm            mpicomm;
  int                 mpiret;
  int                 size, rank;
  p8est_connectivity_t *conn;

  mpiret = MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpicomm = MPI_COMM_WORLD;
  mpiret = MPI_Comm_size (mpicomm, &size);
  SC_CHECK_MPI (mpiret);
  mpiret = MPI_Comm_rank (mpicomm, &rank);
  SC_CHECK_MPI (mpiret);

  sc_init (mpicomm, true, true, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  for (i = 1; i <= 5; i++) {
    for (j = 1; j <= 5; j++) {
      for (k = 1; k <= 5; k++) {
        for (l = 0; l < 2; l++) {
          for (m = 0; m < 2; m++) {
            for (n = 0; n < 2; n++) {
              conn = p8est_connectivity_new_brick (i, j, k, l, m, n);
              check_brick (conn, i, j, k, l, m, n);
              p8est_connectivity_destroy (conn);
            }
          }
        }
      }
    }
  }

  /* clean up and exit */
  sc_finalize ();

  mpiret = MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}

/* EOF test_brick3.c */
