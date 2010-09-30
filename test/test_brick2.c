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

#include <p4est.h>

static void
check_brick (p4est_connectivity_t * conn, int mi, int ni,
             int periodic_a, int periodic_b)
{
  p4est_topidx_t      m = (p4est_topidx_t) mi;
  p4est_topidx_t      n = (p4est_topidx_t) ni;
  int                 i;
  p4est_topidx_t      ti, tj;
  p4est_topidx_t     *tree_to_vertex = conn->tree_to_vertex;
  p4est_topidx_t     *tree_to_corner = conn->tree_to_corner;
  p4est_topidx_t     *tree_to_tree = conn->tree_to_tree;
  int8_t             *tree_to_face = conn->tree_to_face;
  p4est_topidx_t     *ctt_offset = conn->ctt_offset;
  p4est_topidx_t     *corner_to_tree = conn->corner_to_tree;
  int8_t             *corner_to_corner = conn->corner_to_corner;
  p4est_topidx_t      num_trees = conn->num_trees;
  p4est_topidx_t      num_vertices = conn->num_vertices;
  p4est_topidx_t      num_corners = conn->num_corners;
  double             *vertices = conn->vertices;
  double             *vertex[4];
  int8_t             *vert_counter, *corn_counter;
  p4est_topidx_t    **quad_counter;
  int8_t              total, face1, face2, face3, corn1;
  p4est_topidx_t      tx, ty, ttree1, ttree2, ttree3, tcorn1;
  p4est_topidx_t      diffx, diffy;

  SC_CHECK_ABORT (num_trees > 0, "no trees");
  SC_CHECK_ABORT (num_trees == m * n, "bad dimensions");
  SC_CHECK_ABORT (num_vertices == (m + 1) * (n + 1),
                  "wrong number of vertices");
  quad_counter = P4EST_ALLOC (p4est_topidx_t *, m);
  vert_counter = P4EST_ALLOC_ZERO (int8_t, num_vertices);
  corn_counter = NULL;
  if (num_corners > 0) {
    corn_counter = P4EST_ALLOC_ZERO (int8_t, num_corners);
  }

  for (ti = 0; ti < m; ti++) {
    quad_counter[ti] = P4EST_ALLOC (p4est_topidx_t, n);
    for (tj = 0; tj < n; tj++) {
      quad_counter[ti][tj] = -1;
    }
  }

  for (ti = 0; ti < num_trees; ti++) {
    for (i = 0; i < 4; i++) {
      vertex[i] = vertices + 3 * tree_to_vertex[ti * 4 + i];
      vert_counter[tree_to_vertex[ti * 4 + i]]++;
      if (num_corners > 0 && tree_to_corner[ti * 4 + i] != -1) {
        corn_counter[tree_to_corner[ti * 4 + i]]++;
      }
    }
    tx = (p4est_topidx_t) vertex[0][0];
    ty = (p4est_topidx_t) vertex[0][1];
    SC_CHECK_ABORT (tx < m, "vertex coordinates out of range");
    SC_CHECK_ABORT (ty < n, "vertex coordinates out of range");
    quad_counter[tx][ty] = ti;
    for (i = 1; i < 4; i++) {
      tx = (p4est_locidx_t) (vertex[i][0] - vertex[0][0]);
      ty = (p4est_locidx_t) (vertex[i][1] - vertex[0][1]);
      if ((i & 1) == 1) {
        SC_CHECK_ABORT (tx == 1, "non-unit vertex difference");
      }
      else {
        SC_CHECK_ABORT (tx == 0, "non-unit vertex difference");
      }
      if (((i >> 1) & 1) == 1) {
        SC_CHECK_ABORT (ty == 1, "non-unit vertex difference");
      }
      else {
        SC_CHECK_ABORT (ty == 0, "non-unit vertex difference");
      }
    }
  }

  for (ti = 0; ti < m; ti++) {
    for (tj = 0; tj < n; tj++) {
      SC_CHECK_ABORT (quad_counter[ti][tj] != -1, "grid points has no tree");
    }
  }

  for (ti = 0; ti < num_vertices; ti++) {
    tx = (p4est_topidx_t) vertices[ti * 3];
    ty = (p4est_topidx_t) vertices[ti * 3 + 1];
    total = 4;
    if (tx == m || tx == 0) {
      total /= 2;
    }
    if (ty == n || ty == 0) {
      total /= 2;
    }
    SC_CHECK_ABORT (vert_counter[ti] == total,
                    "vertex has too many or too few trees");
  }

  if (num_corners > 0) {
    for (ti = 0; ti < num_corners; ti++) {
      SC_CHECK_ABORT (corn_counter[ti] == 4,
                      "corner has too many or too few trees");
      SC_CHECK_ABORT (ctt_offset[ti] == 4 * ti, "corner offset incorrect");
    }
    SC_CHECK_ABORT (ctt_offset[ti] == 4 * ti, "corner offset incorrect");
  }

  for (ti = 0; ti < m; ti++) {
    for (tj = 0; tj < n; tj++) {
      ttree1 = quad_counter[ti][tj];
      for (face1 = 0; face1 < 4; face1++) {
        ttree2 = tree_to_tree[ttree1 * 4 + face1];
        face2 = tree_to_face[ttree1 * 4 + face1];
        if (!periodic_a &&
            ((face1 == 0 && ti == 0) || (face1 == 1 && ti == m - 1))) {
          SC_CHECK_ABORT (ttree2 == ttree1 && face2 == face1,
                          "boundary tree without boundary face");
        }
        else if (!periodic_b &&
                 ((face1 == 2 && tj == 0) || (face1 == 3 && tj == n - 1))) {
          SC_CHECK_ABORT (ttree2 == ttree1 && face2 == face1,
                          "boundary tree without boundary face");
        }
        else {
          switch (face1) {
          case 0:
            ttree3 = quad_counter[(ti + m - 1) % m][tj];
            break;
          case 1:
            ttree3 = quad_counter[(ti + 1) % m][tj];
            break;
          case 2:
            ttree3 = quad_counter[ti][(tj + n - 1) % n];
            break;
          case 3:
            ttree3 = quad_counter[ti][(tj + 1) % n];
            break;
          default:
            SC_ABORT_NOT_REACHED ();
          }
          face3 = face1 ^ 1;
          SC_CHECK_ABORT (ttree3 == ttree2 && face2 == face3,
                          "tree has incorrect neighbor");
          ttree3 = tree_to_tree[ttree2 * 4 + face2];
          SC_CHECK_ABORT (ttree1 == ttree3, "tree mismatch");
          face3 = tree_to_face[ttree2 * 4 + face2];
          SC_CHECK_ABORT (face1 == face3, "face mismatch");
        }
      }
      if (num_corners > 0) {
        for (corn1 = 0; corn1 < 4; corn1++) {
          if ((!periodic_a &&
               (((corn1 & 1) == 0 && ti == 0) ||
                ((corn1 & 1) == 1 && ti == m - 1))) ||
              (!periodic_b &&
               ((((corn1 >> 1) & 1) == 0 && tj == 0) ||
                (((corn1 >> 1) & 1) == 1 && tj == n - 1)))) {
            SC_CHECK_ABORT (tree_to_corner[ttree1 * 4 + corn1] == -1,
                            "boundary tree without boundary corner");
          }
          else {
            tcorn1 = tree_to_corner[ttree1 * 4 + corn1];
            SC_CHECK_ABORT (corner_to_tree[tcorn1 * 4 + (3 - corn1)] ==
                            ttree1, "corner_to_tree mismatch");
            SC_CHECK_ABORT (corner_to_corner[tcorn1 * 4 + 3 - corn1] ==
                            corn1, "corner_to_corner mismatch");
            for (i = 0; i < 4; i++) {
              ttree2 = corner_to_tree[tcorn1 * 4 + i];
              tx = (p4est_topidx_t) vertices[3 * tree_to_vertex[ttree2 * 4]];
              ty = (p4est_topidx_t)
                vertices[3 * tree_to_vertex[ttree2 * 4] + 1];
              diffx = (i & 1) - ((3 - corn1) & 1);
              diffy = ((i >> 1) & 1) - (((3 - corn1) >> 1) & 1);
              SC_CHECK_ABORT ((ti + diffx + m) % m == tx,
                              "unexpected trees around corner");
              SC_CHECK_ABORT ((tj + diffy + n) % n == ty,
                              "unexpected trees around corner");
            }
          }
        }
      }
    }
  }

  P4EST_FREE (vert_counter);
  if (num_corners > 0) {
    P4EST_FREE (corn_counter);
  }
  for (ti = 0; ti < m; ti++) {
    P4EST_FREE (quad_counter[ti]);
  }
  P4EST_FREE (quad_counter);

}

int
main (int argc, char **argv)
{
  int                 i, j;
  int                 l, m;
  MPI_Comm            mpicomm;
  int                 mpiret;
  int                 size, rank;
  p4est_connectivity_t *conn;

  mpiret = MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpicomm = MPI_COMM_WORLD;
  mpiret = MPI_Comm_size (mpicomm, &size);
  SC_CHECK_MPI (mpiret);
  mpiret = MPI_Comm_rank (mpicomm, &rank);
  SC_CHECK_MPI (mpiret);

  sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  for (i = 1; i <= 5; i++) {
    for (j = 1; j <= 5; j++) {
      for (l = 0; l < 2; l++) {
        for (m = 0; m < 2; m++) {
          conn = p4est_connectivity_new_brick (i, j, l, m);
          check_brick (conn, i, j, l, m);
          p4est_connectivity_destroy (conn);
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
