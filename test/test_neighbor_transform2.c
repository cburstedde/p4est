/*
  This file is part of p4est.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2010 The University of Texas System
  Additional copyright (C) 2011 individual authors
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
#include <p4est_connectivity.h>
#include <p4est_bits.h>
#else
#include <p8est_connectivity.h>
#include <p8est_bits.h>
#endif

static void
quad_coords (const p4est_quadrant_t * quad, p4est_qcoord_t coords[P4EST_DIM],
             int corner)
{
  p4est_qcoord_t      h = P4EST_QUADRANT_LEN (quad->level);

  coords[0] = quad->x;
  coords[1] = quad->y;
#ifdef P4_TO_P8
  coords[2] = quad->z;
#endif
  for (int d = 0; d < P4EST_DIM; d++) {
    coords[d] += (corner & (1 << d)) ? h : 0;
  }
}

#ifdef P4EST_ENABLE_DEBUG

static int
coords_equal (const p4est_qcoord_t A[], p4est_qcoord_t B[])
{
  for (int d = 0; d < P4EST_DIM; d++) {
    if (A[d] != B[d]) {
      return 0;
    }
  }
  return 1;
}

#endif

static void
test_find_self_transform (sc_array_t * neigh_transforms, int t, int found[])
{
  for (size_t iz = 0; iz < neigh_transforms->elem_count; iz++) {
    p4est_neighbor_transform_t *nt =
      (p4est_neighbor_transform_t *) sc_array_index (neigh_transforms, iz);

    if (found[iz])
      continue;
    if (nt->neighbor_type == P4EST_CONNECT_SELF) {
      P4EST_ASSERT (nt->neighbor == t);
      P4EST_ASSERT (nt->index_self == 0);
      P4EST_ASSERT (nt->index_neighbor == 0);
      for (int d = 0; d < P4EST_DIM; d++) {
        P4EST_ASSERT (nt->origin_self[d] == nt->origin_neighbor[d]);
        P4EST_ASSERT (nt->perm[d] == d);
        P4EST_ASSERT (nt->sign[d] == 1);
      }

      for (uint64_t id = 0; id < (1 << (P4EST_DIM * 2)); id++) {
        p4est_quadrant_t    self_q, neigh_q;
        p4est_quadrant_set_morton (&self_q, 2, id);

        /* test forward transform */
        P4EST_QUADRANT_INIT (&neigh_q);
        p4est_neighbor_transform_quadrant (nt, &self_q, &neigh_q);
        P4EST_ASSERT (p4est_quadrant_compare (&self_q, &neigh_q) == 0);

        /* test reverse transform */
        P4EST_QUADRANT_INIT (&self_q);
        p4est_neighbor_transform_quadrant_reverse (nt, &neigh_q, &self_q);
        P4EST_ASSERT (p4est_quadrant_compare (&self_q, &neigh_q) == 0);

        for (int c = 0; c < P4EST_CHILDREN; c++) {
          p4est_qcoord_t      self_coords[P4EST_DIM];
          p4est_qcoord_t      neigh_coords[P4EST_DIM];

          quad_coords (&self_q, self_coords, c);

          for (int d = 0; d < P4EST_DIM; d++) {
            neigh_coords[d] = -1;
          }
          p4est_neighbor_transform_coordinates (nt, self_coords,
                                                neigh_coords);
          P4EST_ASSERT (coords_equal (self_coords, neigh_coords));

          for (int d = 0; d < P4EST_DIM; d++) {
            self_coords[d] = -1;
          }
          p4est_neighbor_transform_coordinates_reverse (nt, neigh_coords,
                                                        self_coords);
          P4EST_ASSERT (coords_equal (self_coords, neigh_coords));
        }
      }
      P4EST_GLOBAL_LDEBUG ("Self transformation ok\n");
      found[iz] = 1;
      return;
    }
  }
}

static void
test_find_face_transform (sc_array_t * neigh_transforms, int j, int ntree,
                          const int ftransform[9], int found[])
{
  for (size_t iz = 0; iz < neigh_transforms->elem_count; iz++) {
    p4est_neighbor_transform_t *nt =
      (p4est_neighbor_transform_t *) sc_array_index (neigh_transforms, iz);

    if (found[iz])
      continue;
    if (nt->neighbor_type == P4EST_CONNECT_FACE && nt->index_self == j) {
      P4EST_ASSERT (nt->neighbor == ntree);
      P4EST_ASSERT (0 <= nt->index_neighbor
                    && nt->index_neighbor < P4EST_FACES);

      for (uint64_t id = 0; id < (1 << (P4EST_DIM * 2)); id++) {
        p4est_quadrant_t    self_quad[2];
        int                 fs[2];

        fs[0] = j;
        fs[1] = j ^ 1;

        /* only test quadrants that touch face j */
        p4est_quadrant_set_morton (&self_quad[0], 2, id);
        p4est_quadrant_face_neighbor (&self_quad[0], j, &self_quad[1]);
        if (p4est_quadrant_is_inside_root (&self_quad[1])) {
          continue;
        }

        for (int s = 0; s < 2; s++) {
          p4est_quadrant_t    neigh_quad, neigh_quad_f, self_quad_rev;

          /* compare neighbor transform and face transform */
          p4est_quadrant_transform_face (&self_quad[s], &neigh_quad_f,
                                         ftransform);
          p4est_neighbor_transform_quadrant (nt, &self_quad[s], &neigh_quad);
          P4EST_ASSERT (p4est_quadrant_compare (&neigh_quad, &neigh_quad_f) ==
                        0);

          /* reverse is a left inverse */
          p4est_neighbor_transform_quadrant_reverse (nt, &neigh_quad,
                                                     &self_quad_rev);
          P4EST_ASSERT (p4est_quadrant_compare (&self_quad[s], &self_quad_rev)
                        == 0);

          for (int fc = 0; fc < P4EST_CHILDREN / 2; fc++) {
            p4est_qcoord_t      self_coords[P4EST_DIM];
            p4est_qcoord_t      neigh_coords[P4EST_DIM];
            p4est_qcoord_t      neigh_coords_f[P4EST_DIM];
            p4est_qcoord_t      self_coords_rev[P4EST_DIM];
            int                 c = p4est_face_corners[fs[s]][fc];

            quad_coords (&self_quad[s], self_coords, c);
            p4est_coordinates_transform_face (self_coords, neigh_coords_f,
                                              ftransform);
            p4est_neighbor_transform_coordinates (nt, self_coords,
                                                  neigh_coords);
            P4EST_ASSERT (coords_equal (neigh_coords_f, neigh_coords));

            p4est_neighbor_transform_coordinates_reverse (nt, neigh_coords,
                                                          self_coords_rev);
            P4EST_ASSERT (coords_equal (self_coords, self_coords_rev));
          }
        }

      }
      P4EST_GLOBAL_LDEBUGF ("Face transformation across %d ok\n", j);
      found[iz] = 1;
      return;
    }
  }
}

static void
test_self_transform (p4est_connectivity_t * conn, p4est_topidx_t t, int j,
                     sc_array_t * neigh_transforms)
{
  int                 found = 0;

  P4EST_ASSERT (neigh_transforms->elem_count == 1);
  test_find_self_transform (neigh_transforms, t, &found);
  P4EST_ASSERT (found);
}

static void
test_face_transform (p4est_connectivity_t * conn, p4est_topidx_t t, int j,
                     sc_array_t * neigh_transforms)
{
  int                *found =
    P4EST_ALLOC_ZERO (int, neigh_transforms->elem_count);
  p4est_topidx_t      ntree;
  int                 ftransform[9];

  P4EST_ASSERT (neigh_transforms->elem_count == 1
                || neigh_transforms->elem_count == 2);
  test_find_self_transform (neigh_transforms, t, found);
  ntree = p4est_find_face_transform (conn, t, j, ftransform);
  if (ntree >= 0) {
    test_find_face_transform (neigh_transforms, j, ntree, ftransform, found);
  }
  for (size_t iz = 0; iz < neigh_transforms->elem_count; iz++) {
    P4EST_ASSERT (found[iz]);
  }
  P4EST_FREE (found);
}

#ifdef P4_TO_P8
static void
test_find_edge_transform (sc_array_t * neigh_transforms, p4est_topidx_t t,
                          int j, p8est_edge_info_t * ei,
                          p8est_edge_transform_t * et, int found[])
{
  for (size_t iz = 0; iz < neigh_transforms->elem_count; iz++) {
    p4est_neighbor_transform_t *nt =
      (p4est_neighbor_transform_t *) sc_array_index (neigh_transforms, iz);

    if (found[iz])
      continue;
    if (nt->neighbor_type == P8EST_CONNECT_EDGE && nt->index_self == j) {
      P4EST_ASSERT (nt->neighbor == et->ntree);
      P4EST_ASSERT (nt->index_neighbor == et->nedge);

      for (uint64_t id = 0; id < (1 << (P4EST_DIM * 2)); id++) {
        p4est_quadrant_t    self_quad[2];
        int                 es[2];

        es[0] = j;
        es[1] = j ^ 3;

        /* only test quadrants that touch face j */
        p4est_quadrant_set_morton (&self_quad[0], 2, id);
        p8est_quadrant_edge_neighbor (&self_quad[0], j, &self_quad[1]);
        if (!p8est_quadrant_is_outside_edge (&self_quad[1])) {
          continue;
        }

        for (int s = 0; s < 2; s++) {
          p4est_quadrant_t    neigh_quad, neigh_quad_e, self_quad_rev;

          /* compare neighbor transform and face transform */
          p8est_quadrant_transform_edge (&self_quad[s], &neigh_quad_e, ei, et,
                                         s);
          p4est_neighbor_transform_quadrant (nt, &self_quad[s], &neigh_quad);
          P4EST_ASSERT (p4est_quadrant_compare (&neigh_quad, &neigh_quad_e) ==
                        0);

          /* reverse is a left inverse */
          p4est_neighbor_transform_quadrant_reverse (nt, &neigh_quad,
                                                     &self_quad_rev);
          P4EST_ASSERT (p4est_quadrant_compare (&self_quad[s], &self_quad_rev)
                        == 0);

          for (int ec = 0; ec < 2; ec++) {
            p4est_qcoord_t      self_coords[P4EST_DIM];
            p4est_qcoord_t      neigh_coords[P4EST_DIM];
            p4est_qcoord_t      neigh_coords_e[P4EST_DIM];
            p4est_qcoord_t      self_coords_rev[P4EST_DIM];
            int                 c = p8est_edge_corners[es[s]][ec];

            quad_coords (&self_quad[s], self_coords, c);
            p8est_coordinates_transform_edge (self_coords, neigh_coords_e, ei,
                                              et);
            p4est_neighbor_transform_coordinates (nt, self_coords,
                                                  neigh_coords);
            P4EST_ASSERT (coords_equal (neigh_coords_e, neigh_coords));

            p4est_neighbor_transform_coordinates_reverse (nt, neigh_coords,
                                                          self_coords_rev);
            P4EST_ASSERT (coords_equal (self_coords, self_coords_rev));
          }
        }

      }
      P4EST_GLOBAL_LDEBUGF ("Edge transformation across %d to %d ok\n", j,
                            et->ntree);
      found[iz] = 1;
      return;
    }
  }
}

static void
test_find_all_edge_transforms (p4est_connectivity_t * conn, p4est_topidx_t t,
                               int j, sc_array_t * neigh_transforms,
                               int found[])
{
  p8est_edge_info_t   ei;
  sc_array_t         *eta = &ei.edge_transforms;

  sc_array_init (eta, sizeof (p8est_edge_transform_t));
  p8est_find_edge_transform (conn, t, j, &ei);
  for (size_t iz = 0; iz < eta->elem_count; iz++) {
    p8est_edge_transform_t *et =
      (p8est_edge_transform_t *) sc_array_index (eta, iz);

    test_find_edge_transform (neigh_transforms, t, j, &ei, et, found);
  }
  sc_array_reset (eta);
}

static void
test_edge_transform (p4est_connectivity_t * conn, p4est_topidx_t t, int j,
                     sc_array_t * neigh_transforms)
{
  int                *found =
    P4EST_ALLOC_ZERO (int, neigh_transforms->elem_count);

  P4EST_ASSERT (neigh_transforms->elem_count >= 1);
  test_find_self_transform (neigh_transforms, t, found);
  for (int ef = 0; ef < 2; ef++) {
    int                 f = p8est_edge_faces[j][ef];
    int                 ftransform[9];
    p4est_topidx_t      ntree;

    ntree = p4est_find_face_transform (conn, t, f, ftransform);
    if (ntree >= 0) {
      test_find_face_transform (neigh_transforms, f, ntree, ftransform,
                                found);
    }
  }
  test_find_all_edge_transforms (conn, t, j, neigh_transforms, found);
  for (size_t iz = 0; iz < neigh_transforms->elem_count; iz++) {
    P4EST_ASSERT (found[iz]);
  }
  P4EST_FREE (found);
}
#endif

static void
test_find_corner_transform (sc_array_t * neigh_transforms, p4est_topidx_t t,
                            int j, p4est_corner_transform_t * ct, int found[])
{
  for (size_t iz = 0; iz < neigh_transforms->elem_count; iz++) {
    p4est_neighbor_transform_t *nt =
      (p4est_neighbor_transform_t *) sc_array_index (neigh_transforms, iz);

    if (found[iz])
      continue;
    if (nt->neighbor_type == P4EST_CONNECT_CORNER && nt->index_self == j) {
      p4est_quadrant_t    root, self_quad[2];
      int                 cs[2], ncs[2];

      cs[0] = j;
      cs[1] = j ^ (P4EST_CHILDREN - 1);
      ncs[0] = (ct->ncorner) ^ (P4EST_CHILDREN - 1);
      ncs[1] = ct->ncorner;
      p4est_qcoord_t      self_coords[P4EST_DIM];
      p4est_qcoord_t      neigh_coords_c[P4EST_DIM];
      p4est_qcoord_t      neigh_coords[P4EST_DIM];
      p4est_qcoord_t      self_coords_rev[P4EST_DIM];

      P4EST_ASSERT (nt->neighbor == ct->ntree);
      P4EST_ASSERT (nt->index_neighbor == ct->ncorner);

      memset (&root, 0, sizeof (root));
      p4est_quadrant_corner_descendant (&root, &self_quad[0], j, 2);
      p4est_quadrant_corner_neighbor (&self_quad[0], j, &self_quad[1]);
      for (int s = 0; s < 2; s++) {
        p4est_quadrant_t    neigh_quad_c, neigh_quad, self_quad_rev;

        neigh_quad_c = self_quad[s];
        p4est_quadrant_transform_corner (&neigh_quad_c, ct->ncorner, s);
        p4est_neighbor_transform_quadrant (nt, &self_quad[s], &neigh_quad);
        P4EST_ASSERT (p4est_quadrant_compare (&neigh_quad_c, &neigh_quad) ==
                      0);

        p4est_neighbor_transform_quadrant_reverse (nt, &neigh_quad,
                                                   &self_quad_rev);
        P4EST_ASSERT (p4est_quadrant_compare (&self_quad[s], &self_quad_rev)
                      == 0);

        quad_coords (&self_quad[s], self_coords, cs[s]);
        quad_coords (&neigh_quad, neigh_coords_c, ncs[s]);
        p4est_neighbor_transform_coordinates (nt, self_coords, neigh_coords);
        P4EST_ASSERT (coords_equal (neigh_coords_c, neigh_coords));

        p4est_neighbor_transform_coordinates_reverse (nt, neigh_coords,
                                                      self_coords_rev);
        P4EST_ASSERT (coords_equal (self_coords, self_coords_rev));
      }
      P4EST_GLOBAL_LDEBUGF ("Corner transformation across %d to %d ok\n", j,
                            ct->ntree);
      found[iz] = 1;
      return;
    }
  }
}

static void
test_corner_transform (p4est_connectivity_t * conn, p4est_topidx_t t, int j,
                       sc_array_t * neigh_transforms)
{
  int                *found =
    P4EST_ALLOC_ZERO (int, neigh_transforms->elem_count);
  p4est_corner_info_t ci;
  sc_array_t         *cta = &ci.corner_transforms;;

  P4EST_ASSERT (neigh_transforms->elem_count >= 1);
  test_find_self_transform (neigh_transforms, t, found);
  for (int cf = 0; cf < P4EST_DIM; cf++) {
    int                 f = p4est_corner_faces[j][cf];
    int                 ftransform[9];
    p4est_topidx_t      ntree;

    ntree = p4est_find_face_transform (conn, t, f, ftransform);
    if (ntree >= 0) {
      test_find_face_transform (neigh_transforms, f, ntree, ftransform,
                                found);
    }
  }
#ifdef P4_TO_P8
  for (int ce = 0; ce < P4EST_DIM; ce++) {
    int                 e = p8est_corner_edges[j][ce];

    test_find_all_edge_transforms (conn, t, e, neigh_transforms, found);
  }
#endif
  sc_array_init (cta, sizeof (p4est_corner_transform_t));
  p4est_find_corner_transform (conn, t, j, &ci);
  for (size_t iz = 0; iz < cta->elem_count; iz++) {
    p4est_corner_transform_t *ct =
      (p4est_corner_transform_t *) sc_array_index (cta, iz);

    test_find_corner_transform (neigh_transforms, t, j, ct, found);
  }
  sc_array_reset (cta);
  for (size_t iz = 0; iz < neigh_transforms->elem_count; iz++) {
    P4EST_ASSERT (found[iz]);
  }
  P4EST_FREE (found);
}

int
main (int argc, char **argv)
{
  sc_MPI_Init (&argc, &argv);
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_DEBUG);
  p4est_init (NULL, SC_LP_DEBUG);
#ifndef P4_TO_P8
#define N_NAMES 13
  const char         *names[N_NAMES] =
    { "brick23", "corner", "cubed", "disk", "icosahedron", "moebius",
    "periodic", "pillow", "rotwrap", "star", "shell2d", "disk2d", "unit"
  };
  const p4est_connect_type_t neigh_types[3] =
    { P4EST_CONNECT_SELF, P4EST_CONNECT_FACE, P4EST_CONNECT_CORNER };
  const char         *neigh_type_names[3] = { "self", "face", "corner" };
  const int           neigh_sizes[3] = { 1, P4EST_FACES, P4EST_CHILDREN };
#else
#define N_NAMES 9
  const char         *names[] =
    { "brick235", "periodic", "rotcubes", "rotwrap", "shell", "sphere",
    "twocubes", "twowrap", "unit"
  };
  const p4est_connect_type_t neigh_types[4] =
    { P4EST_CONNECT_SELF, P4EST_CONNECT_FACE, P8EST_CONNECT_EDGE,
    P4EST_CONNECT_CORNER
  };
  const char         *neigh_type_names[4] =
    { "self", "face", "edge", "corner" };
  const int           neigh_sizes[4] =
    { 1, P4EST_FACES, P8EST_EDGES, P4EST_CHILDREN };
#endif

  for (int i = 0; i < N_NAMES; i++) {
    p4est_connectivity_t *conn = p4est_connectivity_new_byname (names[i]);
    P4EST_ASSERT (conn != NULL);

    for (p4est_topidx_t t = 0; t < conn->num_trees; t++) {
      for (int n = 0; n <= P4EST_DIM; n++) {
        p4est_connect_type_t neigh_type = neigh_types[n];

        for (int j = 0; j < neigh_sizes[n]; j++) {
          sc_array_t          neigh_transforms;

          P4EST_GLOBAL_LDEBUGF ("Testing connectivity %s, tree %d, %s %d\n",
                                names[i], t, neigh_type_names[n], j);
          p4est_log_indent_push ();

          sc_array_init (&neigh_transforms,
                         sizeof (p4est_neighbor_transform_t));
          p4est_connectivity_get_neighbor_transforms (conn, t, neigh_type, j,
                                                      &neigh_transforms);

          switch (n) {
          case 0:
            test_self_transform (conn, t, j, &neigh_transforms);
            break;
          case 1:
            test_face_transform (conn, t, j, &neigh_transforms);
            break;
#ifdef P4_TO_P8
          case 2:
            test_edge_transform (conn, t, j, &neigh_transforms);
            break;
#endif
          case P4EST_DIM:
            test_corner_transform (conn, t, j, &neigh_transforms);
            break;
          default:
            SC_ABORT_NOT_REACHED ();
          }

          sc_array_reset (&neigh_transforms);
          p4est_log_indent_pop ();
        }
      }
    }

    p4est_connectivity_destroy (conn);
  }
  sc_MPI_Finalize ();
  return 0;
}
