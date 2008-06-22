/*
  This file is part of p4est.
  p4est is a C library to manage a parallel collection of quadtrees and/or
  octrees.

  Copyright (C) 2008 Carsten Burstedde, Lucas Wilcox.

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

#include <p8est_connectivity.h>

/* *INDENT-OFF* */
const int           p8est_face_vertices[6][4] =
{{ 0, 2, 4, 6 },
 { 1, 3, 5, 7 },
 { 0, 1, 4, 5 },
 { 2, 3, 6, 7 },
 { 0, 1, 2, 3 },
 { 4, 5, 6, 7 }};

const int           p8est_face_edges[6][4] =
{{ 4, 6,  8, 10 },
 { 5, 7,  9, 11 },
 { 0, 2,  8,  9 },
 { 1, 3, 10, 11 },
 { 0, 1,  4,  5 },
 { 2, 3,  6,  7 }};

const int           p8est_face_permutations[8][4] =
{{ 0, 1, 2, 3 },                /* no.  0 of 0..23 */
 { 0, 2, 1, 3 },                /* no.  2 of 0..23 */
 { 1, 0, 3, 2 },                /* no.  7 of 0..23 */
 { 1, 3, 0, 2 },                /* no. 10 of 0..23 */
 { 2, 0, 3, 1 },                /* no. 13 of 0..23 */
 { 2, 3, 0, 1 },                /* no. 16 of 0..23 */
 { 3, 1, 2, 0 },                /* no. 21 of 0..23 */
 { 3, 2, 1, 0 }};               /* no. 23 of 0..23 */

const int           p8est_face_permutation_sets[3][4] =
{{ 1, 2, 5, 6 },
 { 0, 3, 4, 7 },
 { 0, 4, 3, 7 }};

const int           p8est_face_permutation_refs[6][6] =
{{ 0, 1, 1, 0, 0, 1 },
 { 2, 0, 0, 1, 1, 0 },
 { 2, 0, 0, 1, 1, 0 },
 { 0, 2, 2, 0, 0, 1 },
 { 0, 2, 2, 0, 0, 1 },
 { 2, 0, 0, 2, 2, 0 }};

const int           p8est_edge_vertices[12][2] =
{{ 0, 1 },
 { 2, 3 },
 { 4, 5 },
 { 6, 7 },
 { 0, 2 },
 { 1, 3 },
 { 4, 6 },
 { 5, 7 },
 { 0, 4 },
 { 1, 5 },
 { 2, 6 },
 { 3, 7 }};
/* *INDENT-ON* */

static p8est_connectivity_t *
p8est_connectivity_new_copy (p4est_topidx_t num_trees,
                             p4est_topidx_t num_vertices,
                             p4est_topidx_t num_vtt,
                             const p4est_topidx_t * ttv,
                             const p4est_topidx_t * ttt,
                             const int8_t * ttf,
                             const double *vertices,
                             const p4est_topidx_t * voff,
                             const p4est_topidx_t * vtt,
                             const p4est_topidx_t * vtv)
{
  const bool          alloc_vertices = (vertices != NULL);
  p8est_connectivity_t *conn;

  conn = p8est_connectivity_new (num_trees, num_vertices, num_vtt,
                                 alloc_vertices);

  memcpy (conn->tree_to_vertex, ttv, sizeof (p4est_topidx_t) * 8 * num_trees);
  memcpy (conn->tree_to_tree, ttt, sizeof (p4est_topidx_t) * 6 * num_trees);
  memcpy (conn->tree_to_face, ttf, sizeof (int8_t) * 6 * num_trees);
  if (alloc_vertices) {
    memcpy (conn->vertices, vertices, sizeof (double) * 3 * num_vertices);
  }
  memcpy (conn->vtt_offset, voff,
          sizeof (p4est_topidx_t) * (num_vertices + 1));
  memcpy (conn->vertex_to_tree, vtt, sizeof (p4est_topidx_t) * num_vtt);
  memcpy (conn->vertex_to_vertex, vtv, sizeof (p4est_topidx_t) * num_vtt);

  P4EST_ASSERT (p8est_connectivity_is_valid (conn));

  return conn;
}

p8est_connectivity_t *
p8est_connectivity_new (p4est_topidx_t num_trees, p4est_topidx_t num_vertices,
                        p4est_topidx_t num_vtt, bool alloc_vertices)
{
  p8est_connectivity_t *conn;

  conn = P4EST_ALLOC_ZERO (p8est_connectivity_t, 1);

  conn->num_trees = num_trees;
  conn->num_vertices = num_vertices;

  conn->tree_to_vertex = P4EST_ALLOC (p4est_topidx_t, 8 * num_trees);
  conn->tree_to_tree = P4EST_ALLOC (p4est_topidx_t, 6 * num_trees);
  conn->tree_to_face = P4EST_ALLOC (int8_t, 6 * num_trees);

  if (alloc_vertices)
    conn->vertices = P4EST_ALLOC (double, 3 * num_vertices);
  else
    conn->vertices = NULL;

  conn->vtt_offset = P4EST_ALLOC (p4est_topidx_t, num_vertices + 1);
  conn->vtt_offset[num_vertices] = -1;  /* catch bugs */
  conn->vertex_to_tree = P4EST_ALLOC (p4est_topidx_t, num_vtt);
  conn->vertex_to_vertex = P4EST_ALLOC (p4est_topidx_t, num_vtt);

  return conn;
}

void
p8est_connectivity_destroy (p8est_connectivity_t * conn)
{
  P4EST_FREE (conn->tree_to_vertex);
  P4EST_FREE (conn->tree_to_tree);
  P4EST_FREE (conn->tree_to_face);
  P4EST_FREE (conn->vertices);
  P4EST_FREE (conn->vtt_offset);
  P4EST_FREE (conn->vertex_to_tree);
  P4EST_FREE (conn->vertex_to_vertex);

  P4EST_FREE (conn);
}

bool
p8est_connectivity_is_valid (p8est_connectivity_t * conn)
{
#if 0
  int                 num_found;
#endif
  int                 face, rface, nface, orientation;  /*, corner; */
  int                 fvert, face_ref, face_perm, nvert;
  p4est_topidx_t      tree, ntree;      /*, ctree; */
  p4est_topidx_t      vertex, nvertex;
#if 0
  /* , cvertex, corner_trees; */
  p4est_topidx_t      v1, v2, w1, w2;
#endif
  p4est_topidx_t      nvtt;
  const p4est_topidx_t num_trees = conn->num_trees;
  const p4est_topidx_t num_vertices = conn->num_vertices;
  const p4est_topidx_t num_vtt = conn->vtt_offset[num_vertices];
  const p4est_topidx_t *ttv = conn->tree_to_vertex;
  const p4est_topidx_t *ttt = conn->tree_to_tree;
  const int8_t       *ttf = conn->tree_to_face;
  const p4est_topidx_t *vtt = conn->vertex_to_tree;
  const p4est_topidx_t *vtv = conn->vertex_to_vertex;
#if 0
  const p4est_topidx_t *voff = conn->vtt_offset;
#endif

  if (num_trees < 1 || num_vertices < 8) {
    fprintf (stderr, "Invalid numbers of trees or vertices");
    return false;
  }

  for (tree = 0; tree < num_trees; ++tree) {
    for (face = 0; face < 6; ++face) {
      ntree = ttt[tree * 6 + face];
      if (ntree < 0 || ntree >= num_trees) {
        fprintf (stderr, "Tree range A in %lld %d\n", (long long) tree, face);
        return false;
      }
      rface = (int) ttf[tree * 6 + face];
      if (rface < 0 || rface >= 24) {
        fprintf (stderr, "Face range in %lld %d\n", (long long) tree, face);
        return false;
      }
      nface = rface % 6;        /* clamp to a real face index */
      orientation = rface / 6;  /* 0..3 for relative rotation */
      if (ntree == tree) {
        /* no neighbor across this face or self-periodic */
        if (nface == face && orientation != 0) {
          fprintf (stderr, "Face invalid in %lld %d\n",
                   (long long) tree, face);
          return false;
        }
      }
      if (ntree != tree || nface != face) {
        /* check reciprocity */
        if (ttt[ntree * 6 + nface] != tree) {
          fprintf (stderr, "Tree reciprocity in %lld %d\n",
                   (long long) tree, face);
          return false;
        }
        if ((int) ttf[ntree * 6 + nface] != face + 6 * orientation) {
          fprintf (stderr, "Face reciprocity in %lld %d\n",
                   (long long) tree, face);
          return false;
        }

        /* check vertex consistency across faces */
        for (fvert = 0; fvert < 4; ++fvert) {
          vertex = ttv[tree * 8 + p8est_face_vertices[face][fvert]];
          if (vertex < 0 || vertex >= num_vertices) {
            fprintf (stderr, "Invalid vertex in %lld %d %d\n",
                     (long long) tree, face, fvert);
            return false;
          }
          for (nvert = fvert + 1; nvert < 4; ++nvert) {
            nvertex = ttv[tree * 8 + p8est_face_vertices[face][nvert]];
            if (vertex == nvertex) {
              fprintf (stderr, "Duplicate vertex in %lld %d %d %d",
                       (long long) tree, face, fvert, nvert);
              return false;
            }
          }
          face_ref = p8est_face_permutation_refs[face][nface];
          face_perm = p8est_face_permutation_sets[face_ref][orientation];
          nvert = p8est_face_permutations[face_perm][fvert];
          nvertex = ttv[ntree * 8 + p8est_face_vertices[nface][nvert]];
#if 0
          P4EST_VERBOSEF ("Me %lld %d %d %lld Other %lld %d %d %lld",
                          (long long) tree, face, fvert,
                          (long long) vertex,
                          (long long) ntree, nface, nvert,
                          (long long) nvertex);
          P4EST_VERBOSEF (" Perm %d %d %d\n",
                          orientation, face_ref, face_perm);
#endif
          if (ntree != tree && nvertex != vertex) {
            fprintf (stderr, "Vertex reciprocity in %lld %d %d\n",
                     (long long) tree, face, fvert);
            return false;
          }
        }

        /* TODO: everything below here needs to be adapted for 3D */
#if 0
        /* a neighbor across this face */
        v1 = ttv[tree * 8 + face];
        v2 = ttv[tree * 8 + (face + 1) % 4];
        w1 = ttv[ntree * 8 + nface];
        w2 = ttv[ntree * 8 + (nface + 1) % 4];
        if (v1 == v2 || w1 == w2) {
          fprintf (stderr, "Vertex invalid in %lld %d\n",
                   (long long) tree, face);
          return false;
        }
        if ((v1 == w2 && v2 == w1) && orientation != 0) {
          fprintf (stderr, "Orientation mismatch A in %lld %d\n",
                   (long long) tree, face);
          return false;
        }
        if ((v1 == w1 && v2 == w2) && orientation != 1) {
          fprintf (stderr, "Orientation mismatch B in %lld %d\n",
                   (long long) tree, face);
          return false;
        }
#endif
      }
    }
  }

  for (nvtt = 0; nvtt < num_vtt; ++nvtt) {
    if (vtt[nvtt] < 0 || vtt[nvtt] >= num_trees) {
      fprintf (stderr, "Vertex to tree %lld out of range\n",
               (long long) nvtt);
      return false;
    }
    if (vtv[nvtt] < 0 || vtv[nvtt] >= num_vertices) {
      fprintf (stderr, "Vertex to vertex %lld out of range\n",
               (long long) nvtt);
      return false;
    }
  }

  return true;
}

p8est_connectivity_t *
p8est_connectivity_new_unitcube (void)
{
  const p4est_topidx_t num_trees = 1;
  const p4est_topidx_t num_vertices = 8;
  const p4est_topidx_t num_vtt = 8;
  const p4est_topidx_t tree_to_vertex[1 * 8] = {
    0, 1, 2, 3, 4, 5, 6, 7,
  };
  const p4est_topidx_t tree_to_tree[1 * 6] = {
    0, 0, 0, 0, 0, 0,
  };
  const int8_t        tree_to_face[1 * 6] = {
    0, 1, 2, 3, 4, 5,
  };
  const double        vertices[3 * 8] = {
    0, 0, 0,
    1, 0, 0,
    0, 1, 0,
    1, 1, 0,
    0, 0, 1,
    1, 0, 1,
    0, 1, 1,
    1, 1, 1,
  };
  const p4est_topidx_t vtt_offset[8 + 1] = {
    0, 1, 2, 3, 4, 5, 6, 7, 8,
  };
  const p4est_topidx_t vertex_to_tree[8] = {
    0, 0, 0, 0, 0, 0, 0, 0,
  };
  const p4est_topidx_t vertex_to_vertex[8] = {
    0, 1, 2, 3, 4, 5, 6, 7,
  };

  return p8est_connectivity_new_copy (num_trees, num_vertices, num_vtt,
                                      tree_to_vertex, tree_to_tree,
                                      tree_to_face, vertices, vtt_offset,
                                      vertex_to_tree, vertex_to_vertex);
}

p8est_connectivity_t *
p8est_connectivity_new_periodic (void)
{
  const p4est_topidx_t num_trees = 1;
  const p4est_topidx_t num_vertices = 8;
  const p4est_topidx_t num_vtt = 8;
  const p4est_topidx_t tree_to_vertex[1 * 8] = {
    0, 1, 2, 3, 4, 5, 6, 7,
  };
  const p4est_topidx_t tree_to_tree[1 * 6] = {
    0, 0, 0, 0, 0, 0,
  };
  const int8_t        tree_to_face[1 * 6] = {
    1, 0, 2, 3, 11, 10,
  };
  const double        vertices[3 * 8] = {
    0, 0, 0,
    1, 0, 0,
    0, 1, 0,
    1, 1, 0,
    0, 0, 1,
    1, 0, 1,
    0, 1, 1,
    1, 1, 1,
  };
  const p4est_topidx_t vtt_offset[8 + 1] = {
    0, 1, 2, 3, 4, 5, 6, 7, 8,
  };
  const p4est_topidx_t vertex_to_tree[8] = {
    0, 0, 0, 0, 0, 0, 0, 0,
  };
  const p4est_topidx_t vertex_to_vertex[8] = {
    0, 1, 2, 3, 4, 5, 6, 7,
  };

  return p8est_connectivity_new_copy (num_trees, num_vertices, num_vtt,
                                      tree_to_vertex, tree_to_tree,
                                      tree_to_face, vertices, vtt_offset,
                                      vertex_to_tree, vertex_to_vertex);
}

p8est_connectivity_t *
p8est_connectivity_new_rotcubes (void)
{
  const p4est_topidx_t num_trees = 4;
  const p4est_topidx_t num_vertices = 20;
  const p4est_topidx_t num_vtt = 20;
  const p4est_topidx_t tree_to_vertex[4 * 8] = {
    0, 17, 3, 4, 15, 11, 13, 14,
    7, 2, 6, 17, 9, 12, 8, 11,
    2, 12, 5, 10, 17, 11, 4, 14,
    19, 13, 18, 14, 16, 15, 1, 11,
  };
  const p4est_topidx_t tree_to_tree[4 * 6] = {
    0, 2, 0, 0, 0, 3,
    1, 2, 1, 1, 1, 1,
    2, 2, 1, 2, 2, 0,
    3, 0, 3, 3, 3, 3,
  };
  const int8_t        tree_to_face[4 * 6] = {
    0, 5, 2, 3, 4, 13,
    0, 2, 2, 3, 4, 5,
    0, 1, 1, 3, 4, 1,
    0, 17, 2, 3, 4, 5,
  };
  const double        vertices[3 * 20] = {
    0, 0, 0,
    1, 0, 2,
    2, 0, 0,
    0, 1, 0,
    1, 1, 0,
    2, 1, 0,
    1, -1, 0,
    2, -1, 0,
    1, -1, 1,
    2, -1, 1,
    2, 1, 1,
    1, 0, 1,
    2, 0, 1,
    0, 1, 1,
    1, 1, 1,
    0, 0, 1,
    0, 0, 2,
    1, 0, 0,
    1, 1, 2,
    0, 1, 2,
  };

  const p4est_topidx_t vtt_offset[20 + 1] = {
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
  };
  const p4est_topidx_t vertex_to_tree[20] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  };
  const p4est_topidx_t vertex_to_vertex[20] = {
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
  };

  return p8est_connectivity_new_copy (num_trees, num_vertices, num_vtt,
                                      tree_to_vertex, tree_to_tree,
                                      tree_to_face, vertices, vtt_offset,
                                      vertex_to_tree, vertex_to_vertex);
}

p4est_topidx_t
p8est_find_face_transform (p8est_connectivity_t * connectivity,
                           p4est_topidx_t my_tree, int my_face,
                           int my_axis[3], int target_axis[3],
                           int edge_reverse[3])
{
  int                 i;
  int                 target_code, target_face, orientation;
  int                 face_ref, face_perm;
  int                 target_edge_vertex[3];
  p4est_topidx_t      target_tree;

  target_tree = connectivity->tree_to_tree[6 * my_tree + my_face];
  target_code = (int) connectivity->tree_to_face[6 * my_tree + my_face];
  target_face = target_code % 6;
  orientation = target_code / 6;

  P4EST_ASSERT (0 <= my_face && my_face < 6);
  P4EST_ASSERT (0 <= target_face && target_face < 6);
  P4EST_ASSERT (0 <= orientation && orientation < 4);

  if (target_tree == my_tree && target_face == my_face) {
    P4EST_ASSERT (orientation == 0);
    return -1;
  }

  /* find if my edges 0 and 2 are parallel to the x, y, or z-axis */
  my_axis[0] = p8est_face_edges[my_face][0] / 4;
  my_axis[1] = p8est_face_edges[my_face][2] / 4;

  /* find matching target vertices */
  face_ref = p8est_face_permutation_refs[my_face][target_face];
  face_perm = p8est_face_permutation_sets[face_ref][orientation];
  target_edge_vertex[0] =
    p8est_face_vertices[target_face][p8est_face_permutations[face_perm][0]];
  target_edge_vertex[1] =
    p8est_face_vertices[target_face][p8est_face_permutations[face_perm][1]];
  target_edge_vertex[2] =
    p8est_face_vertices[target_face][p8est_face_permutations[face_perm][2]];

  /* find matching target edges */
  target_axis[0] = target_axis[1] = -1;
  edge_reverse[0] = edge_reverse[1] = -1;
  for (i = 0; i < 12; ++i) {
    if (target_edge_vertex[0] == p8est_edge_vertices[i][0] &&
        target_edge_vertex[1] == p8est_edge_vertices[i][1]) {
      P4EST_ASSERT (target_axis[0] == -1);
      target_axis[0] = i / 4;
      edge_reverse[0] = 0;
#ifndef P4EST_DEBUG
      if (target_axis[1] >= 0)
        break;
      continue;
#endif
    }
    if (target_edge_vertex[0] == p8est_edge_vertices[i][1] &&
        target_edge_vertex[1] == p8est_edge_vertices[i][0]) {
      P4EST_ASSERT (target_axis[0] == -1);
      target_axis[0] = i / 4;
      edge_reverse[0] = 1;
#ifndef P4EST_DEBUG
      if (target_axis[1] >= 0)
        break;
      continue;
#endif
    }
    if (target_edge_vertex[0] == p8est_edge_vertices[i][0] &&
        target_edge_vertex[2] == p8est_edge_vertices[i][1]) {
      P4EST_ASSERT (target_axis[1] == -1);
      target_axis[1] = i / 4;
      edge_reverse[1] = 0;
#ifndef P4EST_DEBUG
      if (target_axis[0] >= 0)
        break;
      continue;
#endif
    }
    if (target_edge_vertex[0] == p8est_edge_vertices[i][1] &&
        target_edge_vertex[2] == p8est_edge_vertices[i][0]) {
      P4EST_ASSERT (target_axis[1] == -1);
      target_axis[1] = i / 4;
      edge_reverse[1] = 1;
#ifndef P4EST_DEBUG
      if (target_axis[0] >= 0)
        break;
      continue;
#endif
    }
  }

  /* find what axis is normal to the faces */
  my_axis[2] = my_face / 2;
  target_axis[2] = target_face / 2;
  edge_reverse[2] = 2 * (my_face % 2) + target_face % 2;

#ifdef P4EST_DEBUG
  for (i = 0; i < 3; ++i) {
    P4EST_ASSERT (0 <= my_axis[i] && my_axis[i] < 3);
    P4EST_ASSERT (0 <= target_axis[i] && target_axis[i] < 3);
    P4EST_ASSERT (my_axis[0] != my_axis[1] &&
                  my_axis[0] != my_axis[2] && my_axis[1] != my_axis[2]);
    P4EST_ASSERT (target_axis[0] != target_axis[1] &&
                  target_axis[0] != target_axis[2] &&
                  target_axis[1] != target_axis[2]);
  }
#endif

  return target_tree;
}

/* EOF p8est_connectivity.h */
