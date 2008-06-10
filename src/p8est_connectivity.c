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
  int                 face, rface, nface, orientation; /*, corner; */
  p4est_topidx_t      tree, ntree; /*, ctree; */
#if 0
  p4est_topidx_t      vertex, cvertex, corner_trees;
  p4est_topidx_t      v1, v2, w1, w2;
#endif
  p4est_topidx_t      nvtt;
  const p4est_topidx_t num_trees = conn->num_trees;
  const p4est_topidx_t num_vertices = conn->num_vertices;
  const p4est_topidx_t num_vtt = conn->vtt_offset[num_vertices];
#if 0
  const p4est_topidx_t *ttv = conn->tree_to_vertex;
#endif
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
      nface = rface % 4;        /* clamp to a real face index */
      orientation = rface / 4;  /* 0..3 for relative rotation */
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
        /* TODO: everything below here needs to be adapted for 3D */
#if 0
        if ((int) ttf[ntree * 6 + nface] != face + 4 * orientation) {
          fprintf (stderr, "Face reciprocity in %lld %d\n",
                   (long long) tree, face);
          return false;
        }

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
      fprintf (stderr, "Vertex to tree %d out of range", nvtt);
    }
    if (vtv[nvtt] < 0 || vtv[nvtt] >= num_vertices) {
      fprintf (stderr, "Vertex to vertex %d out of range", nvtt);
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
    0, 0, 1,
    0, 1, 0,
    0, 1, 1,
    1, 0, 0,
    1, 0, 1,
    1, 1, 0,
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
  return NULL;
}

/* EOF p8est_connectivity.h */
