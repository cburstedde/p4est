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
