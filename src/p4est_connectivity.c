/*
  This file is part of p4est.
  p4est is a C library to manage a parallel collection of quadtrees and/or
  octrees.

  Copyright (C) 2007 Carsten Burstedde, Lucas Wilcox.

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

#include <p4est_connectivity.h>
#include <p4est_base.h>
#include <p4est_memory.h>

/* *INDENT-OFF* */

const int p4est_transform_table[4][4][2] =
{{{2, 4}, {3, 5}, {0, 6}, {1, 7}},
 {{1, 5}, {2, 6}, {3, 7}, {0, 4}},
 {{0, 6}, {1, 7}, {2, 4}, {3, 5}},
 {{3, 7}, {0, 4}, {1, 5}, {2, 6}}};

/* *INDENT-ON* */

p4est_connectivity_t *
p4est_connectivity_new (int32_t num_trees, int32_t num_vertices,
                        int32_t num_vtt)
{
  p4est_connectivity_t *connectivity;

  connectivity = P4EST_ALLOC_ZERO (p4est_connectivity_t, 1);
  P4EST_CHECK_ALLOC (connectivity);

  connectivity->num_trees = num_trees;
  connectivity->num_vertices = num_vertices;

  connectivity->tree_to_vertex = P4EST_ALLOC (int32_t, 4 * num_trees);
  P4EST_CHECK_ALLOC (connectivity->tree_to_vertex);

  connectivity->tree_to_tree = P4EST_ALLOC (int32_t, 4 * num_trees);
  P4EST_CHECK_ALLOC (connectivity->tree_to_tree);

  connectivity->tree_to_face = P4EST_ALLOC (int8_t, 4 * num_trees);
  P4EST_CHECK_ALLOC (connectivity->tree_to_face);

  connectivity->vertices = P4EST_ALLOC (double, 3 * num_vertices);
  P4EST_CHECK_ALLOC (connectivity->vertices);

  connectivity->vtt_offset = P4EST_ALLOC (int32_t, num_vertices + 1);
  P4EST_CHECK_ALLOC (connectivity->vtt_offset);

  connectivity->vertex_to_tree = P4EST_ALLOC (int32_t, num_vtt);
  P4EST_CHECK_ALLOC (connectivity->vertex_to_tree);

  return connectivity;
}

void
p4est_connectivity_destroy (p4est_connectivity_t * connectivity)
{
  P4EST_FREE (connectivity->tree_to_face);
  P4EST_FREE (connectivity->tree_to_tree);
  P4EST_FREE (connectivity->tree_to_vertex);
  P4EST_FREE (connectivity->vertices);
  P4EST_FREE (connectivity->vtt_offset);
  P4EST_FREE (connectivity->vertex_to_tree);

  P4EST_FREE (connectivity);
}

p4est_connectivity_t *
p4est_connectivity_new_unitsquare (void)
{
  p4est_connectivity_t *connectivity;

  connectivity = p4est_connectivity_new (1, 4, 4);

  /* assign vertex numbers */
  connectivity->tree_to_vertex[0] = 0;
  connectivity->tree_to_vertex[1] = 1;
  connectivity->tree_to_vertex[2] = 2;
  connectivity->tree_to_vertex[3] = 3;

  /* we do not have neighbors, put in ourself */
  connectivity->tree_to_tree[0] = 0;
  connectivity->tree_to_tree[1] = 0;
  connectivity->tree_to_tree[2] = 0;
  connectivity->tree_to_tree[3] = 0;

  /* we do not share faces, put in our own face */
  connectivity->tree_to_face[0] = 0;
  connectivity->tree_to_face[1] = 1;
  connectivity->tree_to_face[2] = 2;
  connectivity->tree_to_face[3] = 3;

  /* This is a square from [0,0] to [1,1].
   * Note that the z component is zero.
   */
  /* x=0 y=0 z=0 */
  connectivity->vertices[0 * 3 + 0] = 0.0;
  connectivity->vertices[0 * 3 + 1] = 0.0;
  connectivity->vertices[0 * 3 + 2] = 0.0;

  /* x=1 y=0 z=0 */
  connectivity->vertices[1 * 3 + 0] = 1.0;
  connectivity->vertices[1 * 3 + 1] = 0.0;
  connectivity->vertices[1 * 3 + 2] = 0.0;

  /* x=1 y=1 z=0 */
  connectivity->vertices[2 * 3 + 0] = 1.0;
  connectivity->vertices[2 * 3 + 1] = 1.0;
  connectivity->vertices[2 * 3 + 2] = 0.0;

  /* x=0 y=1 z=0 */
  connectivity->vertices[3 * 3 + 0] = 0.0;
  connectivity->vertices[3 * 3 + 1] = 1.0;
  connectivity->vertices[3 * 3 + 2] = 0.0;

  connectivity->vtt_offset[0] = 0;
  connectivity->vtt_offset[1] = 1;
  connectivity->vtt_offset[2] = 2;
  connectivity->vtt_offset[3] = 3;
  connectivity->vtt_offset[4] = 4;

  connectivity->vertex_to_tree[0] = 0;
  connectivity->vertex_to_tree[1] = 0;
  connectivity->vertex_to_tree[2] = 0;
  connectivity->vertex_to_tree[3] = 0;

  return connectivity;
}

p4est_connectivity_t *
p4est_connectivity_new_corner (void)
{
  p4est_connectivity_t *connectivity;
  const int32_t       num_trees = 3;
  const int32_t       num_vertices = 7;
  const int32_t       num_vtt = 12;
  const int32_t       tree_to_vertex[] = {
    0, 1, 3, 2, 0, 2, 5, 6, 2, 3, 4, 5
  };
  const int32_t       tree_to_tree[] = {
    0, 0, 2, 1, 0, 2, 1, 1, 0, 2, 2, 1
  };
  const int8_t        tree_to_face[] = {
    0, 1, 0, 0, 3, 3, 2, 3, 2, 1, 2, 1
  };
  const double        vertices[] = {
    -1.00000000000e+00, -1.00000000000e+00, 0.00000000000e+00,
    0.00000000000e+00, -1.00000000000e+00, 0.00000000000e+00,
    0.00000000000e+00, 0.00000000000e+00, 1.00000000000e+00,
    1.00000000000e+00, 0.00000000000e+00, 1.00000000000e+00,
    1.00000000000e+00, 1.00000000000e+00, 1.00000000000e+00,
    0.00000000000e+00, 1.00000000000e+00, 1.00000000000e+00,
    -1.00000000000e+00, 0.00000000000e+00, 0.00000000000e+00
  };
  const int32_t       vtt_offset[] = {
    0, 2, 3, 6, 8, 9, 11, 12
  };
  const int32_t       vertex_to_tree[] = {
    0, 1, 0, 0, 2, 1, 0, 2, 2, 1, 2, 1
  };

  connectivity = p4est_connectivity_new (num_trees, num_vertices, num_vtt);

  memcpy (connectivity->tree_to_vertex, tree_to_vertex,
          sizeof (int32_t) * 4 * num_trees);
  memcpy (connectivity->tree_to_tree, tree_to_tree,
          sizeof (int32_t) * 4 * num_trees);
  memcpy (connectivity->tree_to_face, tree_to_face,
          sizeof (int8_t) * 4 * num_trees);
  memcpy (connectivity->vertices, vertices,
          sizeof (double) * 3 * num_vertices);
  memcpy (connectivity->vtt_offset, vtt_offset,
          sizeof (int32_t) * num_vertices + 1);
  memcpy (connectivity->vertex_to_tree, vertex_to_tree,
          sizeof (int32_t) * num_vtt);

  return connectivity;
}

int
p4est_find_face_transform (p4est_connectivity_t * connectivity,
                           int32_t itree, int8_t face)
{
  int                 orientation;
  int8_t              neighbor_face;
  int32_t             neighbor_tree;
  int32_t             v1, v2;
  int32_t             n1, n2, w1, w2;

  P4EST_ASSERT (itree >= 0 && itree < connectivity->num_trees);
  P4EST_ASSERT (face >= 0 && face < 4);

  neighbor_tree = connectivity->tree_to_tree[4 * itree + face];
  P4EST_ASSERT (neighbor_tree >= 0 &&
                neighbor_tree < connectivity->num_trees);
  if (neighbor_tree == itree) {
    return -1;
  }

  /* look up my own vertices */
  neighbor_face = connectivity->tree_to_face[4 * itree + face];
  P4EST_ASSERT (neighbor_face >= 0 && neighbor_face < 4);
  v1 = connectivity->tree_to_vertex[4 * itree + face];
  v2 = connectivity->tree_to_vertex[4 * itree + (face + 1) % 4];

  /* look up my neighbor's vertices */
  n1 = 4 * neighbor_tree + neighbor_face;
  n2 = 4 * neighbor_tree + (neighbor_face + 1) % 4;
  w1 = connectivity->tree_to_vertex[n1];
  w2 = connectivity->tree_to_vertex[n2];

  /* derive orientation from vertex relations */
  P4EST_ASSERT (v1 != v2 && w1 != w2);
  P4EST_ASSERT ((w1 == v1 && w2 == v2) || (w1 == v2 && w2 == v1));
  orientation = (w1 == v1 && w2 == v2) ? 1 : 0;

  return p4est_transform_table[face][neighbor_face][orientation];
}

/* EOF p4est_connectivity.c */
