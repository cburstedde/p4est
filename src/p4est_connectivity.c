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

#include <math.h>

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
  connectivity->vtt_offset[num_vertices] = -1;  /* catch bugs */

  connectivity->vertex_to_tree = P4EST_ALLOC (int32_t, num_vtt);
  P4EST_CHECK_ALLOC (connectivity->vertex_to_tree);

  connectivity->vertex_to_vertex = P4EST_ALLOC (int32_t, num_vtt);
  P4EST_CHECK_ALLOC (connectivity->vertex_to_vertex);

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
  P4EST_FREE (connectivity->vertex_to_vertex);

  P4EST_FREE (connectivity);
}

int
p4est_connectivity_is_valid (p4est_connectivity_t * connectivity)
{
  int                 found;
  int8_t              face, rface, nface, orientation;
  int8_t              corner;
  int32_t             tree, ntree, ctree;
  int32_t             vertex, cvertex, corner_trees;
  int32_t             v1, v2, w1, w2;
  const int32_t       num_trees = connectivity->num_trees;
  const int32_t       num_vertices = connectivity->num_vertices;
  const int32_t       num_vtt = connectivity->vtt_offset[num_vertices];
  const int32_t      *ttv = connectivity->tree_to_vertex;
  const int32_t      *ttt = connectivity->tree_to_tree;
  const int8_t       *ttf = connectivity->tree_to_face;
  const int32_t      *vtt = connectivity->vertex_to_tree;
  const int32_t      *vtv = connectivity->vertex_to_vertex;
  const int32_t      *voff = connectivity->vtt_offset;

  if (num_trees < 1 || num_vertices < 4) {
    fprintf (stderr, "Invalid numbers of trees or vertices");
    return 0;
  }

  for (tree = 0; tree < num_trees; ++tree) {
    for (face = 0; face < 4; ++face) {
      ntree = ttt[tree * 4 + face];
      if (ntree < 0 || ntree >= num_trees) {
        fprintf (stderr, "Tree range A in %d %d\n", tree, face);
        return 0;
      }
      rface = ttf[tree * 4 + face];
      if (rface < 0 || rface >= 8) {
        fprintf (stderr, "Face range in %d %d\n", tree, face);
        return 0;
      }
      nface = (int8_t) (rface % 4);     /* clamp to a real face index */
      orientation = (int8_t) (rface / 4);       /* 0 (same) or 1 (opposite) */
      if (ntree == tree) {
        /* no neighbor across this face or self-periodic */
        if (nface == face && orientation != 0) {
          fprintf (stderr, "Face invalid in %d %d\n", tree, face);
          return 0;
        }
      }
      if (ntree != tree || nface != face) {
        /* check reciprocity */
        if (ttt[ntree * 4 + nface] != tree) {
          fprintf (stderr, "Tree reciprocity in %d %d\n", tree, face);
          return 0;
        }
        if (ttf[ntree * 4 + nface] != face + 4 * orientation) {
          fprintf (stderr, "Face reciprocity in %d %d\n", tree, face);
          return 0;
        }

        /* a neighbor across this face */
        v1 = ttv[tree * 4 + face];
        v2 = ttv[tree * 4 + (face + 1) % 4];
        w1 = ttv[ntree * 4 + nface];
        w2 = ttv[ntree * 4 + (nface + 1) % 4];
        if (v1 == v2 || w1 == w2) {
          fprintf (stderr, "Vertex invalid in %d %d\n", tree, face);
          return 0;
        }
        if ((v1 == w2 && v2 == w1) && orientation != 0) {
          fprintf (stderr, "Orientation mismatch A in %d %d\n", tree, face);
          return 0;
        }
        if ((v1 == w1 && v2 == w2) && orientation != 1) {
          fprintf (stderr, "Orientation mismatch B in %d %d\n", tree, face);
          return 0;
        }
      }
    }
  }

  for (vertex = 0; vertex < num_vertices; ++vertex) {
    corner_trees = voff[vertex + 1] - voff[vertex];
    if (corner_trees <= 0 || voff[vertex + 1] > num_vtt) {
      fprintf (stderr, "Vertex offset mismatch %d\n", vertex);
      return 0;
    }
    for (ctree = 0; ctree < corner_trees; ++ctree) {
      ntree = vtt[voff[vertex] + ctree];
      if (ntree < 0 || ntree >= num_trees) {
        fprintf (stderr, "Tree range B in %d %d\n", vertex, ntree);
        return 0;
      }
      found = 0;
      for (corner = 0; corner < 4; ++corner) {
        if (ttv[4 * ntree + corner] == vertex) {
          found = 1;
          break;
        }
      }
      if (!found) {
        fprintf (stderr, "Corner mismatch in %d %d\n", vertex, ntree);
        return 0;
      }
      cvertex = vtv[voff[vertex] + ctree];
      if (cvertex < 0 || cvertex >= num_vertices) {
        fprintf (stderr, "Vertex mismatch in %d %d\n", vertex, ntree);
        return 0;
      }
    }
  }

  for (tree = 0; tree < num_trees; ++tree) {
    for (corner = 0; corner < 4; ++corner) {
      vertex = ttv[tree * 4 + corner];
      corner_trees = voff[vertex + 1] - voff[vertex];
      found = 0;
      for (ctree = 0; ctree < corner_trees; ++ctree) {
        if (vtt[voff[vertex] + ctree] == tree &&
            vtv[voff[vertex] + ctree] == vertex) {
          ++found;
        }
      }
      if (found != 1) {
        fprintf (stderr, "Tree and vertex count in %d %d\n", tree, corner);
        return 0;
      }
    }
  }

  return 1;
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

  connectivity->vertex_to_vertex[0] = 0;
  connectivity->vertex_to_vertex[1] = 1;
  connectivity->vertex_to_vertex[2] = 2;
  connectivity->vertex_to_vertex[3] = 3;

  P4EST_ASSERT (p4est_connectivity_is_valid (connectivity));

  return connectivity;
}

p4est_connectivity_t *
p4est_connectivity_new_corner (void)
{
  p4est_connectivity_t *connectivity;
  const int32_t       num_trees = 3;
  const int32_t       num_vertices = 7;
  const int32_t       num_vtt = 12;
  const int32_t       tree_to_vertex[3 * 4] = {
    0, 1, 3, 2, 0, 2, 5, 6, 2, 3, 4, 5,
  };
  const int32_t       tree_to_tree[3 * 4] = {
    0, 0, 2, 1, 0, 2, 1, 1, 0, 2, 2, 1,
  };
  const int8_t        tree_to_face[3 * 4] = {
    0, 1, 0, 0, 3, 3, 2, 3, 2, 1, 2, 1,
  };
  const double        vertices[7 * 3] = {
    -1, -1, 0,
    0, -1, 0,
    0, 0, 1,
    1, 0, 1,
    1, 1, 1,
    0, 1, 1,
    -1, 0, 0,
  };
  const int32_t       vtt_offset[7 + 1] = {
    0, 2, 3, 6, 8, 9, 11, 12,
  };
  const int32_t       vertex_to_tree[12] = {
    0, 1, 0, 0, 2, 1, 0, 2, 2, 1, 2, 1,
  };
  const int32_t       vertex_to_vertex[12] = {
    0, 0, 1, 2, 2, 2, 3, 3, 4, 5, 5, 6,
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
          sizeof (int32_t) * (num_vertices + 1));
  memcpy (connectivity->vertex_to_tree, vertex_to_tree,
          sizeof (int32_t) * num_vtt);
  memcpy (connectivity->vertex_to_vertex, vertex_to_vertex,
          sizeof (int32_t) * num_vtt);

  P4EST_ASSERT (p4est_connectivity_is_valid (connectivity));

  return connectivity;
}

p4est_connectivity_t *
p4est_connectivity_new_moebius (void)
{
  p4est_connectivity_t *connectivity;
  const int32_t       num_trees = 5;
  const int32_t       num_vertices = 10;
  const int32_t       num_vtt = 20;
  const double        halfsqrt3 = .5 * sqrt (3.);
  const double        vertices[10 * 3] = {
    0, 0, 0,
    0, 1, 0,
    1, 0, 0,
    1, 1, 0,
    1.5, 0, halfsqrt3,
    1.5, 1, halfsqrt3,
    .5, .5, 1.5,
    .5, .5, 2,
    -.5, 0, halfsqrt3,
    -.5, 1, halfsqrt3,
  };
  const int32_t       tree_to_vertex[5 * 4] = {
    0, 2, 3, 1,
    3, 5, 4, 2,                 /* left-handed */
    4, 6, 7, 5,
    6, 7, 8, 9,
    9, 8, 0, 1,
  };
  const int32_t       tree_to_tree[5 * 4] = {
    0, 1, 0, 4,
    1, 2, 1, 0,
    2, 3, 2, 1,
    2, 3, 4, 3,
    3, 4, 0, 4,
  };
  const int8_t        tree_to_face[5 * 4] = {
    0, 7, 2, 2,
    0, 7, 2, 5,
    0, 4, 2, 5,
    5, 1, 0, 3,
    2, 1, 3, 3,
  };
  const int32_t       vtt_offset[10 + 1] = {
    0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20,
  };
  const int32_t       vertex_to_tree[20] = {
    0, 4, 0, 4,
    0, 1, 0, 1,
    1, 2, 1, 2,
    2, 3, 2, 3,
    3, 4, 3, 4,
  };
  const int32_t       vertex_to_vertex[20] = {
    0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9,
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
          sizeof (int32_t) * (num_vertices + 1));
  memcpy (connectivity->vertex_to_tree, vertex_to_tree,
          sizeof (int32_t) * num_vtt);
  memcpy (connectivity->vertex_to_vertex, vertex_to_vertex,
          sizeof (int32_t) * num_vtt);

  P4EST_ASSERT (p4est_connectivity_is_valid (connectivity));

  return connectivity;
}

p4est_connectivity_t *
p4est_connectivity_new_star (void)
{
  int                 i;
  p4est_connectivity_t *connectivity;
  const double        r1 = 1.;
  const double        r2 = 1.5;
  const int32_t       num_trees = 6;
  const int32_t       num_vertices = 13;
  const int32_t       num_vtt = 24;
  const int32_t       tree_to_vertex[6 * 4] = {
    0, 1, 2, 3,
    0, 3, 4, 5,
    5, 6, 7, 0,
    8, 7, 0, 9,                 /* left-handed */
    9, 0, 11, 10,               /* left-handed */
    12, 1, 0, 11,
  };
  const int32_t       tree_to_tree[6 * 4] = {
    5, 0, 0, 1,
    0, 1, 1, 2,
    2, 2, 3, 1,
    3, 2, 4, 3,
    3, 5, 4, 4,
    5, 0, 4, 5,
  };
  const int8_t        tree_to_face[6 * 4] = {
    1, 1, 2, 0,
    3, 1, 2, 3,
    0, 1, 5, 3,
    0, 6, 0, 3,
    2, 6, 2, 3,
    0, 0, 5, 3,
  };
  const int32_t       vtt_offset[13 + 1] = {
    0, 6, 8, 9, 11, 12, 14, 15, 17, 18, 20, 21, 23, 24,
  };
  const int32_t       vertex_to_tree[24] = {
    0, 1, 2, 3, 4, 5, 5, 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5,
  };
  const int32_t       vertex_to_vertex[24] = {
    0, 0, 0, 0, 0, 0, 1, 1, 2, 3, 3, 4,
    5, 5, 6, 7, 7, 8, 9, 9, 10, 11, 11, 12,
  };

  connectivity = p4est_connectivity_new (num_trees, num_vertices, num_vtt);

  connectivity->vertices[0 * 3 + 0] = 0;
  connectivity->vertices[0 * 3 + 1] = 0;
  connectivity->vertices[0 * 3 + 2] = 0;
  for (i = 0; i < 6; ++i) {
    connectivity->vertices[(2 * i + 1) * 3 + 0] = r1 * cos (i * M_PI / 3);
    connectivity->vertices[(2 * i + 1) * 3 + 1] = r1 * sin (i * M_PI / 3);
    connectivity->vertices[(2 * i + 1) * 3 + 2] = 0;
    connectivity->vertices[(2 * i + 2) * 3 + 0] =
      r2 * cos ((i + .5) * M_PI / 3);
    connectivity->vertices[(2 * i + 2) * 3 + 1] =
      r2 * sin ((i + .5) * M_PI / 3);
    connectivity->vertices[(2 * i + 2) * 3 + 2] = 0;
  }

  memcpy (connectivity->tree_to_vertex, tree_to_vertex,
          sizeof (int32_t) * 4 * num_trees);
  memcpy (connectivity->tree_to_tree, tree_to_tree,
          sizeof (int32_t) * 4 * num_trees);
  memcpy (connectivity->tree_to_face, tree_to_face,
          sizeof (int8_t) * 4 * num_trees);
  memcpy (connectivity->vtt_offset, vtt_offset,
          sizeof (int32_t) * (num_vertices + 1));
  memcpy (connectivity->vertex_to_tree, vertex_to_tree,
          sizeof (int32_t) * num_vtt);
  memcpy (connectivity->vertex_to_vertex, vertex_to_vertex,
          sizeof (int32_t) * num_vtt);

  P4EST_ASSERT (p4est_connectivity_is_valid (connectivity));

  return connectivity;
}

p4est_connectivity_t *
p4est_connectivity_new_periodic (void)
{
  p4est_connectivity_t *connectivity;
  const int32_t       num_trees = 1;
  const int32_t       num_vertices = 4;
  const int32_t       num_vtt = 12;
  const int32_t       tree_to_vertex[1 * 4] = {
    0, 1, 2, 3,
  };
  const int32_t       tree_to_tree[1 * 4] = {
    0, 0, 0, 0,
  };
  const int8_t        tree_to_face[1 * 4] = {
    6, 3, 4, 1,
  };
  const double        vertices[4 * 3] = {
    0, 0, 0,
    1, 0, 0,
    1, 1, 0,
    0, 1, 0,
  };
  const int32_t       vtt_offset[4 + 1] = {
    0, 3, 6, 9, 12,
  };
  const int32_t       vertex_to_tree[12] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  };
  const int32_t       vertex_to_vertex[12] = {
    0, 1, 2, 0, 1, 3, 0, 2, 3, 1, 2, 3,
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
          sizeof (int32_t) * (num_vertices + 1));
  memcpy (connectivity->vertex_to_tree, vertex_to_tree,
          sizeof (int32_t) * num_vtt);
  memcpy (connectivity->vertex_to_vertex, vertex_to_vertex,
          sizeof (int32_t) * num_vtt);

  P4EST_ASSERT (p4est_connectivity_is_valid (connectivity));

  return connectivity;
}

int
p4est_find_face_transform (p4est_connectivity_t * connectivity,
                           int32_t itree, int8_t iface)
{
  int8_t              nrface, neighbor_face, orientation;
  int32_t             neighbor_tree;

  P4EST_ASSERT (itree >= 0 && itree < connectivity->num_trees);
  P4EST_ASSERT (iface >= 0 && iface < 4);

  neighbor_tree = connectivity->tree_to_tree[4 * itree + iface];
  P4EST_ASSERT (neighbor_tree >= 0 &&
                neighbor_tree < connectivity->num_trees);

  nrface = connectivity->tree_to_face[4 * itree + iface];
  P4EST_ASSERT (nrface >= 0 && nrface < 8);
  neighbor_face = (int8_t) (nrface % 4);
  orientation = (int8_t) (nrface / 4);

  if (neighbor_tree == itree && neighbor_face == iface) {
    return -1;
  }

  return p4est_transform_table[iface][neighbor_face][orientation];
}

void
p4est_find_corner_info (p4est_connectivity_t * conn,
                        int32_t itree, int8_t icorner,
                        p4est_array_t * corner_info)
{
  int                 incount, num_found;
  int8_t              ncorner;
  int32_t             corner_trees, ctree;
  int32_t             ntree, ntree1, ntree2;
  int32_t             ivertex, nvertex;
  p4est_corner_info_t *ci;

  P4EST_ASSERT (0 <= itree && itree < conn->num_trees);
  P4EST_ASSERT (0 <= icorner && icorner < 4);
  P4EST_ASSERT (corner_info->elem_size == sizeof (p4est_corner_info_t));
  incount = corner_info->elem_count;

  ivertex = conn->tree_to_vertex[4 * itree + icorner];
  P4EST_ASSERT (0 <= ivertex && ivertex < conn->num_vertices);

  ntree1 = conn->tree_to_tree[4 * itree + (icorner + 3) % 4];
  ntree2 = conn->tree_to_tree[4 * itree + icorner];

  corner_trees = conn->vtt_offset[ivertex + 1] - conn->vtt_offset[ivertex];
  num_found = 0;
  for (ctree = 0; ctree < corner_trees; ++ctree) {
    ntree = conn->vertex_to_tree[conn->vtt_offset[ivertex] + ctree];
    nvertex = conn->vertex_to_vertex[conn->vtt_offset[ivertex] + ctree];
    if (nvertex == ivertex &&
        (ntree == itree || ntree == ntree1 || ntree == ntree2)) {
      continue;
    }
    /* else we have a true corner with ntree, find its id */
    for (ncorner = 0; ncorner < 4; ++ncorner) {
      if (nvertex == conn->tree_to_vertex[4 * ntree + ncorner]) {
        break;
      }
    }
    P4EST_ASSERT (ncorner < 4);
    if (num_found >= incount) {
      p4est_array_resize (corner_info, num_found + 1);
    }
    ci = p4est_array_index (corner_info, num_found);
    ci->ntree = ntree;
    ci->ncorner = ncorner;
    ++num_found;
  }
  p4est_array_resize (corner_info, num_found);
  P4EST_ASSERT (corner_trees == num_found +
                1 + (ntree1 != itree) + (ntree2 != itree));
}

/* EOF p4est_connectivity.c */
