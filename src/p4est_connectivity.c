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

/* *INDENT-OFF* */
const int           p4est_corner_to_zorder[5] = { 0, 1, 3, 2, 4 };
const int           p4est_zface_to_rface[4] = { 3, 1, 0, 2 };

const int           p4est_face_dual[4] = { 2, 3, 0, 1 };
const int           p4est_face_child_hang[4][4] =
{{  0,  1, -1, -1 },
 { -1,  0, -1,  1 },
 { -1, -1,  0,  1 },
 {  0, -1,  1, -1 }};

const int           p4est_transform_table[4][4][2] =
{{{2, 4}, {3, 5}, {0, 6}, {1, 7}},
 {{1, 5}, {2, 6}, {3, 7}, {0, 4}},
 {{0, 6}, {1, 7}, {2, 4}, {3, 5}},
 {{3, 7}, {0, 4}, {1, 5}, {2, 6}}};

const int           p4est_hanging_corner[4][2] =
{{ 1, 2 },
 { 0, 3 },
 { 0, 3 },
 { 1, 2 }};
const int           p4est_hanging_face[4][2] =
{{ 0, 3 },
 { 0, 1 },
 { 3, 2 },
 { 1, 2 }};
/* *INDENT-ON* */

p4est_connectivity_t *
p4est_connectivity_new (p4est_topidx_t num_trees, p4est_topidx_t num_vertices,
                        p4est_topidx_t num_vtt, bool alloc_vxyz)
{
  p4est_connectivity_t *connectivity;

  connectivity = P4EST_ALLOC_ZERO (p4est_connectivity_t, 1);

  connectivity->num_trees = num_trees;
  connectivity->num_vertices = num_vertices;

  connectivity->tree_to_vertex = P4EST_ALLOC (p4est_topidx_t, 4 * num_trees);
  connectivity->tree_to_tree = P4EST_ALLOC (p4est_topidx_t, 4 * num_trees);
  connectivity->tree_to_face = P4EST_ALLOC (int8_t, 4 * num_trees);

  if (alloc_vxyz)
    connectivity->vertices = P4EST_ALLOC (double, 3 * num_vertices);
  else
    connectivity->vertices = NULL;

  connectivity->vtt_offset = P4EST_ALLOC (p4est_topidx_t, num_vertices + 1);
  connectivity->vtt_offset[num_vertices] = -1;  /* catch bugs */

  connectivity->vertex_to_tree = P4EST_ALLOC (p4est_topidx_t, num_vtt);
  connectivity->vertex_to_vertex = P4EST_ALLOC (p4est_topidx_t, num_vtt);

  return connectivity;
}

void
p4est_connectivity_destroy (p4est_connectivity_t * connectivity)
{
  P4EST_FREE (connectivity->tree_to_vertex);
  P4EST_FREE (connectivity->tree_to_tree);
  P4EST_FREE (connectivity->tree_to_face);
  P4EST_FREE (connectivity->vertices);
  P4EST_FREE (connectivity->vtt_offset);
  P4EST_FREE (connectivity->vertex_to_tree);
  P4EST_FREE (connectivity->vertex_to_vertex);

  P4EST_FREE (connectivity);
}

bool
p4est_connectivity_is_valid (p4est_connectivity_t * connectivity)
{
  int                 num_found;
  int                 face, rface, nface, orientation, corner;
  p4est_topidx_t      tree, ntree, ctree;
  p4est_topidx_t      vertex, cvertex, corner_trees;
  p4est_topidx_t      v1, v2, w1, w2;
  const p4est_topidx_t num_trees = connectivity->num_trees;
  const p4est_topidx_t num_vertices = connectivity->num_vertices;
  const p4est_topidx_t num_vtt = connectivity->vtt_offset[num_vertices];
  const p4est_topidx_t *ttv = connectivity->tree_to_vertex;
  const p4est_topidx_t *ttt = connectivity->tree_to_tree;
  const int8_t       *ttf = connectivity->tree_to_face;
  const p4est_topidx_t *vtt = connectivity->vertex_to_tree;
  const p4est_topidx_t *vtv = connectivity->vertex_to_vertex;
  const p4est_topidx_t *voff = connectivity->vtt_offset;

  if (num_trees < 1 || num_vertices < 4) {
    P4EST_NOTICE ("Invalid numbers of trees or vertices");
    return false;
  }

  for (tree = 0; tree < num_trees; ++tree) {
    for (face = 0; face < 4; ++face) {
      ntree = ttt[tree * 4 + face];
      if (ntree < 0 || ntree >= num_trees) {
        P4EST_NOTICEF ("Tree range A in %lld %d\n", (long long) tree, face);
        return false;
      }
      rface = (int) ttf[tree * 4 + face];
      if (rface < 0 || rface >= 8) {
        P4EST_NOTICEF ("Face range in %lld %d\n", (long long) tree, face);
        return false;
      }
      nface = rface % 4;        /* clamp to a real face index */
      orientation = rface / 4;  /* 0 (same) or 1 (opposite) */
      if (ntree == tree) {
        /* no neighbor across this face or self-periodic */
        if (nface == face && orientation != 0) {
          P4EST_NOTICEF ("Face invalid in %lld %d\n", (long long) tree, face);
          return false;
        }
      }
      if (ntree != tree || nface != face) {
        /* check reciprocity */
        if (ttt[ntree * 4 + nface] != tree) {
          P4EST_NOTICEF ("Tree reciprocity in %lld %d\n",
                         (long long) tree, face);
          return false;
        }
        if ((int) ttf[ntree * 4 + nface] != face + 4 * orientation) {
          P4EST_NOTICEF ("Face reciprocity in %lld %d\n",
                         (long long) tree, face);
          return false;
        }

        /* a neighbor across this face */
        v1 = ttv[tree * 4 + face];
        v2 = ttv[tree * 4 + (face + 1) % 4];
        w1 = ttv[ntree * 4 + nface];
        w2 = ttv[ntree * 4 + (nface + 1) % 4];
        if (v1 == v2 || w1 == w2) {
          P4EST_NOTICEF ("Vertex invalid in %lld %d\n", (long long) tree,
                         face);
          return false;
        }
        if ((v1 == w2 && v2 == w1) && orientation != 0) {
          P4EST_NOTICEF ("Orientation mismatch A in %lld %d\n",
                         (long long) tree, face);
          return false;
        }
        if ((v1 == w1 && v2 == w2) && orientation != 1) {
          P4EST_NOTICEF ("Orientation mismatch B in %lld %d\n",
                         (long long) tree, face);
          return false;
        }
      }
    }
  }

  for (vertex = 0; vertex < num_vertices; ++vertex) {
    corner_trees = voff[vertex + 1] - voff[vertex];     /* same type */
    if (corner_trees <= 0 || voff[vertex + 1] > num_vtt) {
      P4EST_NOTICEF ("Vertex offset mismatch %lld\n", (long long) vertex);
      return false;
    }
    for (ctree = 0; ctree < corner_trees; ++ctree) {
      ntree = vtt[voff[vertex] + ctree];
      if (ntree < 0 || ntree >= num_trees) {
        P4EST_NOTICEF ("Tree range B in %lld %lld\n",
                       (long long) vertex, (long long) ntree);
        return false;
      }
      num_found = 0;
      for (corner = 0; corner < 4; ++corner) {
        if (ttv[4 * ntree + corner] == vertex) {
          num_found = 1;
          break;
        }
      }
      if (!num_found) {
        P4EST_NOTICEF ("Corner mismatch in %lld %lld\n",
                       (long long) vertex, (long long) ntree);
        return false;
      }
      cvertex = vtv[voff[vertex] + ctree];
      if (cvertex < 0 || cvertex >= num_vertices) {
        P4EST_NOTICEF ("Vertex mismatch in %lld %lld\n",
                       (long long) vertex, (long long) ntree);
        return false;
      }
    }
  }

  for (tree = 0; tree < num_trees; ++tree) {
    for (corner = 0; corner < 4; ++corner) {
      vertex = ttv[tree * 4 + corner];
      corner_trees = voff[vertex + 1] - voff[vertex];   /* same type */
      num_found = 0;
      for (ctree = 0; ctree < corner_trees; ++ctree) {
        if (vtt[voff[vertex] + ctree] == tree &&
            vtv[voff[vertex] + ctree] == vertex) {
          ++num_found;
        }
      }
      if (num_found != 1) {
        P4EST_NOTICEF ("Tree and vertex count in %lld %d\n",
                       (long long) tree, corner);
        return false;
      }
    }
  }

  return true;
}

bool
p4est_connectivity_is_equal (p4est_connectivity_t * conn1,
                             p4est_connectivity_t * conn2)
{
  size_t              topsize, int8size;
  size_t              tcount;
  p4est_topidx_t      num_vertices, num_vtt;

  topsize = sizeof (p4est_topidx_t);
  int8size = sizeof (int8_t);

  if (conn1->num_trees != conn2->num_trees ||
      conn1->num_vertices != conn2->num_vertices)
    return false;
  if ((conn1->vertices == NULL) != (conn2->vertices == NULL))
    return false;

  tcount = (size_t) (4 * conn1->num_trees);
  if (memcmp (conn1->tree_to_vertex, conn2->tree_to_vertex,
              tcount * topsize) ||
      memcmp (conn1->tree_to_tree, conn2->tree_to_tree,
              tcount * topsize) ||
      memcmp (conn1->tree_to_face, conn2->tree_to_face, tcount * int8size))
    return false;

  num_vertices = conn1->num_vertices;
  num_vtt = conn1->vtt_offset[num_vertices];

  if (conn1->vertices != NULL &&
      memcmp (conn1->vertices, conn2->vertices,
              sizeof (double) * 3 * num_vertices))
    return false;

  if (memcmp (conn1->vtt_offset, conn2->vtt_offset,
              topsize * (num_vertices + 1)) ||
      memcmp (conn1->vertex_to_tree, conn2->vertex_to_tree,
              topsize * num_vtt) ||
      memcmp (conn1->vertex_to_vertex, conn2->vertex_to_vertex,
              topsize * num_vtt))
    return false;

  return true;
}

static void
sc_fwrite (const void *ptr, size_t size, size_t nmemb, FILE * file,
           const char *errmsg)
{
  size_t              nwritten;

  nwritten = fwrite (ptr, size, nmemb, file);
  SC_CHECK_ABORT (nwritten == nmemb, errmsg);
}

static void
sc_fread (void *ptr, size_t size, size_t nmemb, FILE * file,
          const char *errmsg)
{
  size_t              nread;

  nread = fread (ptr, size, nmemb, file);
  SC_CHECK_ABORT (nread == nmemb, errmsg);
}

void
p4est_connectivity_save (p4est_connectivity_t * conn, const char *filename)
{
  int                 retval;
  bool                alloc_vxyz;
  size_t              u64z, topsize, int8size;
  size_t              tcount;
  uint64_t            array6[6];
  p4est_topidx_t      num_vertices, num_vtt;
  FILE               *file;

  P4EST_ASSERT (p4est_connectivity_is_valid (conn));

  file = fopen (filename, "wb");
  SC_CHECK_ABORT (file != NULL, "file open");

  alloc_vxyz = (conn->vertices != NULL);
  num_vertices = conn->num_vertices;
  num_vtt = conn->vtt_offset[num_vertices];

  u64z = sizeof (uint64_t);
  topsize = sizeof (p4est_topidx_t);
  int8size = sizeof (int8_t);
  array6[0] = P4EST_ONDISK_FORMAT;
  array6[1] = (uint64_t) topsize;
  array6[2] = (uint64_t) conn->num_trees;
  array6[3] = (uint64_t) num_vertices;
  array6[4] = (uint64_t) num_vtt;
  array6[5] = (uint64_t) alloc_vxyz;
  sc_fwrite (array6, u64z, 6, file, "write header");

  tcount = (size_t) (4 * conn->num_trees);
  sc_fwrite (conn->tree_to_vertex, topsize, tcount, file, "write ttv");
  sc_fwrite (conn->tree_to_tree, topsize, tcount, file, "write ttt");
  sc_fwrite (conn->tree_to_face, int8size, tcount, file, "write ttf");

  if (alloc_vxyz) {
    sc_fwrite (conn->vertices, sizeof (double), 3 * num_vertices, file,
               "write vertices");
  }

  sc_fwrite (conn->vtt_offset, topsize, num_vertices + 1, file,
             "write vtt_offset");
  sc_fwrite (conn->vertex_to_tree, topsize, num_vtt, file, "write vtt");
  sc_fwrite (conn->vertex_to_vertex, topsize, num_vtt, file, "write vtv");

  retval = fclose (file);
  SC_CHECK_ABORT (retval == 0, "file close");
}

p4est_connectivity_t *
p4est_connectivity_load (const char *filename)
{
  int                 retval;
  bool                alloc_vxyz;
  size_t              u64z, topsize, int8size;
  size_t              tcount;
  uint64_t            array6[6];
  p4est_topidx_t      num_trees, num_vertices, num_vtt;
  FILE               *file;
  p4est_connectivity_t *conn = NULL;

  file = fopen (filename, "rb");
  SC_CHECK_ABORT (file != NULL, "file open");

  u64z = sizeof (uint64_t);
  topsize = sizeof (p4est_topidx_t);
  int8size = sizeof (int8_t);
  sc_fread (array6, u64z, 6, file, "read header");

  SC_CHECK_ABORT (array6[0] == P4EST_ONDISK_FORMAT,
                  "on-disk format mismatch");
  SC_CHECK_ABORT (array6[1] == (uint64_t) topsize,
                  "p4est_topidx_t size mismatch");

  num_trees = (p4est_topidx_t) array6[2];
  num_vertices = (p4est_topidx_t) array6[3];
  num_vtt = (p4est_topidx_t) array6[4];
  alloc_vxyz = (bool) array6[5];
  SC_CHECK_ABORT (num_trees >= 0, "negative num_trees");
  SC_CHECK_ABORT (num_vertices >= 0, "negative num_vertices");
  SC_CHECK_ABORT (num_vtt >= 0, "negative num_vtt");

  conn =
    p4est_connectivity_new (num_trees, num_vertices, num_vtt, alloc_vxyz);

  tcount = (size_t) (4 * conn->num_trees);
  sc_fread (conn->tree_to_vertex, topsize, tcount, file, "read ttv");
  sc_fread (conn->tree_to_tree, topsize, tcount, file, "read ttt");
  sc_fread (conn->tree_to_face, int8size, tcount, file, "read ttf");

  if (alloc_vxyz) {
    sc_fread (conn->vertices, sizeof (double), 3 * num_vertices, file,
              "read vertices");
  }

  sc_fread (conn->vtt_offset, topsize, num_vertices + 1, file,
            "read vtt_offset");
  SC_CHECK_ABORT (num_vtt == conn->vtt_offset[num_vertices],
                  "num_vtt mismatch");
  sc_fread (conn->vertex_to_tree, topsize, num_vtt, file, "read vtt");
  sc_fread (conn->vertex_to_vertex, topsize, num_vtt, file, "read vtv");

  retval = fclose (file);
  SC_CHECK_ABORT (retval == 0, "file close");

  SC_CHECK_ABORT (p4est_connectivity_is_valid (conn), "invalid connectivity");

  return conn;
}

p4est_connectivity_t *
p4est_connectivity_new_unitsquare (void)
{
  p4est_connectivity_t *connectivity;

  connectivity = p4est_connectivity_new (1, 4, 4, true);

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
  const p4est_topidx_t num_trees = 3;
  const p4est_topidx_t num_vertices = 7;
  const p4est_topidx_t num_vtt = 12;
  const p4est_topidx_t tree_to_vertex[3 * 4] = {
    0, 1, 3, 2, 0, 2, 5, 6, 2, 3, 4, 5,
  };
  const p4est_topidx_t tree_to_tree[3 * 4] = {
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
  const p4est_topidx_t vtt_offset[7 + 1] = {
    0, 2, 3, 6, 8, 9, 11, 12,
  };
  const p4est_topidx_t vertex_to_tree[12] = {
    0, 1, 0, 0, 2, 1, 0, 2, 2, 1, 2, 1,
  };
  const p4est_topidx_t vertex_to_vertex[12] = {
    0, 0, 1, 2, 2, 2, 3, 3, 4, 5, 5, 6,
  };

  connectivity =
    p4est_connectivity_new (num_trees, num_vertices, num_vtt, true);

  memcpy (connectivity->tree_to_vertex, tree_to_vertex,
          sizeof (p4est_topidx_t) * 4 * num_trees);
  memcpy (connectivity->tree_to_tree, tree_to_tree,
          sizeof (p4est_topidx_t) * 4 * num_trees);
  memcpy (connectivity->tree_to_face, tree_to_face,
          sizeof (int8_t) * 4 * num_trees);
  memcpy (connectivity->vertices, vertices,
          sizeof (double) * 3 * num_vertices);
  memcpy (connectivity->vtt_offset, vtt_offset,
          sizeof (p4est_topidx_t) * (num_vertices + 1));
  memcpy (connectivity->vertex_to_tree, vertex_to_tree,
          sizeof (p4est_topidx_t) * num_vtt);
  memcpy (connectivity->vertex_to_vertex, vertex_to_vertex,
          sizeof (p4est_topidx_t) * num_vtt);

  P4EST_ASSERT (p4est_connectivity_is_valid (connectivity));

  return connectivity;
}

p4est_connectivity_t *
p4est_connectivity_new_moebius (void)
{
  p4est_connectivity_t *connectivity;
  const p4est_topidx_t num_trees = 5;
  const p4est_topidx_t num_vertices = 10;
  const p4est_topidx_t num_vtt = 20;
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
  const p4est_topidx_t tree_to_vertex[5 * 4] = {
    0, 2, 3, 1,
    3, 5, 4, 2,                 /* left-handed */
    4, 6, 7, 5,
    6, 7, 8, 9,
    9, 8, 0, 1,
  };
  const p4est_topidx_t tree_to_tree[5 * 4] = {
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
  const p4est_topidx_t vtt_offset[10 + 1] = {
    0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20,
  };
  const p4est_topidx_t vertex_to_tree[20] = {
    0, 4, 0, 4,
    0, 1, 0, 1,
    1, 2, 1, 2,
    2, 3, 2, 3,
    3, 4, 3, 4,
  };
  const p4est_topidx_t vertex_to_vertex[20] = {
    0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9,
  };

  connectivity =
    p4est_connectivity_new (num_trees, num_vertices, num_vtt, true);

  memcpy (connectivity->tree_to_vertex, tree_to_vertex,
          sizeof (p4est_topidx_t) * 4 * num_trees);
  memcpy (connectivity->tree_to_tree, tree_to_tree,
          sizeof (p4est_topidx_t) * 4 * num_trees);
  memcpy (connectivity->tree_to_face, tree_to_face,
          sizeof (int8_t) * 4 * num_trees);
  memcpy (connectivity->vertices, vertices,
          sizeof (double) * 3 * num_vertices);
  memcpy (connectivity->vtt_offset, vtt_offset,
          sizeof (p4est_topidx_t) * (num_vertices + 1));
  memcpy (connectivity->vertex_to_tree, vertex_to_tree,
          sizeof (p4est_topidx_t) * num_vtt);
  memcpy (connectivity->vertex_to_vertex, vertex_to_vertex,
          sizeof (p4est_topidx_t) * num_vtt);

  P4EST_ASSERT (p4est_connectivity_is_valid (connectivity));

  return connectivity;
}

p4est_connectivity_t *
p4est_connectivity_new_star (void)
{
  int                 i;
  p4est_connectivity_t *connectivity;
  const double        pi = 4.0 * atan (1.0);
  const double        r1 = 1.;
  const double        r2 = 1.5;
  const p4est_topidx_t num_trees = 6;
  const p4est_topidx_t num_vertices = 13;
  const p4est_topidx_t num_vtt = 24;
  const p4est_topidx_t tree_to_vertex[6 * 4] = {
    0, 1, 2, 3,
    0, 3, 4, 5,
    5, 6, 7, 0,
    8, 7, 0, 9,                 /* left-handed */
    9, 0, 11, 10,               /* left-handed */
    12, 1, 0, 11,
  };
  const p4est_topidx_t tree_to_tree[6 * 4] = {
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
  const p4est_topidx_t vtt_offset[13 + 1] = {
    0, 6, 8, 9, 11, 12, 14, 15, 17, 18, 20, 21, 23, 24,
  };
  const p4est_topidx_t vertex_to_tree[24] = {
    0, 1, 2, 3, 4, 5, 5, 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5,
  };
  const p4est_topidx_t vertex_to_vertex[24] = {
    0, 0, 0, 0, 0, 0, 1, 1, 2, 3, 3, 4,
    5, 5, 6, 7, 7, 8, 9, 9, 10, 11, 11, 12,
  };

  connectivity =
    p4est_connectivity_new (num_trees, num_vertices, num_vtt, true);

  connectivity->vertices[0 * 3 + 0] = 0;
  connectivity->vertices[0 * 3 + 1] = 0;
  connectivity->vertices[0 * 3 + 2] = 0;
  for (i = 0; i < 6; ++i) {
    connectivity->vertices[(2 * i + 1) * 3 + 0] = r1 * cos (i * pi / 3);
    connectivity->vertices[(2 * i + 1) * 3 + 1] = r1 * sin (i * pi / 3);
    connectivity->vertices[(2 * i + 1) * 3 + 2] = 0;
    connectivity->vertices[(2 * i + 2) * 3 + 0] =
      r2 * cos ((i + .5) * pi / 3);
    connectivity->vertices[(2 * i + 2) * 3 + 1] =
      r2 * sin ((i + .5) * pi / 3);
    connectivity->vertices[(2 * i + 2) * 3 + 2] = 0;
  }

  memcpy (connectivity->tree_to_vertex, tree_to_vertex,
          sizeof (p4est_topidx_t) * 4 * num_trees);
  memcpy (connectivity->tree_to_tree, tree_to_tree,
          sizeof (p4est_topidx_t) * 4 * num_trees);
  memcpy (connectivity->tree_to_face, tree_to_face,
          sizeof (int8_t) * 4 * num_trees);
  memcpy (connectivity->vtt_offset, vtt_offset,
          sizeof (p4est_topidx_t) * (num_vertices + 1));
  memcpy (connectivity->vertex_to_tree, vertex_to_tree,
          sizeof (p4est_topidx_t) * num_vtt);
  memcpy (connectivity->vertex_to_vertex, vertex_to_vertex,
          sizeof (p4est_topidx_t) * num_vtt);

  P4EST_ASSERT (p4est_connectivity_is_valid (connectivity));

  return connectivity;
}

p4est_connectivity_t *
p4est_connectivity_new_periodic (void)
{
  p4est_connectivity_t *connectivity;
  const p4est_topidx_t num_trees = 1;
  const p4est_topidx_t num_vertices = 4;
  const p4est_topidx_t num_vtt = 16;
  const p4est_topidx_t tree_to_vertex[1 * 4] = {
    0, 1, 2, 3,
  };
  const p4est_topidx_t tree_to_tree[1 * 4] = {
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
  const p4est_topidx_t vtt_offset[4 + 1] = {
    0, 4, 8, 12, 16,
  };
  const p4est_topidx_t vertex_to_tree[16] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  };
  const p4est_topidx_t vertex_to_vertex[16] = {
    0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
  };

  connectivity =
    p4est_connectivity_new (num_trees, num_vertices, num_vtt, true);

  memcpy (connectivity->tree_to_vertex, tree_to_vertex,
          sizeof (p4est_topidx_t) * 4 * num_trees);
  memcpy (connectivity->tree_to_tree, tree_to_tree,
          sizeof (p4est_topidx_t) * 4 * num_trees);
  memcpy (connectivity->tree_to_face, tree_to_face,
          sizeof (int8_t) * 4 * num_trees);
  memcpy (connectivity->vertices, vertices,
          sizeof (double) * 3 * num_vertices);
  memcpy (connectivity->vtt_offset, vtt_offset,
          sizeof (p4est_topidx_t) * (num_vertices + 1));
  memcpy (connectivity->vertex_to_tree, vertex_to_tree,
          sizeof (p4est_topidx_t) * num_vtt);
  memcpy (connectivity->vertex_to_vertex, vertex_to_vertex,
          sizeof (p4est_topidx_t) * num_vtt);

  P4EST_ASSERT (p4est_connectivity_is_valid (connectivity));

  return connectivity;
}

int
p4est_find_face_transform (p4est_connectivity_t * connectivity,
                           p4est_topidx_t itree, int iface)
{
  int                 nrface, neighbor_face, orientation;
  p4est_topidx_t      neighbor_tree;

  P4EST_ASSERT (itree >= 0 && itree < connectivity->num_trees);
  P4EST_ASSERT (iface >= 0 && iface < 4);

  neighbor_tree = connectivity->tree_to_tree[4 * itree + iface];
  P4EST_ASSERT (neighbor_tree >= 0 &&
                neighbor_tree < connectivity->num_trees);

  nrface = (int) connectivity->tree_to_face[4 * itree + iface];
  P4EST_ASSERT (nrface >= 0 && nrface < 8);
  neighbor_face = nrface % 4;
  orientation = nrface / 4;

  if (neighbor_tree == itree && neighbor_face == iface) {
    return -1;
  }

  return p4est_transform_table[iface][neighbor_face][orientation];
}

void
p4est_find_corner_transform (p4est_connectivity_t * conn,
                             p4est_topidx_t itree, int icorner,
                             sc_array_t * ctransforms)
{
  int                 ncorner;
  p4est_topidx_t      corner_trees, ctree;
  p4est_topidx_t      ntree, ntree1, ntree2;
  p4est_topidx_t      ivertex, nvertex;
  p4est_corner_transform_t *ct;

  P4EST_ASSERT (0 <= itree && itree < conn->num_trees);
  P4EST_ASSERT (0 <= icorner && icorner < 4);
  P4EST_ASSERT (ctransforms->elem_size == sizeof (p4est_corner_transform_t));
  sc_array_resize (ctransforms, 0);

  ivertex = conn->tree_to_vertex[4 * itree + icorner];
  P4EST_ASSERT (0 <= ivertex && ivertex < conn->num_vertices);

  ntree1 = conn->tree_to_tree[4 * itree + (icorner + 3) % 4];
  ntree2 = conn->tree_to_tree[4 * itree + icorner];

  corner_trees =                /* same type */
    conn->vtt_offset[ivertex + 1] - conn->vtt_offset[ivertex];

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
    ct = sc_array_push (ctransforms);
    ct->ntree = ntree;
    ct->ncorner = (int8_t) p4est_corner_to_zorder[ncorner];
  }
  P4EST_ASSERT (corner_trees == (p4est_topidx_t) ctransforms->elem_count
                + 1 + (ntree1 != itree) + (ntree2 != itree));
}

/* EOF p4est_connectivity.c */
