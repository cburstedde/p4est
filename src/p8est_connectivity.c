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

/*
 * If num_vertices == 0 or ttv == NULL or vertices == NULL,
 * then there will be no vertex information created.
 */
static p8est_connectivity_t *
p8est_connectivity_new_copy (p4est_topidx_t num_vertices,
                             p4est_topidx_t num_trees,
                             p4est_topidx_t num_edges,
                             p4est_topidx_t num_corners,
                             const double *vertices,
                             const p4est_topidx_t * ttv,
                             const p4est_topidx_t * ttt,
                             const int8_t * ttf,
                             const p4est_topidx_t * tte,
                             const p4est_topidx_t * eoff,
                             const p4est_topidx_t * ett,
                             const int8_t * ete,
                             const p4est_topidx_t * ttc,
                             const p4est_topidx_t * coff,
                             const p4est_topidx_t * ctt, const int8_t * ctc)
{
  const bool          alloc_vertices = (vertices != NULL);
  p4est_topidx_t      num_ett;
  p4est_topidx_t      num_ctt;
  p8est_connectivity_t *conn;

  num_ett = eoff[num_edges];
  num_ctt = coff[num_corners];
  conn = p8est_connectivity_new (num_vertices, num_trees,
                                 num_edges, num_ett, num_corners, num_ctt,
                                 alloc_vertices);

  if (alloc_vertices) {
    memcpy (conn->vertices, vertices, sizeof (double) * 3 * num_vertices);
  }
  else {
    conn->vertices = NULL;
  }
  memcpy (conn->tree_to_vertex, ttv, sizeof (p4est_topidx_t) * 8 * num_trees);
  memcpy (conn->tree_to_tree, ttt, sizeof (p4est_topidx_t) * 6 * num_trees);
  memcpy (conn->tree_to_face, ttf, sizeof (int8_t) * 6 * num_trees);

  if (num_edges > 0) {
    memcpy (conn->tree_to_edge, tte,
            sizeof (p4est_topidx_t) * 12 * num_trees);
    memcpy (conn->edge_to_tree, ett, sizeof (p4est_topidx_t) * num_ett);
    memcpy (conn->edge_to_edge, ete, sizeof (int8_t) * num_ett);
  }
  memcpy (conn->ett_offset, eoff, sizeof (p4est_topidx_t) * (num_edges + 1));

  if (num_corners > 0) {
    memcpy (conn->tree_to_corner, ttc,
            sizeof (p4est_topidx_t) * 8 * num_trees);
    memcpy (conn->corner_to_tree, ctt, sizeof (p4est_topidx_t) * num_ctt);
    memcpy (conn->corner_to_corner, ctc, sizeof (int8_t) * num_ctt);
  }
  memcpy (conn->ctt_offset, coff,
          sizeof (p4est_topidx_t) * (num_corners + 1));

  P4EST_ASSERT (p8est_connectivity_is_valid (conn));

  return conn;
}

p8est_connectivity_t *
p8est_connectivity_new (p4est_topidx_t num_vertices,
                        p4est_topidx_t num_trees,
                        p4est_topidx_t num_edges,
                        p4est_topidx_t num_ett,
                        p4est_topidx_t num_corners, p4est_topidx_t num_ctt,
                        bool alloc_vertices)
{
  p8est_connectivity_t *conn;

  conn = P4EST_ALLOC_ZERO (p8est_connectivity_t, 1);

  conn->num_vertices = num_vertices;
  conn->num_trees = num_trees;
  if (alloc_vertices) {
    conn->vertices = P4EST_ALLOC (double, 3 * num_vertices);
  }
  else {
    conn->vertices = NULL;
  }
  conn->tree_to_vertex = P4EST_ALLOC (p4est_topidx_t, 8 * num_trees);
  conn->tree_to_tree = P4EST_ALLOC (p4est_topidx_t, 6 * num_trees);
  conn->tree_to_face = P4EST_ALLOC (int8_t, 6 * num_trees);

  conn->num_edges = num_edges;
  if (num_edges > 0) {
    conn->tree_to_edge = P4EST_ALLOC (p4est_topidx_t, 12 * num_trees);
    conn->edge_to_tree = P4EST_ALLOC (p4est_topidx_t, num_ett);
    conn->edge_to_edge = P4EST_ALLOC (int8_t, num_ett);
  }
  else {
    conn->tree_to_edge = NULL;
    conn->edge_to_tree = NULL;
    conn->edge_to_edge = NULL;
  }
  conn->ett_offset = P4EST_ALLOC (p4est_topidx_t, num_edges + 1);
  conn->ett_offset[num_edges] = num_ett;

  conn->num_corners = num_corners;
  if (num_corners > 0) {
    conn->tree_to_corner = P4EST_ALLOC (p4est_topidx_t, 8 * num_trees);
    conn->corner_to_tree = P4EST_ALLOC (p4est_topidx_t, num_ctt);
    conn->corner_to_corner = P4EST_ALLOC (int8_t, num_ctt);
  }
  else {
    conn->tree_to_corner = NULL;
    conn->corner_to_tree = NULL;
    conn->corner_to_corner = NULL;
  }
  conn->ctt_offset = P4EST_ALLOC (p4est_topidx_t, num_corners + 1);
  conn->ctt_offset[num_corners] = num_ctt;

  return conn;
}

void
p8est_connectivity_destroy (p8est_connectivity_t * conn)
{
  P4EST_FREE (conn->vertices);
  P4EST_FREE (conn->tree_to_vertex);

  P4EST_FREE (conn->tree_to_tree);
  P4EST_FREE (conn->tree_to_face);

  P4EST_FREE (conn->tree_to_edge);
  P4EST_FREE (conn->ett_offset);
  P4EST_FREE (conn->edge_to_tree);
  P4EST_FREE (conn->edge_to_edge);

  P4EST_FREE (conn->tree_to_corner);
  P4EST_FREE (conn->ctt_offset);
  P4EST_FREE (conn->corner_to_tree);
  P4EST_FREE (conn->corner_to_corner);

  P4EST_FREE (conn);
}

bool
p8est_connectivity_is_valid (p8est_connectivity_t * conn)
{
  int                 fvert, face_ref, face_perm, nvert;
  int                 face, rface, nface, orientation;
  int                 edge, nedge, corner, ncorner;
  int                 flip, nflip, nflip1, nflip2;
  int                 ecode, ecount;
  p4est_topidx_t      v0, v1, nv0, nv1, ntree1, ntree2;
  p4est_topidx_t      vertex, nvertex, tree, ntree;
  p4est_topidx_t      aedge, edge_begin, edge_end;
  p4est_topidx_t      acorner, corner_begin, corner_end;
  p4est_topidx_t      nett, nctt;
  const p4est_topidx_t num_vertices = conn->num_vertices;
  const p4est_topidx_t num_trees = conn->num_trees;
  const p4est_topidx_t num_edges = conn->num_edges;
  const p4est_topidx_t num_corners = conn->num_corners;
  const p4est_topidx_t *ttv = conn->tree_to_vertex;
  const p4est_topidx_t *ttt = conn->tree_to_tree;
  const int8_t       *ttf = conn->tree_to_face;
  const p4est_topidx_t *tte = conn->tree_to_edge;
  const p4est_topidx_t *eoff = conn->ett_offset;
  const p4est_topidx_t *ett = conn->edge_to_tree;
  const int8_t       *ete = conn->edge_to_edge;
  const p4est_topidx_t num_ett = eoff[num_edges];
  const p4est_topidx_t *ttc = conn->tree_to_corner;
  const p4est_topidx_t *coff = conn->ctt_offset;
  const p4est_topidx_t *ctt = conn->corner_to_tree;
  const int8_t       *ctc = conn->corner_to_corner;
  const p4est_topidx_t num_ctt = coff[num_corners];

  for (nett = 0; nett < num_ett; ++nett) {
    if (ett[nett] < 0 || ett[nett] >= num_trees) {
      P4EST_NOTICEF ("Edge to tree %lld out of range\n", (long long) nett);
      return false;
    }
    if (ete[nett] < 0 || ete[nett] >= 24) {
      P4EST_NOTICEF ("Edge to edge %lld out of range\n", (long long) nett);
      return false;
    }
  }

  for (nctt = 0; nctt < num_ctt; ++nctt) {
    if (ctt[nctt] < 0 || ctt[nctt] >= num_trees) {
      P4EST_NOTICEF ("Corner to tree %lld out of range\n", (long long) nctt);
      return false;
    }
    if (ctc[nctt] < 0 || ctc[nctt] >= 8) {
      P4EST_NOTICEF ("Corner to corner %lld out of range\n",
                     (long long) nctt);
      return false;
    }
  }

  for (tree = 0; tree < num_trees; ++tree) {
    for (nvert = 0; nvert < 8; ++nvert) {
      vertex = ttv[tree * 8 + nvert];
      if (vertex < 0 || vertex >= num_vertices) {
        P4EST_NOTICEF ("Tree to vertex out of range %lld %d",
                       (long long) tree, nvert);
        return false;
      }
    }
  }

  for (tree = 0; tree < num_trees; ++tree) {
    for (face = 0; face < 6; ++face) {
      ntree = ttt[tree * 6 + face];
      if (ntree < 0 || ntree >= num_trees) {
        P4EST_NOTICEF ("Tree to tree out of range %lld %d\n",
                       (long long) tree, face);
        return false;
      }
      rface = (int) ttf[tree * 6 + face];
      if (rface < 0 || rface >= 24) {
        P4EST_NOTICEF ("Tree to face out of range %lld %d\n",
                       (long long) tree, face);
        return false;
      }
      nface = rface % 6;        /* clamp to a real face index */
      orientation = rface / 6;  /* 0..3 for relative rotation */
      if (ntree == tree) {
        /* no neighbor across this face or self-periodic */
        if (nface == face && orientation != 0) {
          P4EST_NOTICEF ("Face invalid in %lld %d\n", (long long) tree, face);
          return false;
        }
      }
      if (ntree != tree || nface != face) {
        /* check reciprocity */
        if (ttt[ntree * 6 + nface] != tree) {
          P4EST_NOTICEF ("Tree to tree reciprocity in %lld %d\n",
                         (long long) tree, face);
          return false;
        }
        if ((int) ttf[ntree * 6 + nface] != face + 6 * orientation) {
          P4EST_NOTICEF ("Tree to face reciprocity in %lld %d\n",
                         (long long) tree, face);
          return false;
        }

        /* check vertex consistency across faces */
        for (fvert = 0; fvert < 4; ++fvert) {
          vertex = ttv[tree * 8 + p8est_face_vertices[face][fvert]];
          if (vertex < 0 || vertex >= num_vertices) {
            P4EST_NOTICEF ("Invalid vertex in %lld %d %d\n",
                           (long long) tree, face, fvert);
            return false;
          }
          for (nvert = fvert + 1; nvert < 4; ++nvert) {
            nvertex = ttv[tree * 8 + p8est_face_vertices[face][nvert]];
            if (vertex == nvertex) {
              P4EST_NOTICEF ("Duplicate vertex in %lld %d %d %d",
                             (long long) tree, face, fvert, nvert);
              return false;
            }
          }
          face_ref = p8est_face_permutation_refs[face][nface];
          face_perm = p8est_face_permutation_sets[face_ref][orientation];
          nvert = p8est_face_permutations[face_perm][fvert];
          nvertex = ttv[ntree * 8 + p8est_face_vertices[nface][nvert]];
          if (ntree != tree && nvertex != vertex) {
            P4EST_NOTICEF ("Vertex reciprocity in %lld %d %d\n",
                           (long long) tree, face, fvert);
            return false;
          }
        }
      }
    }

    if (num_edges > 0) {
      for (edge = 0; edge < 12; ++edge) {
        aedge = tte[tree * 12 + edge];
        if (aedge < -1 || aedge >= num_edges) {
          P4EST_NOTICEF ("Tree to edge out of range %lld %d\n",
                         (long long) tree, edge);
          return false;
        }
        if (aedge == -1) {
          continue;
        }
        ecode = ecount = 0;
        flip = nflip1 = nflip2 = -1;
        v0 = ttv[8 * tree + p8est_edge_vertices[edge][0]];
        v1 = ttv[8 * tree + p8est_edge_vertices[edge][1]];
        ntree1 = ttt[6 * tree + ((edge < 4) ? 2 : 0) + edge % 2];
        ntree2 = ttt[6 * tree + ((edge < 8) ? 4 : 2) + (edge / 2) % 2];
        edge_begin = eoff[aedge];
        edge_end = eoff[aedge + 1];
        if (edge_begin < 0 || edge_begin >= num_ett ||
            edge_end < 0 || edge_end > num_ett) {
          P4EST_NOTICEF ("Invalid edge range %lld %d\n",
                         (long long) tree, edge);
          return false;
        }
        for (nett = edge_begin; nett < edge_end; ++nett) {
          ntree = ett[nett];
          nedge = (int) ete[nett] % 12;
          if (tte[ntree * 12 + nedge] != aedge) {
            P4EST_NOTICEF ("Edge to edge reciprocity in %lld %d %lld\n",
                           (long long) tree, edge, (long long) nett);
            return false;
          }
          nflip = (int) ete[nett] / 12;
          if (ntree == tree && nedge == edge) {
            if (flip != -1) {
              ecode = 1;
              break;
            }
            flip = nflip;
            continue;
          }
          if (ntree == ntree1 || ntree == ntree2) {
            nv0 = ttv[8 * ntree + p8est_edge_vertices[nedge][0]];
            nv1 = ttv[8 * ntree + p8est_edge_vertices[nedge][1]];
            if (nv0 == v0 && nv1 == v1) {
              if (ntree == ntree1) {
                if (nflip1 != -1) {
                  ecode = 2;
                  break;
                }
                nflip1 = nflip;
              }
              else {
                if (nflip2 != -1) {
                  ecode = 3;
                  break;
                }
                nflip2 = nflip;
              }
              continue;
            }
            else if (nv0 == v1 && nv1 == v0) {
              if (ntree == ntree1) {
                if (nflip1 != -1) {
                  ecode = 4;
                  break;
                }
                nflip1 = 1 - nflip;
              }
              else {
                if (nflip2 != -1) {
                  ecode = 5;
                  break;
                }
                nflip2 = 1 - nflip;
              }
              continue;
            }
          }
          ++ecount;
        }
        if (ecode > 0) {
          P4EST_NOTICEF ("Shared edge %lld %d %lld inconsistency %d\n",
                         (long long) tree, edge, (long long) nett, ecode);
          return false;
        }
        if ((int) (edge_end - edge_begin) !=
            ecount + 1 + (ntree1 != tree) + (ntree2 != tree)) {
          P4EST_NOTICEF ("Shared edge %lld %d inconsistent count\n",
                         (long long) tree, edge);
          return false;
        }
        if (flip == -1 ||
            !(nflip1 == -1 || nflip1 == flip) ||
            !(nflip2 == -1 || nflip2 == flip)) {
          P4EST_NOTICEF ("Shared edge %lld %d inconsistent flip\n",
                         (long long) tree, edge);
          return false;
        }
      }
    }

    if (num_corners > 0) {
      for (corner = 0; corner < 8; ++corner) {
        acorner = ttc[tree * 8 + corner];
        if (acorner < -1 || acorner >= num_corners) {
          P4EST_NOTICEF ("Tree to corner out of range %lld %d\n",
                         (long long) tree, corner);
          return false;
        }
        if (acorner == -1) {
          continue;
        }
        corner_begin = eoff[acorner];
        corner_end = eoff[acorner + 1];
        if (corner_begin < 0 || corner_begin >= num_ctt ||
            corner_end < 0 || corner_end > num_ctt) {
          P4EST_NOTICEF ("Invalid corner range %lld %d\n",
                         (long long) tree, corner);
          return false;
        }
        for (nctt = corner_begin; nctt < corner_end; ++nctt) {
          ntree = ctt[nctt];
          ncorner = (int) ctc[nctt];
          if (ttc[ntree * 8 + ncorner] != acorner) {
            P4EST_NOTICEF ("Corner to corner reciprocity in %lld %d %lld\n",
                           (long long) tree, corner, (long long) nctt);
          }
        }
      }
    }
  }

  return true;
}

p8est_connectivity_t *
p8est_connectivity_new_unitcube (void)
{
  const p4est_topidx_t num_vertices = 8;
  const p4est_topidx_t num_trees = 1;
  const p4est_topidx_t num_ett = 0;
  const p4est_topidx_t num_ctt = 0;
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
  const p4est_topidx_t tree_to_vertex[1 * 8] = {
    0, 1, 2, 3, 4, 5, 6, 7,
  };
  const p4est_topidx_t tree_to_tree[1 * 6] = {
    0, 0, 0, 0, 0, 0,
  };
  const int8_t        tree_to_face[1 * 6] = {
    0, 1, 2, 3, 4, 5,
  };

  return p8est_connectivity_new_copy (num_vertices, num_trees, 0, 0,
                                      vertices, tree_to_vertex,
                                      tree_to_tree, tree_to_face,
                                      NULL, &num_ett, NULL, NULL,
                                      NULL, &num_ctt, NULL, NULL);
}

p8est_connectivity_t *
p8est_connectivity_new_periodic (void)
{
  const p4est_topidx_t num_vertices = 8;
  const p4est_topidx_t num_trees = 1;
  const p4est_topidx_t num_edges = 4;
  const p4est_topidx_t num_corners = 1;
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
  const p4est_topidx_t tree_to_vertex[1 * 8] = {
    0, 1, 2, 3, 4, 5, 6, 7,
  };
  const p4est_topidx_t tree_to_tree[1 * 6] = {
    0, 0, 0, 0, 0, 0,
  };
  const int8_t        tree_to_face[1 * 6] = {
    1, 0, 2, 3, 11, 10,
  };
  const p4est_topidx_t tree_to_edge[1 * 12] = {
    0, 0, 1, 1, 1, 1, 0, 0, 2, 2, 3, 3,
  };
  const p4est_topidx_t ett_offset[4 + 1] = {
    0, 4, 8, 10, 12,
  };
  const p4est_topidx_t edge_to_tree[12] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  };
  const int8_t        edge_to_edge[12] = {
    0, 7, 1, 6, 2, 4, 3, 5, 8, 9, 10, 11,
  };
  const p4est_topidx_t tree_to_corner[1 * 8] = {
    0, 0, 0, 0, 0, 0, 0, 0,
  };
  const p4est_topidx_t ctt_offset[1 + 1] = {
    0, 8,
  };
  const p4est_topidx_t corner_to_tree[8] = {
    0, 0, 0, 0, 0, 0, 0, 0,
  };
  const int8_t        corner_to_corner[8] = {
    0, 1, 2, 3, 4, 5, 6, 7,
  };

  return p8est_connectivity_new_copy (num_vertices, num_trees,
                                      num_edges, num_corners,
                                      vertices, tree_to_vertex,
                                      tree_to_tree, tree_to_face,
                                      tree_to_edge, ett_offset,
                                      edge_to_tree, edge_to_edge,
                                      tree_to_corner, ctt_offset,
                                      corner_to_tree, corner_to_corner);
}

p8est_connectivity_t *
p8est_connectivity_new_twocubes (void)
{
  const p4est_topidx_t num_vertices = 12;
  const p4est_topidx_t num_trees = 2;
  const p4est_topidx_t num_ett = 0;
  const p4est_topidx_t num_ctt = 0;
  const double        vertices[3 * 12] = {
    0, 0, 0,
    1, 0, 0,
    2, 0, 0,
    0, 1, 0,
    1, 1, 0,
    2, 1, 0,
    0, 0, 1,
    1, 0, 1,
    2, 0, 1,
    0, 1, 1,
    1, 1, 1,
    2, 1, 1,
  };
  const p4est_topidx_t tree_to_vertex[2 * 8] = {
    0, 1, 3, 4, 6, 7, 9, 10,
    1, 2, 4, 5, 7, 8, 10, 11,
  };
  const p4est_topidx_t tree_to_tree[2 * 6] = {
    0, 1, 0, 0, 0, 0,
    0, 1, 1, 1, 1, 1,
  };
  const int8_t        tree_to_face[2 * 6] = {
    0, 0, 2, 3, 4, 5,
    1, 1, 2, 3, 4, 5,
  };

  return p8est_connectivity_new_copy (num_vertices, num_trees, 0, 0,
                                      vertices, tree_to_vertex,
                                      tree_to_tree, tree_to_face,
                                      NULL, &num_ett, NULL, NULL,
                                      NULL, &num_ctt, NULL, NULL);
}

p8est_connectivity_t *
p8est_connectivity_new_rotcubes (void)
{
  const p4est_topidx_t num_vertices = 26;
  const p4est_topidx_t num_trees = 6;
  const p4est_topidx_t num_edges = 3;
  const p4est_topidx_t num_corners = 0;
  const double        vertices[3 * 26] = {
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
    2.5, 1.5, 2,
    2, 1.5, 2,
    2, 1.5, 2.5,
    2, .5, 2.5,
    2.5, .5, 2,
    2, .5, 2,
  };
  const p4est_topidx_t tree_to_vertex[6 * 8] = {
    0, 17, 3, 4, 15, 11, 13, 14,
    7, 2, 6, 17, 9, 12, 8, 11,
    2, 12, 5, 10, 17, 11, 4, 14,
    19, 13, 18, 14, 16, 15, 1, 11,
    14, 11, 21, 25, 18, 1, 22, 23,
    21, 20, 25, 24, 14, 10, 11, 12,
  };
  const p4est_topidx_t tree_to_tree[6 * 6] = {
    0, 2, 0, 0, 0, 3,
    1, 2, 1, 1, 1, 1,
    2, 5, 1, 2, 2, 0,
    3, 0, 3, 4, 3, 3,
    4, 4, 3, 4, 5, 4,
    4, 5, 5, 5, 5, 2,
  };
  const int8_t        tree_to_face[6 * 6] = {
    0, 5, 2, 3, 4, 13,
    0, 2, 2, 3, 4, 5,
    0, 23, 1, 3, 4, 1,
    0, 17, 2, 8, 4, 5,
    0, 1, 9, 3, 12, 5,
    16, 1, 2, 3, 4, 19,
  };
  const p4est_topidx_t tree_to_edge[6 * 12] = {
    -1, -1, -1, -1, -1, -1, -1, 0, -1, 2, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, 2,
    -1, -1, 2, -1, -1, -1, -1, 0, -1, 1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0,
    0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, 1, -1, -1, 0, -1, -1, -1, -1, -1,
  };
  const p4est_topidx_t ett_offset[3 + 1] = { 0, 5, 8, 11 };
  const p4est_topidx_t edge_to_tree[11] = {
    0, 2, 3, 4, 5, 1, 2, 5, 0, 1, 2
  };
  const int8_t        edge_to_edge[11] = {
    7, 7, 23, 12, 18, 7, 9, 15, 9, 11, 2
  };
  const p4est_topidx_t tree_to_corner[6 * 8] = {
    -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1,
  };
  const p4est_topidx_t ctt_offset = 0;

  return p8est_connectivity_new_copy (num_vertices, num_trees,
                                      num_edges, num_corners,
                                      vertices, tree_to_vertex,
                                      tree_to_tree, tree_to_face,
                                      tree_to_edge, ett_offset,
                                      edge_to_tree, edge_to_edge,
                                      tree_to_corner, &ctt_offset,
                                      NULL, NULL);
}

p4est_topidx_t
p8est_find_face_transform (p8est_connectivity_t * connectivity,
                           p4est_topidx_t my_tree, int my_face,
                           int ftransform[])
{
  int                 i;
  int                 target_code, target_face, orientation;
  int                 face_ref, face_perm;
  int                 low[2], high[2], swap;
  int                *my_axis = &ftransform[0];
  int                *target_axis = &ftransform[3];
  int                *edge_reverse = &ftransform[6];
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
  target_axis[0] = target_axis[1] = -1;
  edge_reverse[0] = edge_reverse[1] = 0;

  /* find matching target vertices. TODO: precompute this */
  face_ref = p8est_face_permutation_refs[my_face][target_face];
  face_perm = p8est_face_permutation_sets[face_ref][orientation];
  low[0] = low[1] =
    p8est_face_vertices[target_face][p8est_face_permutations[face_perm][0]];
  high[0] =
    p8est_face_vertices[target_face][p8est_face_permutations[face_perm][1]];
  high[1] =
    p8est_face_vertices[target_face][p8est_face_permutations[face_perm][2]];
  if (low[0] > high[0]) {
    swap = low[0];
    low[0] = high[0];
    high[0] = swap;
    edge_reverse[0] = 1;
  }
  if (low[1] > high[1]) {
    swap = low[1];
    low[1] = high[1];
    high[1] = swap;
    edge_reverse[1] = 1;
  }

  /* find matching target edges */
  for (i = 0; i < 12; ++i) {
    if (low[0] == p8est_edge_vertices[i][0] &&
        high[0] == p8est_edge_vertices[i][1]) {
      P4EST_ASSERT (target_axis[0] == -1);
      target_axis[0] = i / 4;
#ifndef P4EST_DEBUG
      if (target_axis[1] >= 0)
        break;
#endif
    }
    else if (low[1] == p8est_edge_vertices[i][0] &&
             high[1] == p8est_edge_vertices[i][1]) {
      P4EST_ASSERT (target_axis[1] == -1);
      target_axis[1] = i / 4;
#ifndef P4EST_DEBUG
      if (target_axis[0] >= 0)
        break;
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
  }
  P4EST_ASSERT (my_axis[0] != my_axis[1] &&
                my_axis[0] != my_axis[2] && my_axis[1] != my_axis[2]);
  P4EST_ASSERT (target_axis[0] != target_axis[1] &&
                target_axis[0] != target_axis[2] &&
                target_axis[1] != target_axis[2]);
#endif

  return target_tree;
}

void
p8est_find_edge_transform (p8est_connectivity_t * conn,
                           p4est_topidx_t itree, int iedge,
                           p8est_edge_info_t * ei)
{
  int                 redge, nedge, nflip;
  p4est_topidx_t      edge_trees, etree;
  p4est_topidx_t      aedge, ntree, ntree1, ntree2;
  p4est_topidx_t      v0, v1, nv0, nv1;
  p8est_edge_transform_t *et;
  sc_array_t         *ta = &ei->edge_transforms;

  P4EST_ASSERT (0 <= itree && itree < conn->num_trees);
  P4EST_ASSERT (0 <= iedge && iedge < 12);
  P4EST_ASSERT (ta->elem_size == sizeof (p8est_edge_transform_t));

  ei->iedge = (int8_t) iedge;
  ei->iflip = -1;
  sc_array_resize (ta, 0);
  if (conn->num_edges == 0) {
    return;
  }
  aedge = conn->tree_to_edge[12 * itree + iedge];
  if (aedge == -1) {
    return;
  }

  v0 = conn->tree_to_vertex[8 * itree + p8est_edge_vertices[iedge][0]];
  v1 = conn->tree_to_vertex[8 * itree + p8est_edge_vertices[iedge][1]];
  ntree1 = conn->tree_to_tree[6 * itree + ((iedge < 4) ? 2 : 0) + iedge % 2];
  ntree2 =
    conn->tree_to_tree[6 * itree + ((iedge < 8) ? 4 : 2) + (iedge / 2) % 2];
  edge_trees =                  /* same type */
    conn->ett_offset[aedge + 1] - conn->ett_offset[aedge];

  for (etree = 0; etree < edge_trees; ++etree) {
    ntree = conn->edge_to_tree[conn->ett_offset[aedge] + etree];
    redge = (int) conn->edge_to_edge[conn->ett_offset[aedge] + etree];
    nedge = redge % 12;
    nflip = redge / 12;
    if (nedge == iedge && ntree == itree) {
      P4EST_ASSERT (ei->iflip == -1);
      ei->iflip = (int8_t) nflip;
      continue;
    }
    if (ntree == ntree1 || ntree == ntree2) {
      nv0 = conn->tree_to_vertex[8 * ntree + p8est_edge_vertices[nedge][0]];
      nv1 = conn->tree_to_vertex[8 * ntree + p8est_edge_vertices[nedge][1]];
      if ((nv0 == v0 && nv1 == v1) || (nv0 == v1 && nv1 == v0)) {
        continue;
      }
    }

    /* else we have a true diagonal edge with ntree */
    et = sc_array_push (ta);
    et->ntree = ntree;
    et->nedge = (int8_t) nedge;
    et->naxis[0] = (int8_t) (nedge / 4);
    et->naxis[1] = (int8_t) (nedge < 4 ? 1 : 0);
    et->naxis[2] = (int8_t) (nedge < 8 ? 2 : 1);
    et->nflip = (int8_t) nflip;
    et->corners = (int8_t) (2 * ((nedge / 2) % 2) + nedge % 2);
  }
  P4EST_ASSERT (edge_trees == (p4est_topidx_t) ta->elem_count
                + 1 + (ntree1 != itree) + (ntree2 != itree));
  P4EST_ASSERT (ei->iflip >= 0);
}

/* EOF p8est_connectivity.h */
