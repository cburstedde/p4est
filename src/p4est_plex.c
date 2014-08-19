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

#ifndef P4_TO_P8
#include <p4est_lnodes.h>
#else
#include <p8est_lnodes.h>
#endif

/** Decode the information from p{4,8}est_lnodes_t for a given element.
 *
 * \see p4est_lnodes.h for an in-depth discussion of the encoding.
 * \param [in] face_code         Bit code as defined in p{4,8}est_lnodes.h.
 * \param [out] hanging_corner   Undefined if no node is hanging.
 *                               If any node is hanging, this contains
 *                               one integer per corner, which is -1
 *                               for corners that are not hanging,
 *                               and the number of the non-hanging
 *                               corner on the hanging face/edge otherwise.
 *                               For faces in 3D, it is diagonally opposite.
 * \return true if any node is hanging, false otherwise.
 */
static int
lnodes_decode2 (p4est_lnodes_code_t face_code,
                int hanging_corner[P4EST_CHILDREN])
{
  if (face_code) {
    const int           c = (int) (face_code & ones);
    int                 i, h;
    int                 work = (int) (face_code >> P4EST_DIM);

    /* These two corners are never hanging by construction. */
    hanging_corner[c] = hanging_corner[c ^ ones] = -1;
    for (i = 0; i < P4EST_DIM; ++i) {
      /* Process face hanging corners. */
      h = c ^ (1 << i);
      hanging_corner[h ^ ones] = (work & 1) ? c : -1;
#ifdef P4_TO_P8
      /* Process edge hanging corners. */
      hanging_corner[h] = (work & P4EST_CHILDREN) ? c : -1;
#endif
      work >>= 1;
    }
    return 1;
  }
  return 0;
}

static void
mark_parent (p4est_locidx_t qid, int ctype_int, p4est_lnodes_code_t *F,
             p4est_locidx_t *quad_to_local, int8_t *is_parent,
             int8_t *node_dim)
{
#ifndef P4_TO_P8
  int                 dim_limits[3] = {0, 4, 8};
#else
  int                 dim_limits[4] = {0, 6, 18, 26};
#endif
  int hanging[3][12];
  int has_hanging;
  int V = dim_limits[ctype_int];

#ifndef P4_TO_P8
  has_hanging = p4est_lnodes_decode (F[qid], &hanging[0][0]);
#else
  has_hanging = p8est_lnodes_decode (F[qid], &hanging[0][0], &hanging[1][0]);
#endif
  if (has_hanging) {
    /* no corners */
    climit = SC_MIN (P4EST_DIM - 1, ctype_int);
    for (c = 0; c < climit; c++) {
      int v, vstart = dim_limits[c];
      int vend = dim_limits[c+1];

      for (v = vstart; v < vend; v++) {
        is_parent[quad_to_local[qid * V + v]] = 1;
      }
    }
  }
  for (c = 0; c < ctype_int; c++) {
    int v, vstart = dim_limits[c];
    int vend = dim_limits[c+1];

    for (v = vstart; v < vend; v++) {
      node_dim[quad_to_local[qid * V + v]] = P4EST_DIM - 1 - c;
    }
  }
}

static void
parent_to_child (p4est_quadrant_t *q, p4est_topidx_t, p4est_locidx_t qid, int ctype_int,
                 p4est_lnodes_code_t *F, p4est_locidx_t *quad_to_local,
                 p4est_gloidx_t *quad_to_gobal,
                 int8_t *is_parent, int8_t *node_dim,
                 p4est_locidx_t * child_offsets, int64_t *child_to_parent,
                 int64_t *child_to_id, p4est_connectivity_t *conn)
{
#ifndef P4_TO_P8
  int                 dim_limits[3] = {0, 4, 8};
#else
  int                 dim_limits[4] = {0, 6, 18, 26};
#endif
  int hanging[3][12];
  int has_hanging;
  int V = dim_limits[ctype_int];

#ifndef P4_TO_P8
  has_hanging = p4est_lnodes_decode (F[qid], &hanging[0][0]);
#else
  has_hanging = p8est_lnodes_decode (F[qid], &hanging[0][0], &hanging[1][0]);
#endif
  has_hanging |= lnodes_decode2 (F[qid], &hanging[P4EST_DIM-1][0]);
  if (has_hanging) {
    /* no corners */
    int cid = p4est_quadrant_child_id (q);
    for (c = 0; c < ctype_int; c++) {
      int vstart = dim_limits[c];
      int vend = dim_limits[c+1];

      if (!c) {

        for (v = vstart; v < vend; v++) {
          if (hanging[0][v] >= 0) {
            p4est_locidx_t child = child_offsets[quad_to_local[qid * V + v]] + hanging[0][v];
            quad_to_local[qid * V + v] = child;
            quad_to_global[qid * V + v] = -1;
          }
        }
      }
      else if (c == P4EST_DIM - 1) {
        for (v = vstart; v < vend; v++) {
          int corner = v-vstart;
          if (hanging[P4EST_DIM - 1][corner] >= 0) {
            p4est_locidx_t child;
            int f;

            f = p4est_child_corner_faces[cid][corner];
#ifndef P4_TO_P8
            P4EST_ASSERT (f >= 0);
#endif
            if (f >= 0) {
              child = child_offsets[quad_to_local[qid * V + f]] + P4EST_HALF;
              quad_to_local[qid * V + v] = child;
              quad_to_global[qid * V + v] = -1;
            }
#ifdef P4_TO_P8
            else {
              int e = p8est_child_corner_edges[cid][corner];

              P4EST_ASSERT (e >= 0);

              child = child_offsets[quad_to_local[qid * V + e + P4EST_FACES]] + 2;
              quad_to_local[qid * V + v] = child;
              quad_to_global[qid * V + v] = -1;
            }
#endif
          }
        }
      }
#ifdef P4_TO_P8
      else {
        for (v = vstart; v < vend; v++) {
          int edge = v-vstart;
          if (hanging[1][edge] >= 0) {
            p4est_locidx_t child;

            if (hanging[1][edge] < 4) {
              int h = hanging[1][edge] % 2;

              child = child_offsets[quad_to_local[qid * V + v]] + h;
              quad_to_local[qid * V + v] = child;
              quad_to_global[qid * V + v] = -1;
            }
            else {
              int i;

              for (i = 0; i < 2; i++) {
                int f = p8est_face_edges[edge][i];
                int c = p8est_corner_face_corners[cid][f];
                int e, j;
                p4est_quadrant_t tempq;

                if (c >= 0) {
                  int dir;
                  int hc = hanging[0][f];
                  int he[2];
                  int child_edge;

                  P4EST_ASSERT (hc >= 0);

                  for (j = 0; j < 4; j++) {
                    e = p8est_face_edges[f][j];

                    if (e == edge) {
                      break;
                    }
                  }
                  dir = j / 2;

                  he[0] = (hc & 1);
                  he[1] = 2 + ((hc & 2) >> 1);

                  /* flip dir if the orientation is different */
                  p4est_quadrant_face_neighbor (q, f, &tempq);
                  if (p4est_quadrant_is_outside_face (&tempq)) {
                    p4est_topidx_t nt;
                    int nf;

                    nt = conn->tree_to_tree[P4EST_FACES * t + f];
                    nf = conn->tree_to_face[P4EST_FACES * t + f];
                    if (nt != t && nf != f) {
                      if (nt < t || (nt == t && nf < f)) {
                        int o = nf / P4EST_FACES;
                        int ref = p8est_face_permutation_refs[nf][f];
                        int set = p8est_face_permutation_sets[ref][o];
                        if (set == 1 || set == 3 || set == 4 || set == 6) {
                          dir = (dir ^ 1);
                        }
                      }
                    }
                  }
                  if (!dir) {
                    /* first direction */
                    child_edge = he[0];
                  }
                  else {
                    /* second direction */
                    child_edge = he[1];
                  }

                  child = child_offsets[quad_to_local[qid * V + f]] + child_edge;
                  quad_to_local[qid * V + v] = child;
                  quad_to_global[qid * V + v] = -1;
                  break;
                }
              }
              P4EST_ASSERT (i < 2);
            }
          }
        }
      }
#endif
    }
  }
}

static void
p4est_get_plex_data_local_64bit (p4est_t *p4est, p4est_ghost_t *ghost,
                                 p4est_lnodes_t *lnodes,
                                 int overlap,
                                 int local_first)
{
  sc_array_t *cell_to_quad;
  sc_array_t *cell_to_proc;
  sc_array_t *cell_cones;
  sc_array_t *cell_orientations;
  sc_array_t *child_to_parent, *child_to_id;
  p4est_gloidx_t *quad_to_global;
  p4est_lnodes_code_t *F;
  p4est_locidx_t G = (p4est_locidx_t) ghost->ghosts.elem_count;
  p4est_locidx_t Klocal = p4est->local_num_quadrants;
  p4est_locidx_t K = Klocal + overlap ? G : 0;
  p4est_locidx_t Gpre, Gpost;
  p4est_locidx_t local_quad_offset;
  int mpirank = p4est->mpirank;
  int mpisize = p4est->mpisize;
  int qid, p;
  int V = lnodes->vnodes;
  sc_array_t *is_parent, *node_dim;
  int hanging[3][12];
  int ctype_int = p4est_connect_type_int (ghost->btype);
  p4est_locidx_t *child_offsets;

  Gpre = (overlap && local_first) ? (ghost->proc_offsets[mpirank] -
                                     ghost->proc_offsets[0]) : 0;
  Gpost = overlap ?
          (local_first ? ghost->proc_offsets[mpisize] -
           ghost->proc_offsets[mpirank + 1] : G) : 0;

  local_quad_offset = Gpre;

  cell_to_quad = sc_array_new_size (sizeof (int64_t), (size_t) K);
  cell_to_proc = sc_array_new_size (sizeof (int64_t), (size_t) K);

  for (qid = 0, p = 0; qid < Gpre; qid++) {
    p4est_quadrant_t *q;
    int64_t *gidx = (int64_t *) sc_array_index (cell_to_quad, (size_t) qid);
    int64_t *proc = (int64_t *) sc_array_index (cell_to_proc, (size_t) qid);
    q = p4est_quadrant_array_index (&ghost->ghosts, (size_t) qid);

    while (p < mpisize && ghost->proc_offsets[p+1] <= qid) {
      p++;
    }
    P4EST_ASSERT (ghost->proc_offsets[p] <= qid &&
                  ghost->proc_offstes[p+1] > qid);

    *gidx = (int64_t) (p4est->global_first_quadrant[p] + q->p.piggy3.local_num);
    *proc = (int64_t) p;
  }
  for (qid = 0; qid < Klocal; qid++) {
    int64_t *gidx = (int64_t *) sc_array_index (cell_to_quad, (size_t) (qid + local_quad_offset));
    int64_t *proc = (int64_t *) sc_array_index (cell_to_proc, (size_t) (qid + local_quad_offset));

    *gidx = (int64_t) (p4est->global_first_quadrant[mpirank] + qid);
    *proc = (int64_t) mpirank;
  }
  for (qid = 0, p = mpirank + 1; qid < Gpost; qid++) {
    p4est_quadrant_t *q;
    int64_t *gidx = (int64_t *) sc_array_index (cell_to_quad, (size_t) (qid + Gpre + Klocal));
    int64_t *proc = (int64_t *) sc_array_index (cell_to_proc, (size_t) (qid + Gpre + Klocal));
    q = p4est_quadrant_array_index (&ghost->ghosts, (size_t) qid + Gpre);

    while (p < mpisize && ghost->proc_offsets[p+1] <= qid + Gpre) {
      p++;
    }
    P4EST_ASSERT (ghost->proc_offsets[p] <= qid + Gpre &&
                  ghost->proc_offstes[p+1] > qid + Gpre);

    *gidx = (int64_t) (p4est->global_first_quadrant[p] + q->p.piggy3.local_num);
    *proc = (int64_t) p;
  }

  quad_to_global = P4EST_ALLOC (p4est_gloidx_t, K * V);
  if (overlap) {
    F = P4EST_ALLOC (p4est_lnodes_codes_t, K);
  }
  else {
    F = lnodes->face_code;
  }

  for (qid = 0; qid < Klocal * V; qid++) {
    quad_to_global[qid] = p4est_lnodes_global_index (lnodes, qid);
  }

  if (overlap) {
    p4est_gloidx_t    **mirror_data;
    p4est_locidx_t il, num_mirrors = (p4est_locidx_t) ghost->mirrors.elem_count;

    mirror_data = P4EST_ALLOC (p4est_gloidx_t *, ghost->mirrors.elem_count);
    for (il = 0; (size_t) il < num_mirrors; il++) {
      p4est_quadrant_t   *q;

      q = p4est_quadrant_array_index (&ghost->mirrors, il);

      qid = q->p.piggy3.local_num;
      mirror_data[il] = &quad_to_global[qid * vnodes];
    }
    p4est_ghost_exchange_custom (p4est, ghost,
                                 (size_t) V * sizeof (p4est_gloidx_t),
                                 (void **) mirror_data,
                                 &quad_to_global[K * vnodes]);
    P4EST_FREE (mirror_data);

    memcpy (F, p4est_lnodes, K * sizeof (p4est_lnodes_code_t));

    mirror_data = P4EST_ALLOC (p4est_lnodes_code_t *, ghost->mirrors.elem_count);
    for (il = 0; (size_t) il < num_mirrors; il++) {
      mirror_data[il] = &F[qid];
    }
    p4est_ghost_exchange_custom (p4est, ghost,
                                 sizeof (p4est_lnodes_code_t),
                                 (void **) mirror_data,
                                 &F[K]);
    P4EST_FREE (mirror_data);
  }

  /* create a list of all global nodes seen */
  all_global = sc_array_new_size (2 * sizeof (p4est_gloidx_t), V * K);

  for (qid = 0; qid < K; qid++) {
    for (v = 0; v < V; v++) {
      p4est_gloidx_t *pair = (p4est_gloidx_t *) sc_array_index (all_global, (size_t) qid * V + v);

      pair[0] = quad_to_global[qid * V + v];
      pair[1] = qid * V + v;
    }
  }

  sc_array_sort (all_global, p4est_gloidx_compare);

  num_global = 0;
  last_global = -1;

  quad_to_local = P4EST_ALLOC (p4est_locidx_t, K * V);
  for (zz = 0; zz < K * V; zz++) {
    p4est_gloidx_t *pair = (p4est_gloidx_t *) sc_array_index (all_global, zz);
    p4est_gloidx_t gidx;

    gidx = pair[0];
    if (gidx != last_global) {
      num_global++;
      last_global = gidx;
    }
    quad_to_local[pair[1]] = num_global - 1;
  }
#ifdef P4EST_DEBUG
  sc_array_uniq (all_global, p4est_gloidx_compare);
  P4EST_ASSERT (all_global->elem_count == (size_t) num_global);
#endif
  sc_array_destroy (all_global);

  /* figure out which nodes are parents */
  is_parent = sc_array_new_size (sizeof (int8_t), num_global);
  node_dim = sc_array_new_size (sizeof (int8_t), num_global);
  for (qid = 0; qid < K; qid++) {
    mark_parent (qid, ctype_int, F, quad_to_local, (int8_t *) is_parent->array,
                 (int8_t *) node_dim->array);
  }
  child_offsets = P4EST_ALLOC (p4est_locidx_t, num_global+1);

  child_offsets[0] = num_global;
  for (il = 0; il < num_global; il++) {
    int8_t parent = *((int8_t *) sc_array_index (is_parent, (size_t) il));
    int8_t dim = *((int8_t *) sc_array_index (node_dim, (size_t) il));
    int count;
    if (!parent) {
      count = 0;
    }
    else {
      P4EST_ASSERT (dim == 1 || dim == 2);
      if (dim == 1) {
        count = 3;
      }
      else {
        count = 9;
      }
    }
    child_offsets[il+1] = child_offsets[il] + count;
  }

  num_global_plus_children = child_offsets[num_global];

  sc_array_resize (is_parent, num_global_plus_children);
  sc_array_resize (node_dim, num_global_plus_children);
  child_to_parent = sc_array_new (sizeof (int64_t), num_global_plus_children);
  memset (child_to_parent->array, -1, child_to_parent->elem_count * child_to_parent->elem_size);
  child_to_id = sc_array_new (sizeof (int64_t), num_global_plus_children);

  for (il = num_global; il < num_global_plus_children; il++) {
    p4est_locidx_t jstart, jend, jl;
    int8_t pdim = *((int8_t *) sc_array_index (node_dim, (size_t) il));

    jstart = child_offsets[il];
    jend = child_offsets[il+1];

    for (jl = jstart; jl < jend; jl++) {
      int8_t *parent = (int8_t *) sc_array_index (is_parent, (size_t) jl);
      int8_t *dim = (int8_t *) sc_array_index (node_dim, (size_t) jl);
      int64_t *parentidx = (int64_t *) sc_array_index (child_to_parent, (size_t) jl);
      int64_t *childid = (int64_t *) sc_array_index (child_to_id, (size_t) jl);

      *parentidx = (int64_t) il;
      *parent = 0;
      *child_to_id = (int64_t) (jl - jstart);
      if (pdim == 1) {
        if (jl - jstart < 2) {
          *dim = 1;
        }
        else {
          *dim = 0;
        }
      }
      else {
        if (jl - jstart < 4) {
          *dim = 2;
        }
        else if (jl - jstart < 8) {
          *dim = 1;
        }
        else {
          *dim = 0;
        }
      }
    }
  }

  /* put children in quad_to_local */
  for (qid = 0, t = flt; t <= llt; t++) {
    p4est_tree_t * tree = p4est_tree_array_index (p4est->trees, t);
    sc_array_t *quadrants = &(tree->quadrants);
    p4est_locidx_t il, num_quads = (p4est_locidx_t) quadrants->elem_count;

    for (il = 0; il < num_quads; il++, qid++) {
      p4est_quadrant_t *q = p4est_quadrant_array_index (quadrants, (size_t) il);

      parent_to_child (q, t, qid, ctype_int, F, quad_to_local, quad_to_global,
                       (int8_t *) is_parent->array,
                       (int8_t *) node_dim->array, child_offsets,
                       (int64_t *) child_to_parent->array,
                       (int64_t *) child_to_id->array,
                       p4est->connectivity);
    }
  }
  for (t = 0; t < p4est->connectivity->num_trees; t++) {
    p4est_locidx_t il, istart = ghost->tree_offsets[t];
    p4est_locidx_t iend = ghost->tree_offsets[t+1];

    for (il = istart; il < iend; il++) {
      p4est_quadrant_t *q = p4est_quadrant_array_index (&ghost->ghosts, (size_t) il);

      parent_to_child (q, t, qid + Klocal, ctype_int, F, quad_to_local, quad_to_global,
                       (int8_t *) is_parent->array,
                       (int8_t *) node_dim->array, child_offsets,
                       (int64_t *) child_to_parent->array,
                       (int64_t *) child_to_id->array,
                       p4est->connectivity);
    }
  }

  /* now we pop up the cones */

  /* now we permute to stratify */
}

static int
p4est_gloidx_compare_duplex (const void *A, const void *B)
{
  const p4est_gloidx_t *a = (const p4est_gloidx_t *)A;
  const p4est_gloidx_t *b = (const p4est_gloidx_t *)B;
  p4est_gloidx_t diff;

  int ret = p4est_gloidx_compare (A, B);
  if (ret) return ret;
  diff = a[1] - b[1];
  if (!diff) return 0;
  return diff > 0 ? 1 : -1;
}

static void
p4est_get_plex_data_64bit (p4est_t *p4est, p4est_connect_type_t ctype, int overlap)
{
  p4est_ghost_t *ghost;
  p4est_lnodes_t *lnodes;
  int ctype_int = p4est_connect_type_int (ctype);
  int o;
  p4est_lnodes_code_t *F;
  p4est_locidx_t     *quad_to_node_local;
  p4est_gloidx_t     *combined_offset;
  p4est_locidx_t     *dimension_offsets;
  p4est_gloidx_t     *node_global;
  p4est_locidx_t      G, N, K, KG;
  int mpisize = p4est->mpisize;
  int mpirank = p4est->mpirank;

  /* get the p4est data */
  ghost = p4est_ghost_new (p4est, ctype);
  lnodes = p4est_lnodes_new (p4est, ghost, -ctype_int);
  if (overlap > 0) {
    p4est_ghost_support_lnodes (p4est, lnodes, ghost);
  }
  for (o = 1; o < overlap; o++) {
    p4est_ghost_expand_by_lnodes (p4est, lnodes, ghost);
  }

  G = (p4est_locidx_t) ghost->ghosts.elem_count;
  K = lnodes->num_local_elements;
  N = lnodes->num_local_nodes;

  /* get the local to global numbering for the nodes (stratified),
   * and the element to local numbering that takes into account new nodes
   * brought in with the ghost layer */
  {
    p4est_gloidx_t     *quad_to_node_global;
    p4est_gloidx_t    **mirror_data;
    p4est_locidx_t      il, num_mirrors, qid;
    int                 v, nid;
    int                 vnodes = lnodes->vnodes;
    int                *node_counted;
    p4est_locidx_t     *dimension_counts;
#ifndef P4_TO_P8
    int                 dim_limits[3] = {0, 4, 8};
#else
    int                 dim_limits[4] = {0, 6, 18, 26};
#endif
    int p, d;

    combined_offset = P4EST_ALLOC (p4est_gloidx_t, p4est->mpisize + 1);
    combined_offset[0] = 0;
    for (p = 0; p < mpisize; p++) {
      combined_offset[p + 1] = combined_offset[p] +
                               (p4est->global_first_quadrant[p + 1] -
                                p4est->global_first_quadrant[p]) +
                               lnodes->global_owned_count[p];
    }

    dimension_counts = P4EST_ALLOC_ZERO (p4est_locidx_t, ctype_int);
    node_counted = P4EST_ALLOC_ZERO (int, lnodes->owned_count);

    for (qid = 0; qid < K; qid++) {
      for (d = 0; d < ctype_int; d++) {
        int vstart, vend;

        vstart = dim_limits[d];
        vend = dim_limits[d + 1];
        for (v = vstart; v < vend; v++) {
          nid = lnodes->element_nodes[qid * vnodes + v];
          if (nid < lnodes->owned_count && !node_counted[nid]) {
            node_counted[nid] = 1;
            dimension_counts[d]++;
          }
        }
      }
    }
    P4EST_FREE (node_counted);

    dimension_offsets = P4EST_ALLOC_ZERO (p4est_locidx_t, ctype_int + 1);
    dimension_offsets[0] = p4est->local_num_quadrants;
    for (d = 0; d < ctype_int; d++) {
      dimension_offsets[d + 1] = dimension_offsets[d] + dimension_counts[d];
    }

    node_global = P4EST_ALLOC (p4est_gloidx_t, N);

    for (qid = 0; qid < K; qid++) {
      for (d = 0; d < ctype_int; d++) {
        int vstart, vend;

        vstart = dim_limits[d];
        vend = dim_limits[d + 1];
        for (v = vstart; v < vend; v++) {
          nid = lnodes->element_nodes[qid * vnodes + v];
          if (nid < lnodes->owned_count) {
            node_global[nid] = d;
          }
        }
      }
    }

    for (d = 0; d < ctype_int; d++) {
      dimension_counts[d] = 0;
    }

    for (nid = 0; nid < lnodes->owned_count; nid++) {
      d = node_global[nid];
      node_global[nid] = combined_offset[mpirank] +
                         dimension_offsets[d] +
                         dimension_counts[d]++;
    }
    P4EST_FREE (dimension_counts);

    /* communicate */
    {
      sc_array_t ngarray;

      sc_array_init_data (&ngarray, node_global, sizeof (p4est_gloidx_t),
                          (size_t) N);
      p4est_lnodes_share_owned (&ngarray, lnodes);
    }

    if (!overlap) {
      KG = K;
      F = lnodes->face_code;
    }
    else {
      KG = K + G;
      F = P4eST_ALLOC (p4est_lnodes_code_t, KG);
    }
    quad_to_node_global = P4EST_ALLOC (p4est_gloidx_t, KG * vnodes);

    /* copy into the element arrays */
    for (qid = 0; qid < K * vnodes; qid++) {
      nid = lnodes->element_nodes[qid];
      quad_to_node_global[qid] = node_global[nid];
    }

    if (overlap) {
      num_mirrors = (p4est_locidx_t) ghost->mirrors.elem_count;

      mirror_data = P4EST_ALLOC (p4est_gloidx_t *, ghost->mirrors.elem_count);
      for (il = 0; (size_t) il < num_mirrors; il++) {
        p4est_quadrant_t   *q;

        q = p4est_quadrant_array_index (&ghost->mirrors, il);

        qid = q->p.piggy3.local_num;
        mirror_data[il] = &quad_to_node_global[qid * vnodes];
      }
      p4est_ghost_exchange_custom (p4est, ghost,
                                   (size_t) vnodes * sizeof (p4est_gloidx_t),
                                   (void **) mirror_data,
                                   &quad_to_node_global[K * vnodes]);
      P4EST_FREE (mirror_data);

      memcpy (F, p4est_lnodes, K * sizeof (p4est_lnodes_code_t));

      mirror_data = P4EST_ALLOC (p4est_lnodes_code_t *, ghost->mirrors.elem_count);
      for (il = 0; (size_t) il < num_mirrors; il++) {
        mirror_data[il] = &F[qid];
      }
      p4est_ghost_exchange_custom (p4est, ghost,
                                   sizeof (p4est_lnodes_code_t),
                                   (void **) mirror_data,
                                   &F[K]);
      P4EST_FREE (mirror_data);
    }

    {
      sc_array_t *all_nodes, view;

      all_nodes = sc_array_new_size (2 * sizeof (p4est_gloidx_t), KG * vnodes);
      for (qid = 0; qid < KG; qid++) {
        for (d = 0; d < ctype_int; d++) {
          int vstart = dim_limits[d];
          int vend = dim_limits[d+1];

          for (v = vstart; v < vend; v++) {
            p4est_gloidx_t *idxandtype = (p4est_gloidx_t *) sc_array_index (all_nodes, qid * vnodes + v);

            idxandtype[0] = (p4est_gloidx_t) d;
            idxandtype[1] = quad_to_node_global[qid * vnodes + v];
          }
        }
      }
      sc_array_sort (all_nodes, p4est_gloidx_compare_duplex);
      sc_array_uniq (all_nodes, p4est_gloidx_compare_duplex);

      P4EST_FREE (node_global);

      P4EST_ASSERT (all_nodes->elem_count >= N);
      N = all_nodes->elem_count;
      node_global = P4EST_ALLOC (p4est_gloidx_t, N);
      for (nid = 0; nid < N; nid++) {
        p4est_gloidx_t *idxandtype = (p4est_gloidx_t *) sc_array_index (all_nodes, nid);

        node_global[nid] = idxandtype[1];
      }
      sc_array_destroy (all_nodes);
      quad_to_node_local = P4EST_ALLOC (p4est_locidx_t, N);

      sc_array_init_data (&view, node_global, sizeof (p4est_gloidx_t), N);
      for (qid = 0; qid < KG * vnodes; qid++) {
        ssize_t idx;
        p4est_gloidx_t gid = quad_to_node_global[qid];

        idx = sc_array_bsearch (&view, &gid, p4est_gloidx_compare);
        P4EST_ASSERT (idx >= 0);
        quad_to_node_local[qid] = (p4est_locidx_t) idx;
      }
    }
    P4EST_FREE (quad_to_node_global);
  }

  /* destroy the p4est data */
  p4est_lnodes_destroy (lnodes);
  p4est_ghost_destroy (ghost);
}
