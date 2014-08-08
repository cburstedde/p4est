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
