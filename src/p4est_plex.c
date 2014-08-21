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
#include <p4est_bits.h>
#include <p4est_lnodes.h>
#include <p4est_plex.h>
#else
#include <p8est_bits.h>
#include <p8est_lnodes.h>
#include <p8est_plex.h>
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
  int ones = P4EST_CHILDREN - 1;
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
  int hanging[2][12];
  int has_hanging;
  int c, V = dim_limits[ctype_int];

#ifndef P4_TO_P8
  has_hanging = p4est_lnodes_decode (F[qid], &hanging[0][0]);
#else
  has_hanging = p8est_lnodes_decode (F[qid], &hanging[0][0], &hanging[1][0]);
#endif
  if (has_hanging) {
    int climit;

    /* no corners */
    climit = SC_MIN (P4EST_DIM - 1, ctype_int);
    for (c = 0; c < climit; c++) {
      int v, vstart = dim_limits[c];
      int vend = dim_limits[c+1];

      for (v = vstart; v < vend; v++) {
        if (hanging[c][v - vstart] >= 0) {
          is_parent[quad_to_local[qid * V + v]] = 1;
        }
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
parent_to_child (p4est_quadrant_t *q, p4est_topidx_t t, p4est_locidx_t qid, int ctype_int,
                 p4est_lnodes_code_t *F, p4est_locidx_t *quad_to_local,
                 int8_t *quad_to_orientations,
                 int8_t *referenced,
                 int8_t *node_dim,
                 p4est_locidx_t *child_offsets,
                 p4est_locidx_t *child_to_id, p4est_connectivity_t *conn)
{
#ifndef P4_TO_P8
  int                 dim_limits[3] = {0, 4, 8};
  int                 no = P4EST_FACES;
#else
  int                 dim_limits[4] = {0, 6, 18, 26};
  int                 e, no = P4EST_FACES + P8EST_EDGES;
#endif
  int hanging[3][12];
  int has_hanging;
  int f, V = dim_limits[ctype_int];

#ifndef P4_TO_P8
  has_hanging = p4est_lnodes_decode (F[qid], &hanging[0][0]);
#else
  has_hanging = p8est_lnodes_decode (F[qid], &hanging[0][0], &hanging[1][0]);
#endif
  has_hanging |= lnodes_decode2 (F[qid], &hanging[P4EST_DIM-1][0]);
  if (has_hanging) {
    int c, cid = p4est_quadrant_child_id (q), v;

    /* no corners */
    for (c = ctype_int - 1; c >= 0; c--) {
      int vstart = dim_limits[c];
      int vend = dim_limits[c+1];

      if (!c) {
        for (v = vstart; v < vend; v++) {
          if (hanging[0][v] >= 0) {
            p4est_locidx_t child = child_offsets[quad_to_local[qid * V + v]] + hanging[0][v];
            quad_to_local[qid * V + v] = child;
            referenced[qid * V + v] = 1;
          }
        }
      }
      else if (c == P4EST_DIM - 1) {
        for (v = vstart; v < vend; v++) {
          int corner = v-vstart;
          if (hanging[P4EST_DIM - 1][corner] >= 0) {
            p4est_locidx_t child = -1;
            int dim;

            f = p4est_child_corner_faces[cid][corner];
            P4EST_ASSERT (P4EST_DIM == 3 || f >= 0);
            if (f >= 0) {
              dim = P4EST_DIM - 1;
              child = child_offsets[quad_to_local[qid * V + f]];
            }
#ifdef P4_TO_P8
            else {
              int e = p8est_child_corner_edges[cid][corner];

              P4EST_ASSERT (e >= 0);
              dim = 1;
              child = child_offsets[quad_to_local[qid * V + e + P4EST_FACES]];
            }
#endif
            P4EST_ASSERT (dim == 1 || dim == 2);
            child += (dim == 1) ? 2 : 8;
            quad_to_local[qid * V + v] = child;
            referenced[qid * V + v] = 1;
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
              referenced[qid * V + v] = 1;
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
                  P4EST_ASSERT (j < 4);
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
                  referenced[qid * V + v] = 1;
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
  for (f = 0; f < P4EST_FACES; f++) {
    p4est_quadrant_t tempq;

    p4est_quadrant_face_neighbor (q, f, &tempq);
    quad_to_orientations[qid * no + f] = 0;
    if (p4est_quadrant_is_outside_face (&tempq)) {
      p4est_topidx_t nt;
      int nf;

      nt = conn->tree_to_tree[P4EST_FACES * t + f];
      nf = conn->tree_to_face[P4EST_FACES * t + f];
      if (nt != t && nf != f) {
        if (nt < t || (nt == t && nf < f)) {
          int set, o = nf / P4EST_FACES;
#ifdef P4_TO_P8
          int ref;
#endif
          o = nf / P4EST_FACES;

#ifndef P4_TO_P8
          set = o;
#else
          ref = p8est_face_permutation_refs[nf][f];
          set = p8est_face_permutation_sets[ref][o];
#endif
          quad_to_orientations[qid * no + f] = set;
        }
      }
    }
  }
#ifdef P4_TO_P8
  for (e = 0; e < P4EST_FACES; e++) {
    p4est_quadrant_t tempq;

    p4est_quadrant_face_neighbor (q, e, &tempq);
    quad_to_orientations[qid * no + P4EST_FACES + e] = 0;
    if (p4est_quadrant_is_outside_face (&tempq)) {
      int set;
      int i, f = -1;
      int cid[2];

      for (i = 0; i < 2; i++) {
        int dir;
        p4est_qcoord_t d = 0;

        f = p8est_edge_faces[e][i];
        dir = f / 2;
        switch (dir) {
        case 0:
          d = tempq.x;
          break;
        case 1:
          d = tempq.y;
          break;
        case 2:
          d = tempq.z;
          break;
        default:
          SC_ABORT_NOT_REACHED();
          break;
        }
        if (d < 0 || d >= P4EST_ROOT_LEN) {
          break;
        }
      }
      P4EST_ASSERT (f >= 0);

      set = quad_to_orientations[qid * no + f];
      for (i = 0; i < 2; i++) {
        int c, face_ex, face_in;

        c = p8est_edge_corners[e][0];
        face_ex = p8est_corner_face_corners[c][f];
        P4EST_ASSERT (face_ex >= 0);
        face_in = p8est_face_permutations[set][face_ex];
        cid[i] = face_in;
      }

      if (cid[0] < cid[1]) {
        quad_to_orientations[qid * no + P4EST_FACES + e] = 0;
      }
      else {
        quad_to_orientations[qid * no + P4EST_FACES + e] = 1;
      }
    }
    else if (p8est_quadrant_is_outside_edge (&tempq)) {
      int edge = conn->tree_to_edge ? conn->tree_to_edge[t * P8EST_EDGES + e] : -1;

      if (edge >= 0) {
        int estart, eend, i;

        estart = conn->ett_offset[edge];
        eend = conn->ett_offset[edge+1];
        for (i = estart; i < eend; i++) {
          p4est_topidx_t nt;
          int8_t te;


          nt = conn->edge_to_tree[i];
          te = conn->edge_to_edge[i];
          if (nt == t && (te % P8EST_EDGES == e)) {
            quad_to_orientations[qid * no + P4EST_FACES + e] = te / P8EST_EDGES;
            break;
          }
        }
        P4EST_ASSERT (i < eend);
      }
      else {
        p4est_locidx_t ownt = t;
        int owne = e;
        int i, j, o = 0;
        for (i = 0; i < 2; i++) {
          p4est_locidx_t nt;
          int8_t nf;
          int fo;
          int ref, set;
          int cid[2];
          int ne;

          f = p8est_face_edges[e][i];
          nt = conn->tree_to_tree[P4EST_FACES * t + f];
          nf = conn->tree_to_face[P4EST_FACES * t + f];
          fo = nf / P8EST_FACES;
          nf = nf % P8EST_FACES;

          ref = p8est_face_permutation_refs[f][nf];
          set = p8est_face_permutation_sets[ref][fo];
          for (j = 0; j < 2; j++) {
            int c, face_ex, face_in;

            c = p8est_edge_corners[e][0];
            face_ex = p8est_corner_face_corners[c][f];
            P4EST_ASSERT (face_ex >= 0);
            face_in = p8est_face_permutations[set][face_ex];
            cid[j] = p8est_face_corners[nf][face_in];
          }
          ne = p8est_child_corner_edges[cid[0]][cid[1]];
          P4EST_ASSERT (ne >= 0);
          if (nt < ownt || (nt == ownt && ne < owne)) {
            ownt = nt;
            owne = e;
            o = (cid[0] < cid[1]) ? 0 : 1;
          }
        }
        quad_to_orientations[qid * no + P4EST_FACES + e] = o;
      }
    }
    {
      p4est_locidx_t childid = child_to_id[quad_to_local[qid * V + P4EST_FACES
                               + e]];
      if (childid >= 4 && childid & 1) {
        quad_to_orientations[qid * no + P4EST_FACES + e] ^= 1;
      }
    }
  }
#endif
}

static int
p4est_locidx_compare_double (const void *A, const void *B)
{
  const p4est_locidx_t *a = (const p4est_locidx_t *) A;
  const p4est_locidx_t *b = (const p4est_locidx_t *) A;
  int diff;

  diff = p4est_locidx_compare (A, B);
  if (diff) {
    return diff;
  }
  else {
    return p4est_locidx_compare (&a[1], &b[1]);
  }
}

#ifndef P4_TO_P8
static int p4est_to_plex_child_id[1][3] = {{9, 10, 25}};
static int p4est_to_plex_face_orientation[4][2] = {{-2,  0},
                                                   { 0, -2},
                                                   { 0, -2},
                                                   {-2,  0}};
static int p4est_to_plex_position[1][4] = {{3, 2, 0, 1}};
#else
static int p4est_to_plex_child_id[2][9] = {{15, 16, 18, 17, 90, 88, 87, 89, 137},
                                           {63, 64, 125, -1, -1, -1, -1, -1, -1}};
static int p4est_to_plex_face_orientation[6][8] = {{-1,  0,  3, -2, -4,  1,  2, -3},
                                                   { 0, -1, -2,  3,  1, -4, -3,  2},
                                                   { 0, -1, -2,  3,  1, -4, -3,  2},
                                                   {-2,  1,  0, -3, -1,  2,  3, -4},
                                                   {-1,  0,  3, -2, -4,  1,  2, -3},
                                                   { 0, -1, -2,  3,  1, -4, -3,  2}};
static int p4est_to_plex_edge_orientation[4][2] = {{-2,  0},
                                                   { 0, -2},
                                                   { 0, -2},
                                                   {-2,  0}};
static int p4est_to_plex_position[2][6] = {{4, 5, 5, 4, 0, 1},
                                           {3, 2, 0, 1, -1, -1}};
#endif

static void
p4est_get_plex_data_int (p4est_t *p4est, p4est_ghost_t *ghost,
                         p4est_lnodes_t *lnodes, int overlap, int local_first,
                         p4est_locidx_t *first_local_quad,
                         sc_array_t *out_points_per_dim,
                         sc_array_t *out_cone_sizes,
                         sc_array_t *out_cones,
                         sc_array_t *out_cone_orientations,
                         sc_array_t *out_vertex_coords,
                         sc_array_t *out_children,
                         sc_array_t *out_parents,
                         sc_array_t *out_childids,
                         sc_array_t *out_leaves,
                         sc_array_t *out_remotes)
{
#ifndef P4_TO_P8
  int                 dim_limits[3] = {0, 4, 8};
  int                 no = P4EST_FACES;
#else
  int                 dim_limits[4] = {0, 6, 18, 26};
  int                 no = P4EST_FACES + P8EST_EDGES;
#endif
  p4est_locidx_t *cones;
  int *orientations;
  sc_array_t *child_to_parent, *child_to_id;
  p4est_locidx_t *quad_to_local;
  int8_t *quad_to_orientations;
  int8_t *referenced;
  p4est_lnodes_code_t *F;
  p4est_topidx_t t, flt = p4est->first_local_tree;
  p4est_topidx_t llt = p4est->last_local_tree;
  p4est_locidx_t il, G = (p4est_locidx_t) ghost->ghosts.elem_count;
  p4est_locidx_t Klocal = p4est->local_num_quadrants;
  p4est_locidx_t K = Klocal + overlap ? G : 0;
  p4est_locidx_t Gpre, Gpost;
  int p, mpirank = p4est->mpirank;
  int mpisize = p4est->mpisize;
  p4est_locidx_t qid;
  int v, V = lnodes->vnodes;
  sc_array_t *is_parent, *node_dim;
  int ctype_int = p4est_connect_type_int (P4EST_CONNECT_FULL);
  p4est_locidx_t num_global, num_global_plus_children, last_global, *child_offsets;
  sc_array_t *all_global;
  p4est_locidx_t num_mirrors = (p4est_locidx_t) ghost->mirrors.elem_count;

  P4EST_ASSERT (lnodes->degree == -ctype_int);

  if (overlap) {
    /* get the face codes for ghosts */
    p4est_lnodes_code_t **mirror_data_F;

    F = P4EST_ALLOC (p4est_lnodes_code_t, K);
    memcpy (F, lnodes, K * sizeof (p4est_lnodes_code_t));

    mirror_data_F = P4EST_ALLOC (p4est_lnodes_code_t *, num_mirrors);
    for (il = 0; (size_t) il < num_mirrors; il++) {
      p4est_quadrant_t   *q;

      q = p4est_quadrant_array_index (&ghost->mirrors, il);
      qid = q->p.piggy3.local_num;
      mirror_data_F[il] = &F[qid];
    }
    p4est_ghost_exchange_custom (p4est, ghost,
                                 sizeof (p4est_lnodes_code_t),
                                 (void **) mirror_data_F,
                                 &F[Klocal]);
    P4EST_FREE (mirror_data_F);
  }
  else {
    F = lnodes->face_code;
  }

  /* we use the existing lnodes global indices as our anchors to keep all of
   * the new indices straight */
  /* create a list of all global nodes seen */
  {
    p4est_gloidx_t *quad_to_global = P4EST_ALLOC (p4est_gloidx_t, K * V);

    /* fill the local portion of quad_to_global */
    for (qid = 0; qid < Klocal * V; qid++) {
      quad_to_global[qid] = p4est_lnodes_global_index (lnodes, qid);
    }
    if (overlap) {
      /* get the ghost portion of quad_to_global */
      p4est_gloidx_t    **mirror_data;

      mirror_data = P4EST_ALLOC (p4est_gloidx_t *, num_mirrors);
      for (il = 0; (size_t) il < num_mirrors; il++) {
        p4est_quadrant_t   *q;

        q = p4est_quadrant_array_index (&ghost->mirrors, il);
        qid = q->p.piggy3.local_num;
        mirror_data[il] = &quad_to_global[qid * V];
      }
      p4est_ghost_exchange_custom (p4est, ghost,
                                   (size_t) V * sizeof (p4est_gloidx_t),
                                   (void **) mirror_data,
                                   &quad_to_global[Klocal * V]);
      P4EST_FREE (mirror_data);
    }

    all_global = sc_array_new_size (2 * sizeof (p4est_gloidx_t), V * K);
    for (qid = 0; qid < K; qid++) {
      for (v = 0; v < V; v++) {
        p4est_gloidx_t *pair = (p4est_gloidx_t *) sc_array_index (all_global, (size_t) qid * V + v);

        pair[0] = quad_to_global[qid * V + v];
        pair[1] = qid * V + v;
      }
    }
    P4EST_FREE (quad_to_global);
  }

  /* assign a local index to each global node referenced */
  sc_array_sort (all_global, p4est_gloidx_compare);
  quad_to_local = P4EST_ALLOC (p4est_locidx_t, K * V);
  num_global = 0;
  last_global = -1;
  for (il = 0; il < K * V; il++) {
    p4est_gloidx_t *pair = (p4est_gloidx_t *) sc_array_index (all_global, (size_t) il);
    p4est_gloidx_t gidx;

    gidx = pair[0];
    if (gidx != last_global) {
      num_global++;
      last_global = gidx;
    }
    quad_to_local[pair[1]] = num_global - 1;
  }
  /* remove duplicates from all_global so we can use it to search for the
   * local index of a global index */
  sc_array_uniq (all_global, p4est_gloidx_compare);
  P4EST_ASSERT (all_global->elem_count == (size_t) num_global);

  /* in lnodes, hanging faces/edges/corners are not assigned indices, they are
   * simply references to the anchor points with a hanging type arrow.  In
   * plex, however, all points have indices, so we have to expand these points
   * at the end of the list of local nodes (they do not, however, get global
   * indices) */
  /* figure out which nodes are parents and mark node dimensions */
  is_parent = sc_array_new_size (sizeof (int8_t), num_global);
  node_dim = sc_array_new_size (sizeof (int8_t), num_global);
  memset (is_parent->array, 0, is_parent->elem_count * is_parent->elem_size);
  for (qid = 0; qid < K; qid++) {
    mark_parent (qid, ctype_int, F, quad_to_local, (int8_t *) is_parent->array,
                 (int8_t *) node_dim->array);
  }
  child_offsets = P4EST_ALLOC (p4est_locidx_t, num_global+1);
  /* childredn are appended to the list of global nodes */
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
  sc_array_destroy (is_parent);
  num_global_plus_children = child_offsets[num_global];

  /* expand the node_dim array so that they also have entries
   * for children */
  sc_array_resize (node_dim, num_global_plus_children);
  child_to_id = sc_array_new_size (sizeof (p4est_locidx_t), num_global_plus_children);
  child_to_parent = sc_array_new_size (sizeof (p4est_locidx_t), num_global_plus_children);
  memset (child_to_id->array, -1, child_to_parent->elem_count * child_to_parent->elem_size);
  memset (child_to_parent->array, -1, child_to_parent->elem_count * child_to_parent->elem_size);
  /* fill the child_to_parent, child_to_id, and node_dim arrays for
   * the children */
  for (il = 0; il < num_global; il++) {
    p4est_locidx_t jstart, jend, jl;
    int8_t pdim = *((int8_t *) sc_array_index (node_dim, (size_t) il));

    jstart = child_offsets[il];
    jend = child_offsets[il+1];
    for (jl = jstart; jl < jend; jl++) {
      int8_t *dim = (int8_t *) sc_array_index (node_dim, (size_t) jl);
      p4est_locidx_t *parentidx = (p4est_locidx_t *) sc_array_index (child_to_parent, (size_t) jl);
      p4est_locidx_t *childid = (p4est_locidx_t *) sc_array_index (child_to_id, (size_t) jl);

      *parentidx = (p4est_locidx_t) il;
      *childid = (p4est_locidx_t) (jl - jstart);
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

  /* loop over quads:
   * - where quad_to_local refers to a parent, make it refer to the correct
   *   child
   * - mark children that are actually referenced (some may not be)
   * - fill quad_to_orientations
   */
  quad_to_orientations = P4EST_ALLOC (int8_t, K * no);
  referenced = P4EST_ALLOC_ZERO (int8_t, num_global_plus_children);
  for (il = 0; il < num_global; il++) {
    referenced[il] = 1;
  }
  for (qid = 0, t = flt; t <= llt; t++) {
    p4est_tree_t * tree = p4est_tree_array_index (p4est->trees, t);
    sc_array_t *quadrants = &(tree->quadrants);
    p4est_locidx_t num_quads = (p4est_locidx_t) quadrants->elem_count;

    for (il = 0; il < num_quads; il++, qid++) {
      p4est_quadrant_t *q = p4est_quadrant_array_index (quadrants, (size_t) il);

      parent_to_child (q, t, qid, ctype_int, F, quad_to_local,
                       quad_to_orientations, referenced,
                       (int8_t *) node_dim->array, child_offsets,
                       (p4est_locidx_t *) child_to_id->array,
                       p4est->connectivity);
    }
  }
  if (overlap) {
    for (t = 0; t < p4est->connectivity->num_trees; t++) {
      p4est_locidx_t il, istart = ghost->tree_offsets[t];
      p4est_locidx_t iend = ghost->tree_offsets[t+1];

      for (il = istart; il < iend; il++) {
        p4est_quadrant_t *q = p4est_quadrant_array_index (&ghost->ghosts, (size_t) il);

        parent_to_child (q, t, il + Klocal, ctype_int, F, quad_to_local,
                         quad_to_orientations, referenced,
                         (int8_t *) node_dim->array, child_offsets,
                         (p4est_locidx_t *) child_to_id->array,
                         p4est->connectivity);
      }
    }
    P4EST_FREE (F);
  }
  P4EST_FREE (child_offsets);

  /* compress unreferenced children out of the local list */
  {
    p4est_locidx_t *old_to_new, *new_to_old;
    p4est_locidx_t new_count, diff;

    old_to_new = P4EST_ALLOC (p4est_locidx_t, num_global_plus_children);
    new_to_old = P4EST_ALLOC (p4est_locidx_t, num_global_plus_children);

    memset (old_to_new, -1, num_global_plus_children * sizeof (*old_to_new));
    memset (new_to_old, -1, num_global_plus_children * sizeof (*new_to_old));
    new_count = 0;
    for (il = 0; il < num_global_plus_children; il++) {
      if (referenced[il]) {
        p4est_locidx_t newidx;

        newidx = new_count++;
        old_to_new[il] = newidx;
        new_to_old[newidx] = il;
      }
    }
    P4EST_ASSERT (new_count >= num_global && new_count < num_global_plus_children);
    diff = num_global_plus_children - new_count;
    num_global_plus_children -= diff;
    for (il = 0; il < num_global_plus_children; il++) {
      p4est_locidx_t oldidx = new_to_old[il];
      p4est_locidx_t *pidold, *pidnew, *cidold, *cidnew;
      int8_t *dimold, *dimnew;

      if (oldidx != il) {
        cidold = (p4est_locidx_t *) sc_array_index (child_to_id, (size_t) oldidx);
        cidnew = (p4est_locidx_t *) sc_array_index (child_to_id, (size_t) il);
        *cidnew = *cidold;
        pidold = (p4est_locidx_t *) sc_array_index (child_to_parent, (size_t) oldidx);
        pidnew = (p4est_locidx_t *) sc_array_index (child_to_parent, (size_t) il);
        *pidnew = *pidold;
        dimold = (int8_t *) sc_array_index (node_dim, (size_t) oldidx);
        dimnew = (int8_t *) sc_array_index (node_dim, (size_t) il);
        *dimnew = *dimold;
      }
    }
    sc_array_resize (child_to_id, num_global_plus_children);
    sc_array_resize (child_to_parent, num_global_plus_children);
    sc_array_resize (node_dim, num_global_plus_children);
    for (il = 0; il < K * V; il++) {
      quad_to_local[il] = old_to_new[quad_to_local[il]];
    }
    P4EST_FREE (old_to_new);
    P4EST_FREE (new_to_old);
    P4EST_FREE (referenced);
  }

  /* now that we have the list of local nodes with children,
   * we have to assign each local node a plex index: plex indices are
   * stratified, and include the cells */
  {
    int c;
    p4est_locidx_t dim_counts[4] = {0};
    p4est_locidx_t dim_offsets[5];
    p4est_locidx_t dim_cone_offsets[5];
    p4est_locidx_t *plex_to_local;
    p4est_locidx_t *local_to_plex;
    p4est_locidx_t Nplex = num_global_plus_children + K;
    int *plex_to_proc;
    p4est_locidx_t point_count, cone_count;
    p4est_gloidx_t *lnode_global_offset;

    /* count the number of points of each dimension */
    dim_counts[0] = K;
    for (il = 0; il < num_global_plus_children; il++) {
      int dim = *((int8_t *) p4est_quadrant_array_index (node_dim, (size_t) il));

      P4EST_ASSERT (0 <= dim && dim <= P4EST_DIM - 1);
      dim_counts[P4EST_DIM - dim]++;
    }
    /* get the offsets and the cone offsets of the dimensions */
    /* fill out_points_per_dim */
    sc_array_resize (out_points_per_dim, ctype_int + 1);
    dim_offsets[0] = 0;
    dim_cone_offsets[0] = 0;
    for (c = 0; c <= ctype_int; c++) {
      int dim = P4EST_DIM - c;
      p4est_locidx_t * ppd = sc_array_index (out_points_per_dim, dim);

      *ppd = dim_counts[c];
      dim_offsets[c + 1] = dim_offsets[c] + dim_counts[c];
      dim_cone_offsets[c + 1] = dim_cone_offsets[c] + 2 * dim * dim_counts[c];
    }
    P4EST_ASSERT (dim_offsets[ctype_int+1] == Nplex);
    for (c = 0; c <= ctype_int; c++) {
      int dim = P4EST_DIM - c;

      dim_offsets[c + 1] = dim_offsets[c] + dim_counts[c];
      dim_cone_offsets[c + 1] = dim_cone_offsets[c] + 2 * dim * dim_counts[c];
    }
    /* fill out_cone_sizes */
    sc_array_resize (out_cone_sizes, Nplex);
    for (c = 0; c <= ctype_int; c++) {
      p4est_locidx_t pstart, pend;
      int dim = P4EST_DIM - c;

      pstart = dim_cone_offsets[c];
      pend = dim_cone_offsets[c];
      for (il = pstart; il < pend; il++) {
        p4est_locidx_t *size = (p4est_locidx_t *) sc_array_index (out_cone_sizes, (size_t) il);

        *size = 2 * dim;
      }
    }
    /* construct the local_to_plex, plex_to_local, and plex_to_proc arrays */
    plex_to_local = P4EST_ALLOC (p4est_locidx_t, Nplex);
    local_to_plex = P4EST_ALLOC (p4est_locidx_t, Nplex);
    plex_to_proc = P4EST_ALLOC (int, Nplex);
    sc_array_resize (out_cones, dim_cone_offsets[ctype_int+1]);
    cones = (p4est_locidx_t *) out_cones->array;
    sc_array_resize (out_cone_orientations, dim_cone_offsets[ctype_int+1]);
    orientations = (p4est_locidx_t *) out_cone_orientations->array;
#ifdef P4EST_DEBUG
    memset (plex_to_local, -1, Nplex * sizeof (p4est_locidx_t));
    memset (local_to_plex, -1, Nplex * sizeof (p4est_locidx_t));
    memset (plex_to_proc, -1, Nplex * sizeof (int));
#endif
    point_count = 0;
    cone_count = 0;
    /* figure out the locations of ghosts within the base */
    Gpre = (overlap && local_first) ? 0 : (ghost->proc_offsets[mpirank] -
                                           ghost->proc_offsets[0]);
    Gpost = overlap ?
            (local_first ? G : (ghost->proc_offsets[mpisize] -
                                ghost->proc_offsets[mpirank + 1])) : 0;
    for (qid = 0, p = 0; qid < Gpre; qid++) {
      while (p < mpisize && ghost->proc_offsets[p+1] <= qid) {
        p++;
      }
      P4EST_ASSERT (ghost->proc_offsets[p] <= qid &&
                    ghost->proc_offsets[p+1] > qid);
      il = Klocal + qid;
      plex_to_local[point_count] = il;
      plex_to_proc[point_count] = p;
      local_to_plex[il] = point_count++;
    }
    for (il = 0; il < Klocal; il++) {
      plex_to_local[point_count] = il;
      plex_to_proc[point_count] = mpirank;
      local_to_plex[il] = point_count++;
    }
    for (qid = 0, p = mpirank + 1; qid < Gpost; qid++) {
      while (p < mpisize && ghost->proc_offsets[p+1] <= qid + Gpre) {
        p++;
      }
      P4EST_ASSERT (ghost->proc_offsets[p] <= qid + Gpre &&
                    ghost->proc_offsets[p+1] > qid + Gpre);
      il = Klocal + qid + Gpre;
      plex_to_local[point_count] = il;
      plex_to_proc[point_count] = p;
      local_to_plex[il] = point_count++;
    }
    /* reset dim counts so that we can us them for the offsets */
    for (c = 1; c <= ctype_int; c++) {
      dim_counts[c] = 0;
    }
    /* compute the lnodes partition */
    lnode_global_offset = P4EST_ALLOC (p4est_gloidx_t, mpisize + 1);
    lnode_global_offset[0] = 0;
    for (p = 0; p < mpisize; p++) {
      lnode_global_offset[p+1] = lnode_global_offset[p] + lnodes->global_owned_count[p];
    }
    for (il = 0, p = 0; il < num_global_plus_children; il++) {
      int dim = *((int8_t *) p4est_quadrant_array_index (node_dim, (size_t) il));
      int offset;
      p4est_locidx_t pid, nid;
      p4est_gloidx_t gid = -1;

      P4EST_ASSERT (0 <= dim && dim <= P4EST_DIM - 1);
      if (il < num_global) {
        gid = *((p4est_gloidx_t *) sc_array_index (all_global, il));
      }
      else {
        while (p < mpisize && lnode_global_offset[p+1] <= gid) {
          p++;
        }
        P4EST_ASSERT (lnode_global_offset[p] <= gid &&
                      gid < lnode_global_offset[p + 1]);
      }
      c = P4EST_DIM - dim;
      offset = dim_counts[c]++;
      pid = dim_offsets[c] + offset;
      nid = il + K;
      plex_to_local[pid] = nid;
      local_to_plex[nid] = pid;
      if (gid >= 0) {
        plex_to_proc[pid] = p;
      }
      else {
        plex_to_proc[pid] = mpirank;
      }
    }
#ifdef P4EST_DEBUG
    for (il = 0; il < Nplex; il++) {
      P4EST_ASSERT (plex_to_local[il] >= 0);
      P4EST_ASSERT (local_to_plex[il] >= 0);
      P4EST_ASSERT (plex_to_proc[il] >= 0);
    }
#endif
    /* now that we have the plex indices, we can compute out_children,
     * out_parents, out_childids */
    for (il = 0; il < Nplex; il++) {
      p4est_locidx_t nid = plex_to_local[il] - K;
      p4est_locidx_t parent;

      if (nid < 0) {
        continue;
      }
      parent = *((p4est_locidx_t *) sc_array_index (child_to_parent, nid));
      if (parent >= 0) {
        p4est_locidx_t *outc = (p4est_locidx_t *) sc_array_push (out_children);
        p4est_locidx_t *outp = (p4est_locidx_t *) sc_array_push (out_parents);
        p4est_locidx_t *outid = (p4est_locidx_t *) sc_array_push (out_childids);
        int8_t pdim = *((int8_t *) sc_array_index (node_dim, (size_t) parent));
        p4est_locidx_t id;

        *outc = il;
        *outp = local_to_plex[parent + K];
        id = *((p4est_locidx_t *) sc_array_index (child_to_id, nid));
        *outid = p4est_to_plex_child_id[P4EST_DIM - 1 - pdim][id];
      }
    }
    sc_array_destroy (child_to_parent);
    sc_array_destroy (child_to_id);

    /* compute cones and orientations */
    for (qid = 0; qid < K; il++) {
      il = plex_to_local[qid];
      for (c = 0; c < ctype_int; c++) {
        int v, vstart, vend;
        int dim = P4EST_DIM - 1 - c;

        vstart = dim_limits[c];
        vend = dim_limits[c+1];
        if (!c) {
          /* cell to face cones */
          p4est_locidx_t pid;
          p4est_locidx_t cone_off;

          pid = qid;
          cone_off = P4EST_FACES * pid;
          for (v = vstart; v < vend; v++) {
            int pos;

            P4EST_ASSERT (cones[cone_off + v] == -1);
            pos = p4est_to_plex_position[0][v - vstart];

            cones[cone_off + pos] = local_to_plex[quad_to_local[il * V + v] + K];
            orientations[cone_off + pos] = p4est_to_plex_face_orientation[v][quad_to_orientations[il * V + v]];
          }
        }
        else if (c == P4EST_DIM - 1) {
          /* x to vertex cones */
          for (v = vstart; v < vend; v++) {
            int corner = v - vstart;
            int j;
            p4est_locidx_t vid;

            vid = local_to_plex[quad_to_local[il * V + v] + K];
            for (j = 0; j < P4EST_DIM; j++) {
              p4est_locidx_t pid;
              p4est_locidx_t cone_off;
              int k, pos, o;

#ifndef P4_TO_P8
              k = p4est_corner_faces[corner][j];
              o = quad_to_orientations[il * V + k];
              pos = p4est_corner_face_corners[corner][k];
              P4EST_ASSERT (pos >= 0);
              pid = local_to_plex[quad_to_local[il * V + k] + K];
#else
              k = p8est_corner_edges[corner][j];
              o = quad_to_orientations[il * V + P4EST_FACES + k];
              if (p8est_edge_corners[k][0] == corner) {
                pos = 0;
              }
              else {
                pos = 1;
              }
              pid = local_to_plex[quad_to_local[il * V + P4EST_FACES + k] + K];
#endif
              if (o) {
                pos = (pos ^ 1);
              }
              cone_off = 2 * dim * (pid - dim_offsets[c]) + dim_cone_offsets[c];
              cone_off += pos;
              /* another cell may have already computed this cell, but we want
               * to make sure they agree */
              P4EST_ASSERT (cones[cone_off] == -1 ||
                            cones[cone_off] == vid);
              cones[cone_off] = vid;
              orientations[cone_off] = 0;
            }
          }
        }
#ifdef P4_TO_P8
        else {
          /* compute face to edge cones */
          for (v = vstart; v < vend; v++) {
            int edge = v - vstart;
            int j, o;
            p4est_locidx_t vid;

            vid = local_to_plex[quad_to_local[il * V + v] + K];
            o = quad_to_orientations[il * V + P4EST_FACES + edge];
            for (j = 0; j < 2; j++) {
              p4est_locidx_t pid;
              p4est_locidx_t cone_off;
              int k, f, pos, fo, cid[2], l, minc, maxc, or, faceor;

              k = p8est_edge_faces[edge][j];
              fo = quad_to_orientations[il * V + k];
              pid = local_to_plex[quad_to_local[il * V + k] + K];
              for (l = 0; l < 2; l++) {
                int cor, face_ex, face_in;

                cor = p8est_edge_corners[edge][l ^ o];
                face_ex = p8est_corner_face_corners[cor][k];
                P4EST_ASSERT (face_ex >= 0);
                face_in = p8est_face_permutations[fo][face_ex];
                cid[l] = face_in;
              }
              minc = SC_MIN (cid[0], cid[1]);
              maxc = SC_MAX (cid[0], cid[1]);
              /* p4est convention is to number x edges before y edges before z
               * edges, but to be consistent across dimensions, we treat edges
               * as faces, so the order of the dimensions is reversed */
              if ((maxc - minc) == 2) {
                f = (minc & 1);
              }
              else {
                f = 2 + ((minc & 2) >> 1);
              }
              pos = p4est_to_plex_position[1][f];
              cone_off = 2 * dim * (pid - dim_offsets[c]) + dim_cone_offsets[c];
              cone_off += pos;
              P4EST_ASSERT (cones[cone_off] == -1 ||
                            cones[cone_off] == vid);
              cones[cone_off] = vid;
              faceor = cid[0] < cid[1] ? 0 : 1;
              or = p4est_to_plex_edge_orientation[f][faceor];
              P4EST_ASSERT (orientations[cone_off] == -1 ||
                            orientations[cone_off] == or);
              orientations[cone_off] = or;
            }
          }
        }
#endif
      }
    }

    /* cleanup */
    P4EST_FREE (quad_to_local);
    P4EST_FREE (quad_to_orientations);
    sc_array_destroy (node_dim);

    {
      sc_array_t *quad_to_plex;
      sc_array_t *lnodes_to_plex;

      /* communicate local_to_plex for quads to compute leaves and remotes */
      if (overlap) {
        p4est_locidx_t    **mirror_data;
        sc_array_t *quad_plex;

        mirror_data = P4EST_ALLOC (p4est_locidx_t *, num_mirrors);
        quad_plex = sc_array_new_size (sizeof (p4est_locidx_t), K);
        for (il = 0; il < Klocal; il++) {
          p4est_locidx_t *ql = (p4est_locidx_t *) sc_array_index (quad_plex, il);
          *ql = local_to_plex[il];
        }
        for (il = 0; (size_t) il < num_mirrors; il++) {
          p4est_quadrant_t   *q;

          q = p4est_quadrant_array_index (&ghost->mirrors, il);
          qid = q->p.piggy3.local_num;
          mirror_data[il] = (p4est_locidx_t *) sc_array_index (quad_plex, qid);
        }
        p4est_ghost_exchange_custom (p4est, ghost,
                                     (size_t) sizeof (p4est_locidx_t),
                                     (void **) mirror_data,
                                     (p4est_locidx_t *)
                                     sc_array_index (quad_plex, Klocal));
        P4EST_FREE (mirror_data);
        for (il = 0; il < K; il++) {
          p = plex_to_proc[il];
          if (p != mpirank) {
            p4est_locidx_t ql;
            p4est_locidx_t *leaf = (p4est_locidx_t *) sc_array_push (out_leaves);
            p4est_locidx_t *remote = (p4est_locidx_t *) sc_array_push (out_remotes);

            qid = plex_to_local[il];
            ql = *((p4est_locidx_t *) sc_array_index (quad_plex, qid));
            *leaf = il;
            remote[0] = p;
            remote[1] = ql;
          }
        }
        sc_array_destroy (quad_plex);
      }
      /* communicate all_global * local_to_plex to build leaves and remotes */
      lnodes_to_plex = sc_array_new_size (sizeof (p4est_locidx_t), lnodes->num_local_nodes);
      quad_to_plex = sc_array_new_size (sizeof (p4est_locidx_t), V * K);
      if (lnodes->owned_count) {
        ssize_t firstidx;

        firstidx = sc_array_bsearch (all_global, &lnodes->global_offset, p4est_gloidx_compare);
        P4EST_ASSERT (firstidx >= 0);
        for (il = 0; il < lnodes->owned_count; il++) {
          p4est_locidx_t *lp = (p4est_locidx_t *) sc_array_index (lnodes_to_plex, (size_t) il);

          *lp = local_to_plex[firstidx + il + K];
        }
      }
      p4est_lnodes_share_owned (lnodes_to_plex, lnodes);
      for (il = 0; il < Klocal; il++) {
        for (v = 0; v < V; v++) {
          p4est_locidx_t nid = lnodes->element_nodes[il * V + v];
          p4est_locidx_t lp = *((p4est_locidx_t *) sc_array_index (lnodes_to_plex, (size_t) nid));
          p4est_locidx_t *qp = (p4est_locidx_t *) sc_array_index (quad_to_plex, (size_t) (il * V + v));

          *qp = lp;
        }
      }
      sc_array_destroy (lnodes_to_plex);
      if (overlap) {
        p4est_locidx_t    **mirror_data;

        mirror_data = P4EST_ALLOC (p4est_locidx_t *, num_mirrors);
        for (il = 0; (size_t) il < num_mirrors; il++) {
          p4est_quadrant_t   *q;

          q = p4est_quadrant_array_index (&ghost->mirrors, il);
          qid = q->p.piggy3.local_num;
          mirror_data[il] = (p4est_locidx_t *) sc_array_index (quad_to_plex, qid * V);
        }
        p4est_ghost_exchange_custom (p4est, ghost,
                                     (size_t) V * sizeof (p4est_locidx_t),
                                     (void **) mirror_data,
                                     (p4est_locidx_t *)
                                     sc_array_index (quad_to_plex, Klocal * V));
        P4EST_FREE (mirror_data);
      }
      for (il = 0; il < num_global; il++) {
        p4est_locidx_t localpid = il + K;

        p = plex_to_proc[localpid];
        if (p != mpirank) {
          p4est_gloidx_t *gid = (p4est_gloidx_t *) sc_array_index (all_global, il);
          p4est_locidx_t nid = gid[1];
          p4est_locidx_t pid = *((p4est_locidx_t *) sc_array_index (quad_to_plex, (size_t) nid));
          p4est_locidx_t *leaf = (p4est_locidx_t *) sc_array_push (out_leaves);
          p4est_locidx_t *remote = (p4est_locidx_t *) sc_array_push (out_remotes);

          *leaf = localpid;
          remote[0] = p;
          remote[1] = pid;
        }
      }
      sc_array_destroy (quad_to_plex);
    }
    P4EST_FREE (plex_to_local);
    P4EST_FREE (local_to_plex);
    P4EST_FREE (plex_to_proc);
  }

  /* cleanup */
  sc_array_destroy (all_global);
}

void
p4est_get_plex_data (p4est_t *p4est, p4est_connect_type_t ctype,
                     int overlap,
                     p4est_locidx_t *first_local_quad,
                     sc_array_t *out_points_per_dim,
                     sc_array_t *out_cone_sizes,
                     sc_array_t *out_cones,
                     sc_array_t *out_cone_orientations,
                     sc_array_t *out_vertex_coords,
                     sc_array_t *out_children,
                     sc_array_t *out_parents,
                     sc_array_t *out_childids,
                     sc_array_t *out_leaves,
                     sc_array_t *out_remotes)
{
  p4est_ghost_t *ghost;
  p4est_lnodes_t * lnodes;
  int ctype_int = p4est_connect_type_int (ctype);
  int i;

  ghost = p4est_ghost_new (p4est, ctype);
  lnodes = p4est_lnodes_new (p4est, ghost, -ctype_int);
  if (overlap) {
    p4est_ghost_support_lnodes (p4est, lnodes, ghost);
  }
  for (i = 1; i < overlap; i++) {
    p4est_ghost_expand_by_lnodes (p4est, lnodes, ghost);
  }
  if (ctype != P4EST_CONNECT_FULL) {
    p4est_lnodes_destroy (lnodes);
    lnodes = p4est_lnodes_new (p4est, ghost, -ctype);
  }
  p4est_get_plex_data_int (p4est, ghost, lnodes, overlap, 0,
                           first_local_quad, out_points_per_dim,
                           out_cone_sizes, out_cones, out_cone_orientations,
                           out_vertex_coords, out_children, out_parents,
                           out_childids, out_leaves, out_remotes);
  p4est_ghost_destroy (ghost);
  p4est_lnodes_destroy (lnodes);
}

