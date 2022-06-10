/*
  This file is part of p4est.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2010 The University of Texas System
  Additional copyright (C) 2011 individual authors
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
#include <p4est_extended.h>
#else
#include <p8est_bits.h>
#include <p8est_lnodes.h>
#include <p8est_plex.h>
#include <p8est_extended.h>
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
  int                 ones = P4EST_CHILDREN - 1;
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

static int
p4est_gloidx_pair_compare (const void *a, const void *b)
{
  const p4est_gloidx_t *A = (const p4est_gloidx_t *) a;
  const p4est_gloidx_t *B = (const p4est_gloidx_t *) b;

  if (A[0] < B[0]) {
    return -1;
  }
  else if (A[0] == B[0]) {      /* switch the order for the second value, because we want the value selected by sc_array_uniq (the last, though this is undocumented) to be the smallest */
    if (A[1] < B[1]) {
      return 1;
    }
    else if (A[1] == B[1]) {
      return 0;
    }
    else {
      return -1;
    }
  }
  else {
    return 1;
  }
}

static void
p4est_coordinates_canonicalize (p4est_t * p4est,
                                p4est_locidx_t tree_id,
                                const p4est_qcoord_t these_coords[],
                                p4est_locidx_t * neigh_tree_id,
                                p4est_qcoord_t neigh_coords[])
{
  int                 face_code, n_outside;
  p4est_connect_type_t neigh_type;
  int                 neigh_index;
  sc_array_t          neigh_transforms;

  *neigh_tree_id = tree_id;
  face_code = 0;
  n_outside = 0;
  for (int d = 0; d < P4EST_DIM; d++) {
    neigh_coords[d] = these_coords[d];
    P4EST_ASSERT (0 <= these_coords[d] && these_coords[d] <= P4EST_ROOT_LEN);
    face_code |= (these_coords[d] == 0) ? (1 << (2 * d)) : 0;
    face_code |= (these_coords[d] == P4EST_ROOT_LEN) ? (1 << (2 * d + 1)) : 0;
    n_outside += (these_coords[d] == 0 || these_coords[d] == P4EST_ROOT_LEN);
  }
  switch (n_outside) {
  case 0:
    return;
  case 1:
    neigh_type = P4EST_CONNECT_FACE;
    neigh_index = SC_LOG2_8 (face_code);
    break;
#ifdef P4_TO_P8
  case 2:
    neigh_type = P8EST_CONNECT_EDGE;
    neigh_index = -1;
    for (int e = 0; e < P8EST_EDGES; e++) {
      if ((face_code & (1 << p8est_edge_faces[e][0])) &&
          (face_code & (1 << p8est_edge_faces[e][1]))) {
        neigh_index = e;
        break;
      }
    }
    P4EST_ASSERT (neigh_index >= 0);
    break;
#endif
  case P4EST_DIM:
    neigh_type = P4EST_CONNECT_CORNER;
    neigh_index = 0;
    for (int d = 0; d < P4EST_DIM; d++) {
      neigh_index += (face_code & (1 << (2 * d + 1))) ? (1 << d) : 0;
    }
    break;
  default:
    SC_ABORT_NOT_REACHED ();
  }
  sc_array_init (&neigh_transforms, sizeof (p4est_neighbor_transform_t));
  p4est_connectivity_get_neighbor_transforms (p4est->connectivity, tree_id,
                                              neigh_type, neigh_index,
                                              &neigh_transforms);
  /* skip the first neighbor which is self */
  for (size_t iz = 1; iz < neigh_transforms.elem_count; iz++) {
    p4est_neighbor_transform_t *nt =
      (p4est_neighbor_transform_t *) sc_array_index (&neigh_transforms, iz);
    p4est_qcoord_t      trans_coords[P4EST_DIM];

    if (nt->neighbor > *neigh_tree_id) {
      continue;
    }
    p4est_neighbor_transform_coordinates (nt, these_coords, trans_coords);
    if (nt->neighbor < *neigh_tree_id
        || p4est_coordinates_compare (trans_coords, neigh_coords) < 0) {
      *neigh_tree_id = nt->neighbor;
      for (int d = 0; d < P4EST_DIM; d++) {
        neigh_coords[d] = trans_coords[d];
      }
    }
  }
  sc_array_reset (&neigh_transforms);
}

/* Given a quadrant, local or ghost, and its

   - face code and
   - quad_to_local map,

   compute the vertex coordinates at the center of each
   local interface.  The vertex coordinates can be ambiguous
   on a periodic or mesh with disconnected vertices, so
   this function canonicalizes the results so that all
   quadrants that share a node will agree on its coordinates.

   We even compute the coordinates of mesh points that are not vertices (at
   their centers), because the mesh point may be a parent to a vertex,
   and that vertex will inherit its coordinates
*/
static void
quadrant_get_local_coordinates (p4est_t * p4est, p4est_locidx_t tree_id,
                                p4est_quadrant_t * quad,
                                p4est_lnodes_code_t face_code,
                                double coords[][3])
{
#ifndef P4_TO_P8
  int                 dim_limits[3] = { 0, 4, 8 };
#else
  int                 dim_limits[4] = { 0, 6, 18, 26 };
#endif
  int                 has_hanging, hanging[3][12];
  p4est_quadrant_t    p;
  /* initialize each local dof with (tree_id, qcoord) coordinates */
  p4est_qcoord_t      h = P4EST_QUADRANT_LEN (quad->level), H;

  P4EST_QUADRANT_INIT (&p);
#ifndef P4_TO_P8
  has_hanging = p4est_lnodes_decode (face_code, &hanging[0][0]);
#else
  has_hanging =
    p8est_lnodes_decode (face_code, &hanging[0][0], &hanging[1][0]);
#endif
  has_hanging |= lnodes_decode2 (face_code, &hanging[P4EST_DIM - 1][0]);
  H = h;
  if (has_hanging) {
    p4est_quadrant_parent (quad, &p);
    H = P4EST_QUADRANT_LEN (p.level);
  }

  /* compute the coordinates in this tree */
  for (int dim = 0; dim < P4EST_DIM; dim++) {
    int                 vstart = dim_limits[P4EST_DIM - dim - 1];
    int                 vend = dim_limits[P4EST_DIM - dim];
    const int          *vhanging = &(hanging[P4EST_DIM - dim - 1][0]);

    for (int v = 0; v < vend - vstart; v++) {
      p4est_qcoord_t      these_coords[3];
      p4est_qcoord_t      neigh_coords[3];
      p4est_quadrant_t   *q = (has_hanging && (vhanging[v] >= 0)) ? &p : quad;
      p4est_qcoord_t      vh = (has_hanging && (vhanging[v] >= 0)) ? H : h;
      p4est_locidx_t      neigh_tree_id;

      these_coords[0] = q->x;
      these_coords[1] = q->y;
#ifdef P4_TO_P8
      these_coords[2] = q->z;
#endif

      switch (dim) {
      case 0:                  /* corners */
        for (int d = 0; d < P4EST_DIM; d++) {
          these_coords[d] += (v & (1 << d)) ? vh : 0;
        }
        break;
      case (P4EST_DIM - 1):    /* faces */
        for (int d = 0; d < P4EST_DIM; d++) {
          these_coords[d] += (d == v / 2) ? ((v & 1) * vh) : vh / 2;
        }
        break;
#ifdef P4_TO_P8
      case 1:                  /* edges */
        {
          int                 lo_dir[3] = { 1, 0, 0 };
          for (int d = 0; d < P4EST_DIM; d++) {
            these_coords[d] += (d == v / 4) ? vh / 2 :
              vh * ((d == lo_dir[v / 4]) ? (v & 1) : ((v & 2) >> 1));
          }
        }
        break;
#endif
      default:
        SC_ABORT_NOT_REACHED ();
      }
      for (int d = 0; d < P4EST_DIM; d++) {
        P4EST_ASSERT (0 <= these_coords[d]
                      && these_coords[d] <= P4EST_ROOT_LEN);
      }
      p4est_coordinates_canonicalize (p4est, tree_id, these_coords,
                                      &neigh_tree_id, neigh_coords);
      p4est_qcoord_to_vertex (p4est->connectivity, neigh_tree_id,
                              neigh_coords[0], neigh_coords[1],
#ifdef P4_TO_P8
                              neigh_coords[2],
#endif
                              &coords[v + vstart][0]);
    }
  }

}

static void
mark_parent (p4est_locidx_t qid, int ctype_int, p4est_lnodes_code_t * F,
             p4est_locidx_t * quad_to_local, int8_t * is_parent,
             int8_t * node_dim)
{
#ifndef P4_TO_P8
  int                 dim_limits[3] = { 0, 4, 8 };
#else
  int                 dim_limits[4] = { 0, 6, 18, 26 };
#endif
  int                 hanging[2][12];
  int                 has_hanging;
  int                 c, V = dim_limits[ctype_int];

#ifndef P4_TO_P8
  has_hanging = p4est_lnodes_decode (F[qid], &hanging[0][0]);
#else
  has_hanging = p8est_lnodes_decode (F[qid], &hanging[0][0], &hanging[1][0]);
#endif
  if (has_hanging) {
    int                 climit;

    /* no corners */
    climit = SC_MIN (P4EST_DIM - 1, ctype_int);
    for (c = 0; c < climit; c++) {
      int                 v, vstart = dim_limits[c];
      int                 vend = dim_limits[c + 1];

      for (v = vstart; v < vend; v++) {
        if (hanging[c][v - vstart] >= 0) {
          is_parent[quad_to_local[qid * V + v]] = 1;
        }
      }
    }
  }
  for (c = 0; c < ctype_int; c++) {
    int                 v, vstart = dim_limits[c];
    int                 vend = dim_limits[c + 1];

    for (v = vstart; v < vend; v++) {
      node_dim[quad_to_local[qid * V + v]] = P4EST_DIM - 1 - c;
    }
  }
}

static void
fill_orientations (p4est_quadrant_t * q, p4est_topidx_t t,
                   p4est_connectivity_t * conn, int8_t * quad_to_orientations)
{
  int                 f;
#ifdef P4_TO_P8
  int                 e;
#endif

  for (f = 0; f < P4EST_FACES; f++) {
    p4est_quadrant_t    tempq;

    p4est_quadrant_face_neighbor (q, f, &tempq);
    quad_to_orientations[f] = 0;
    if (p4est_quadrant_is_outside_face (&tempq)) {
      p4est_topidx_t      nt;
      int                 nf, o;

      nt = conn->tree_to_tree[P4EST_FACES * t + f];
      nf = conn->tree_to_face[P4EST_FACES * t + f];
      o = nf / P4EST_FACES;
      nf = nf % P4EST_FACES;
      if (nt < t || (nt == t && nf < f)) {
        int                 set;
#ifdef P4_TO_P8
        int                 ref;
#endif

#ifndef P4_TO_P8
        set = o;
#else
        ref = p8est_face_permutation_refs[f][nf];
        set = p8est_face_permutation_sets[ref][o];
#endif
        quad_to_orientations[f] = set;
      }
    }
  }
#ifdef P4_TO_P8
  for (e = 0; e < P8EST_EDGES; e++) {
    p4est_quadrant_t    tempq;

    p8est_quadrant_edge_neighbor (q, e, &tempq);
    quad_to_orientations[P4EST_FACES + e] = 0;
    if (p4est_quadrant_is_outside_face (&tempq)) {
      int                 set;
      int                 i, f = -1;
      int                 cid[2];

      for (i = 0; i < 2; i++) {
        int                 dir;
        p4est_qcoord_t      d = 0;

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
          SC_ABORT_NOT_REACHED ();
          break;
        }
        if (d < 0 || d >= P4EST_ROOT_LEN) {
          break;
        }
      }
      P4EST_ASSERT (f >= 0);

      set = quad_to_orientations[f];
      for (i = 0; i < 2; i++) {
        int                 c, face_ex, face_in;

        c = p8est_edge_corners[e][i];
        face_ex = p8est_corner_face_corners[c][f];
        P4EST_ASSERT (face_ex >= 0);
        face_in = p8est_face_permutations[set][face_ex];
        cid[i] = face_in;
      }

      if (cid[0] < cid[1]) {
        quad_to_orientations[P4EST_FACES + e] = 0;
      }
      else {
        quad_to_orientations[P4EST_FACES + e] = 1;
      }
    }
    else if (p8est_quadrant_is_outside_edge (&tempq)) {
      int                 edge =
        conn->tree_to_edge ? conn->tree_to_edge[t * P8EST_EDGES + e] : -1;

      if (edge >= 0) {
        int                 estart, eend, i;

        estart = conn->ett_offset[edge];
        eend = conn->ett_offset[edge + 1];
        for (i = estart; i < eend; i++) {
          p4est_topidx_t      nt;
          int8_t              te;

          nt = conn->edge_to_tree[i];
          te = conn->edge_to_edge[i];
          if (nt == t && (te % P8EST_EDGES == e)) {
            quad_to_orientations[P4EST_FACES + e] = te / P8EST_EDGES;
            break;
          }
        }
        P4EST_ASSERT (i < eend);
      }
      else {
        p4est_locidx_t      ownt = t;
        int                 owne = e;
        int                 i, j, o = 0;
        for (i = 0; i < 2; i++) {
          p4est_locidx_t      nt;
          int8_t              nf;
          int                 fo;
          int                 ref, set;
          int                 cid[2];
          int                 ne;

          f = p8est_edge_faces[e][i];
          nt = conn->tree_to_tree[P4EST_FACES * t + f];
          nf = conn->tree_to_face[P4EST_FACES * t + f];
          fo = nf / P8EST_FACES;
          nf = nf % P8EST_FACES;

          ref = p8est_face_permutation_refs[f][nf];
          set = p8est_face_permutation_sets[ref][fo];
          for (j = 0; j < 2; j++) {
            int                 c, face_ex, face_in;

            c = p8est_edge_corners[e][j];
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
        quad_to_orientations[P4EST_FACES + e] = o;
      }
    }
  }
#endif
}

static void
parent_to_child (p4est_quadrant_t * q, p4est_topidx_t t, p4est_locidx_t qid,
                 int ctype_int, p4est_lnodes_code_t * F,
                 p4est_locidx_t * quad_to_local,
                 int8_t * quad_to_orientations,
                 int8_t * quad_to_orientations_orig,
                 int8_t * node_dim, p4est_locidx_t * child_offsets,
                 p4est_locidx_t * child_to_id, p4est_connectivity_t * conn,
                 int custom_numbering)
{
#ifndef P4_TO_P8
  int                 dim_limits[3] = { 0, 4, 8 };
  int                 no = P4EST_FACES;
#else
  int                 dim_limits[4] = { 0, 6, 18, 26 };
  int                 no = P4EST_FACES + P8EST_EDGES;
#endif
  int                 hanging[3][12];
  int                 has_hanging;
  int                 f, V = dim_limits[ctype_int];

#ifndef P4_TO_P8
  has_hanging = p4est_lnodes_decode (F[qid], &hanging[0][0]);
#else
  has_hanging = p8est_lnodes_decode (F[qid], &hanging[0][0], &hanging[1][0]);
#endif
  has_hanging |= lnodes_decode2 (F[qid], &hanging[P4EST_DIM - 1][0]);
  fill_orientations (q, t, conn, &quad_to_orientations[qid * no]);
  if (has_hanging) {
    int                 c, cid = p4est_quadrant_child_id (q), v;

    if (quad_to_orientations_orig) {
      p4est_quadrant_t    tempq;

      p4est_quadrant_parent (q, &tempq);
      fill_orientations (&tempq, t, conn,
                         &quad_to_orientations_orig[qid * no]);
      for (f = 0; f < P4EST_FACES; f++) {
        if (hanging[0][f] < 0) {
          quad_to_orientations_orig[qid * no + f] = -1;
        }
      }
#ifdef P4_TO_P8
      for (f = 0; f < P8EST_EDGES; f++) {
        if (hanging[1][f] < 0) {
          quad_to_orientations_orig[qid * no + P4EST_FACES + f] = -1;
        }
      }
#endif
    }
    /* no corners */
    for (c = ctype_int - 1; c >= 0; c--) {
      int                 vstart = dim_limits[c];
      int                 vend = dim_limits[c + 1];

      if (!c) {
        for (v = vstart; v < vend; v++) {
          if (hanging[0][v] >= 0) {
            int                 o = quad_to_orientations[qid * no + v];
            int                 childid = hanging[0][v];
            p4est_locidx_t      child;

#ifndef P4_TO_P8
            childid ^= o;
#else
            childid = p8est_face_permutations[o][childid];
#endif
            child = child_offsets[quad_to_local[qid * V + v]] + childid;
            quad_to_local[qid * V + v] = child;
          }
        }
      }
      else if (c == P4EST_DIM - 1) {
        for (v = vstart; v < vend; v++) {
          int                 corner = v - vstart;
          if (hanging[P4EST_DIM - 1][corner] >= 0) {
            p4est_locidx_t      child = -1;
            int                 dim = 1;

            f = p4est_child_corner_faces[cid][corner];
            P4EST_ASSERT (P4EST_DIM == 3 || f >= 0);
            if (f >= 0) {
              dim = P4EST_DIM - 1;
              child = child_offsets[quad_to_local[qid * V + f]];
            }
#ifdef P4_TO_P8
            else {
              int                 e = p8est_child_corner_edges[cid][corner];

              P4EST_ASSERT (e >= 0);
              P4EST_ASSERT (dim == 1);
              child = child_offsets[quad_to_local[qid * V + e + P4EST_FACES]];
            }
#endif
            P4EST_ASSERT (dim == 1 || dim == 2);
            child += (dim == 1) ? 2 : 8;
            quad_to_local[qid * V + v] = child;
          }
        }
      }
#ifdef P4_TO_P8
      else {
        for (v = vstart; v < vend; v++) {
          int                 edge = v - vstart;
          int                 o =
            quad_to_orientations[qid * no + P4EST_FACES + edge];

          if (hanging[1][edge] >= 0) {
            p4est_locidx_t      child;

            if (hanging[1][edge] < 4) {
              int                 h = hanging[1][edge] % 2;

              /* TODO: reconcile intrinsic/extrinsic order */
              child = child_offsets[quad_to_local[qid * V + v]] + (h ^ o);
              quad_to_local[qid * V + v] = child;
            }
            else {
              int                 i;

              for (i = 0; i < 2; i++) {
                int                 f = p8est_edge_faces[edge][i];
                int                 ch = p8est_corner_face_corners[cid][f];
                int                 e, j;

                if (ch >= 0) {
                  int                 fo = quad_to_orientations[qid * no + f];
                  int                 dir;
                  int                 hc =
                    p8est_face_permutations[fo][hanging[0][f]];
                  int                 he[2];
                  int                 child_edge;
                  int                 cid[2];
                  int                 diff;

                  P4EST_ASSERT (hanging[0][f] >= 0);

                  he[0] = (hc & 1);
                  he[1] = 2 + ((hc & 2) >> 1);

                  for (j = 0; j < 4; j++) {
                    e = p8est_face_edges[f][j];

                    if (e == edge) {
                      break;
                    }
                  }
                  P4EST_ASSERT (j < 4);
                  dir = j / 2;

                  cid[0] = p8est_face_permutations[fo][0];
                  cid[1] = p8est_face_permutations[fo][1];
                  diff = cid[1] - cid[0];
                  diff = (diff < 0) ? -diff : diff;

                  if (diff == 2) {
                    dir ^= 1;
                  }

                  if (!dir) {
                    /* first direction */
                    child_edge = he[0];
                  }
                  else {
                    /* second direction */
                    child_edge = he[1];
                  }

                  child =
                    child_offsets[quad_to_local[qid * V + f]] + P4EST_HALF +
                    child_edge;
                  quad_to_local[qid * V + v] = child;
                  if (!custom_numbering && (child_edge & 1)) {
                    quad_to_orientations[qid * no + P4EST_FACES + edge] ^= 1;
                  }
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

/* *INDENT-OFF* */
#ifndef P4_TO_P8
static const int p4est_to_plex_child_id[1][3] = {{9, 10, 25}};
static const int p4est_to_plex_face_orientation[4][2] = {{-2,  0},
                                                   { 0, -2},
                                                   { 0, -2},
                                                   {-2,  0}};
static const int p4est_to_plex_position[1][4] = {{3, 1, 0, 2}};
#else
static const int p4est_to_plex_child_id_orig[2][9] =
                                     {{15, 16, 18, 17, 90, 88, 87, 89, 137},
                                      {63, 64, 125, -1, -1, -1, -1, -1, -1}};
static const int p4est_to_plex_child_id_custom[2][9] =
                                     {{15, 16, 17, 18, 87, 88, 89, 90, 137},
                                      {63, 64, 125, -1, -1, -1, -1, -1, -1}};
static const int p4est_to_plex_face_orientation[6][8] =
                                           {{-4,  0,  1, -1, -3,  3,  2, -2},
                                            { 0, -4, -1,  1,  3, -3, -2,  2},
                                            { 0, -4, -1,  1,  3, -3, -2,  2},
                                            {-1,  3,  0, -2, -4,  2,  1, -3},
                                            {-4,  0,  1, -1, -3,  3,  2, -2},
                                            { 0, -4, -1,  1,  3, -3, -2,  2}};
static const int p4est_to_plex_edge_orientation[4][2] = {{-2,  0},
                                                   { 0, -2},
                                                   { 0, -2},
                                                   {-2,  0}};
static const int p4est_to_plex_position[2][6] = {{5, 4, 2, 3, 0, 1},
                                           {3, 1, 0, 2, -1, -1}};
#endif
/* *INDENT-ON* */

static void
p4est_get_plex_data_int (p4est_t * p4est, p4est_ghost_t * ghost,
                         p4est_lnodes_t * lnodes, int overlap,
                         int local_first, p4est_locidx_t * first_local_quad,
                         sc_array_t * out_points_per_dim,
                         sc_array_t * out_cone_sizes, sc_array_t * out_cones,
                         sc_array_t * out_cone_orientations,
                         sc_array_t * out_vertex_coords,
                         sc_array_t * out_children, sc_array_t * out_parents,
                         sc_array_t * out_childids, sc_array_t * out_leaves,
                         sc_array_t * out_remotes, int custom_numbering)
{
#ifndef P4_TO_P8
  int                 dim_limits[3] = { 0, 4, 8 };
  int                 no = P4EST_FACES;
#else
  int                 dim_limits[4] = { 0, 6, 18, 26 };
  int                 no = P4EST_FACES + P8EST_EDGES;
  const int           (*p4est_to_plex_child_id)[9] = custom_numbering ?
    p4est_to_plex_child_id_custom : p4est_to_plex_child_id_orig;
#endif
  p4est_locidx_t     *cones;
  p4est_locidx_t     *orientations;
  sc_array_t         *child_to_parent, *child_to_id;
  p4est_locidx_t     *quad_to_local, *quad_to_local_orig = NULL;
  int8_t             *quad_to_orientations, *quad_to_orientations_orig = NULL;
  p4est_lnodes_code_t *F;
  p4est_topidx_t      t, flt = p4est->first_local_tree;
  p4est_topidx_t      llt = p4est->last_local_tree;
  p4est_locidx_t      il, G = (p4est_locidx_t) ghost->ghosts.elem_count;
  p4est_locidx_t      Klocal = p4est->local_num_quadrants;
  p4est_locidx_t      K = Klocal + (overlap ? G : 0);
  p4est_locidx_t      Gpre, Gpost;
  int                 p, mpirank = p4est->mpirank;
  int                 mpisize = p4est->mpisize;
  p4est_locidx_t      qid;
  int                 v, V = lnodes->vnodes;
  sc_array_t         *is_parent, *node_dim;
  sc_array_t         *node_coords;
  int                 ctype_int = p4est_connect_type_int (P4EST_CONNECT_FULL);
  p4est_locidx_t      num_global, num_global_plus_children, last_global,
    *child_offsets;
  sc_array_t         *all_global;
  p4est_locidx_t      num_mirrors =
    (p4est_locidx_t) ghost->mirrors.elem_count;

  P4EST_ASSERT (lnodes->degree == -ctype_int);

  if (overlap) {
    /* get the face codes for ghosts */
    p4est_lnodes_code_t **mirror_data_F;

    F = P4EST_ALLOC (p4est_lnodes_code_t, K);
    memcpy (F, lnodes->face_code, Klocal * sizeof (p4est_lnodes_code_t));

    mirror_data_F = P4EST_ALLOC (p4est_lnodes_code_t *, num_mirrors);
    for (il = 0; il < num_mirrors; il++) {
      p4est_quadrant_t   *q;

      q = p4est_quadrant_array_index (&ghost->mirrors, il);
      qid = q->p.piggy3.local_num;
      mirror_data_F[il] = &F[qid];
    }
    p4est_ghost_exchange_custom (p4est, ghost,
                                 sizeof (p4est_lnodes_code_t),
                                 (void **) mirror_data_F,
                                 Klocal ? &F[Klocal] : NULL);
    P4EST_FREE (mirror_data_F);
  }
  else {
    F = lnodes->face_code;
  }

  /* we use the existing lnodes global indices as our anchors to keep all of
   * the new indices straight */
  /* create a list of all global nodes seen */
  {
    p4est_gloidx_t     *quad_to_global = P4EST_ALLOC (p4est_gloidx_t, K * V);

    /* fill the local portion of quad_to_global */
    for (qid = 0; qid < Klocal * V; qid++) {
      p4est_locidx_t      nid;

      nid = lnodes->element_nodes[qid];
      quad_to_global[qid] = p4est_lnodes_global_index (lnodes, nid);
    }
    if (overlap) {
      /* get the ghost portion of quad_to_global */
      p4est_gloidx_t    **mirror_data;

      mirror_data = P4EST_ALLOC (p4est_gloidx_t *, num_mirrors);
      for (il = 0; il < num_mirrors; il++) {
        p4est_quadrant_t   *q;

        q = p4est_quadrant_array_index (&ghost->mirrors, il);
        qid = q->p.piggy3.local_num;
        mirror_data[il] = &quad_to_global[qid * V];
      }
      p4est_ghost_exchange_custom (p4est, ghost,
                                   (size_t) V * sizeof (p4est_gloidx_t),
                                   (void **) mirror_data,
                                   Klocal ? &quad_to_global[Klocal *
                                                            V] : NULL);
      P4EST_FREE (mirror_data);
    }

    all_global = sc_array_new_size (2 * sizeof (p4est_gloidx_t), V * K);
    for (qid = 0; qid < K; qid++) {
      for (v = 0; v < V; v++) {
        p4est_gloidx_t     *pair =
          (p4est_gloidx_t *) sc_array_index (all_global,
                                             (size_t) qid * V + v);

        pair[0] = quad_to_global[qid * V + v];
        pair[1] = qid * V + v;
      }
    }
    P4EST_FREE (quad_to_global);
  }

  /* assign a local index to each global node referenced */
  sc_array_sort (all_global, p4est_gloidx_pair_compare);
  quad_to_local = P4EST_ALLOC (p4est_locidx_t, K * V);
  num_global = 0;
  last_global = -1;
  for (il = 0; il < K * V; il++) {
    p4est_gloidx_t     *pair =
      (p4est_gloidx_t *) sc_array_index (all_global, (size_t) il);
    p4est_gloidx_t      gidx;

    gidx = pair[0];
    if (gidx != last_global) {
      num_global++;
      last_global = gidx;
    }
    quad_to_local[pair[1]] = num_global - 1;
  }

  if (mpisize > 1) {
    quad_to_local_orig = P4EST_ALLOC (p4est_locidx_t, K * V);
    memcpy (quad_to_local_orig, quad_to_local,
            K * V * sizeof (p4est_locidx_t));
  }
  /* remove duplicates from all_global so we can use it to search for the
   * local index of a global index */
  sc_array_uniq (all_global, p4est_gloidx_compare);
  P4EST_ASSERT (all_global->elem_count == (size_t) num_global);

  /* in lnodes, hanging faces/edges/corners are not assigned indices, they are
   * simply references to the anchor points with a hanging type arrow.  In
   * plex, however, all points have indices, so we have to expand these points
   * at the end of the list of local nodes. */
  /* figure out which nodes are parents and mark node dimensions */
  is_parent = sc_array_new_size (sizeof (int8_t), num_global);
  node_dim = sc_array_new_size (sizeof (int8_t), num_global);
  memset (is_parent->array, 0, is_parent->elem_count * is_parent->elem_size);
  for (qid = 0; qid < K; qid++) {
    mark_parent (qid, ctype_int, F, quad_to_local,
                 (int8_t *) is_parent->array, (int8_t *) node_dim->array);
  }
  /* share the results of marking parents */
  if (mpisize > 1) {
    sc_array_t         *is_parent_lnodes;
    p4est_lnodes_buffer_t *buffer;
    size_t              zz;

    /* share equivalent lnodes */
    is_parent_lnodes =
      sc_array_new_size (sizeof (int8_t), (size_t) lnodes->num_local_nodes);
    memset (is_parent_lnodes->array, 0,
            is_parent_lnodes->elem_count * is_parent_lnodes->elem_size);
    for (il = 0; il < Klocal; il++) {
      for (v = 0; v < V; v++) {
        p4est_locidx_t      lidx = quad_to_local[il * V + v];
        p4est_locidx_t      eidx = lnodes->element_nodes[il * V + v];

        *((int8_t *) sc_array_index (is_parent_lnodes, eidx)) |=
          *((int8_t *) sc_array_index (is_parent, lidx));
      }
    }
    buffer = p4est_lnodes_share_all (is_parent_lnodes, lnodes);
    for (zz = 0; zz < lnodes->sharers->elem_count; zz++) {
      p4est_lnodes_rank_t *rank =
        p4est_lnodes_rank_array_index (lnodes->sharers, zz);
      sc_array_t         *shared_nodes = &(rank->shared_nodes);
      sc_array_t         *recv =
        (sc_array_t *) sc_array_index (buffer->recv_buffers, zz);
      size_t              zy;

      if (rank->rank == mpirank) {
        continue;
      }
      for (zy = 0; zy < shared_nodes->elem_count; zy++) {
        int8_t              val = *((int8_t *) sc_array_index (recv, zy));

        il = *((p4est_locidx_t *) sc_array_index (shared_nodes, zy));
        *((int8_t *) sc_array_index (is_parent_lnodes, il)) |= val;
      }
    }
    p4est_lnodes_buffer_destroy (buffer);
    for (il = 0; il < Klocal; il++) {
      for (v = 0; v < V; v++) {
        p4est_locidx_t      lidx = quad_to_local[il * V + v];
        p4est_locidx_t      eidx = lnodes->element_nodes[il * V + v];

        *((int8_t *) sc_array_index (is_parent, lidx)) |=
          *((int8_t *) sc_array_index (is_parent_lnodes, eidx));
      }
    }
    sc_array_destroy (is_parent_lnodes);
    if (overlap) {              /* share to ghosts */
      int8_t            **mirror_data;
      int8_t             *send_data;
      sc_array_t         *is_parent_quad;

      mirror_data = P4EST_ALLOC (int8_t *, num_mirrors);
      is_parent_quad = sc_array_new_size (V * sizeof (int8_t), K);
      for (il = 0; il < Klocal; il++) {
        int8_t             *vals =
          (int8_t *) sc_array_index (is_parent_quad, il);

        for (v = 0; v < V; v++) {
          p4est_locidx_t      lidx = quad_to_local[il * V + v];

          vals[v] = *((int8_t *) sc_array_index (is_parent, lidx));
        }
      }
      for (il = 0; il < num_mirrors; il++) {
        p4est_quadrant_t   *q;

        q = p4est_quadrant_array_index (&ghost->mirrors, il);
        qid = q->p.piggy3.local_num;
        mirror_data[il] = (int8_t *) sc_array_index (is_parent_quad, qid);
      }
      send_data =
        Klocal ? (int8_t *) sc_array_index (is_parent_quad, Klocal) : NULL;
      p4est_ghost_exchange_custom (p4est, ghost, (size_t) V * sizeof (int8_t),
                                   (void **) mirror_data, send_data);
      P4EST_FREE (mirror_data);
      for (il = Klocal; il < K; il++) {
        int8_t             *vals =
          (int8_t *) sc_array_index (is_parent_quad, il);

        for (v = 0; v < V; v++) {
          p4est_locidx_t      lidx = quad_to_local[il * V + v];

          *((int8_t *) sc_array_index (is_parent, lidx)) |= vals[v];
        }
      }
      sc_array_destroy (is_parent_quad);
    }
  }
  child_offsets = P4EST_ALLOC (p4est_locidx_t, num_global + 1);
  /* children are appended to the list of global nodes */
  child_offsets[0] = num_global;
  for (il = 0; il < num_global; il++) {
    int8_t              parent =
      *((int8_t *) sc_array_index (is_parent, (size_t) il));
    int8_t              dim =
      *((int8_t *) sc_array_index (node_dim, (size_t) il));
    int                 count;
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
    child_offsets[il + 1] = child_offsets[il] + count;
  }
  num_global_plus_children = child_offsets[num_global];

  /* expand the node_dim array so that they also have entries
   * for children */
  sc_array_resize (node_dim, num_global_plus_children);
  child_to_id =
    sc_array_new_size (sizeof (p4est_locidx_t), num_global_plus_children);
  child_to_parent =
    sc_array_new_size (sizeof (p4est_locidx_t), num_global_plus_children);
  memset (child_to_id->array, -1,
          child_to_parent->elem_count * child_to_parent->elem_size);
  memset (child_to_parent->array, -1,
          child_to_parent->elem_count * child_to_parent->elem_size);
  /* fill the child_to_parent, child_to_id, and node_dim arrays for
   * the children */
  for (il = 0; il < num_global; il++) {
    p4est_locidx_t      jstart, jend, jl;
    int8_t              pdim =
      *((int8_t *) sc_array_index (node_dim, (size_t) il));

    jstart = child_offsets[il];
    jend = child_offsets[il + 1];
    for (jl = jstart; jl < jend; jl++) {
      int8_t             *dim =
        (int8_t *) sc_array_index (node_dim, (size_t) jl);
      p4est_locidx_t     *parentidx =
        (p4est_locidx_t *) sc_array_index (child_to_parent, (size_t) jl);
      p4est_locidx_t     *childid =
        (p4est_locidx_t *) sc_array_index (child_to_id, (size_t) jl);

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
   * - compute coordinates for each node
   */
  node_coords = sc_array_new_size (3 * sizeof (double), num_global);
  for (qid = 0, t = flt; t <= llt; t++) {
    p4est_tree_t       *tree = p4est_tree_array_index (p4est->trees, t);
    sc_array_t         *quadrants = &(tree->quadrants);
    p4est_locidx_t      num_quads = (p4est_locidx_t) quadrants->elem_count;

    for (p4est_locidx_t il = 0; il < num_quads; il++, qid++) {
      p4est_quadrant_t   *q =
        p4est_quadrant_array_index (quadrants, (size_t) il);
      p4est_lnodes_code_t fc = F[qid];
      double              q_node_coords[P4EST_INSUL][3];

      quadrant_get_local_coordinates (p4est, t, q, fc, q_node_coords);
      for (int v = 0; v < V; v++) {
        double             *v_coords =
          (double *) sc_array_index (node_coords,
                                     (size_t) quad_to_local[qid * V + v]);

        for (int d = 0; d < 3; d++) {
          v_coords[d] = q_node_coords[v][d];
        }
      }
    }
  }
  if (overlap) {
    for (t = 0; t < p4est->connectivity->num_trees; t++) {
      p4est_locidx_t      il, istart = ghost->tree_offsets[t];
      p4est_locidx_t      iend = ghost->tree_offsets[t + 1];

      for (il = istart; il < iend; il++) {
        p4est_quadrant_t   *q =
          p4est_quadrant_array_index (&ghost->ghosts, (size_t) il);
        p4est_locidx_t      qid = il + Klocal;
        p4est_lnodes_code_t fc = F[qid];
        double              q_node_coords[P4EST_INSUL][3];

        quadrant_get_local_coordinates (p4est, t, q, fc, q_node_coords);
        for (int v = 0; v < V; v++) {
          double             *v_coords =
            (double *) sc_array_index (node_coords,
                                       (size_t) quad_to_local[qid * V + v]);

          for (int d = 0; d < 3; d++) {
            v_coords[d] = q_node_coords[v][d];
          }
        }
      }
    }
  }

  /* loop over quads:
   * - where quad_to_local refers to a parent, make it refer to the correct
   *   child
   * - fill quad_to_orientations
   */
  quad_to_orientations = P4EST_ALLOC (int8_t, K * no);
  if (quad_to_local_orig) {
    quad_to_orientations_orig = P4EST_ALLOC (int8_t, K * no);
    memset (quad_to_orientations_orig, -1, K * no * sizeof (int8_t));
  }
#ifdef P4EST_ENABLE_DEBUG
  memset (quad_to_orientations, -1, K * no * sizeof (int8_t));
#endif
  for (qid = 0, t = flt; t <= llt; t++) {
    p4est_tree_t       *tree = p4est_tree_array_index (p4est->trees, t);
    sc_array_t         *quadrants = &(tree->quadrants);
    p4est_locidx_t      num_quads = (p4est_locidx_t) quadrants->elem_count;

    for (il = 0; il < num_quads; il++, qid++) {
      p4est_quadrant_t   *q =
        p4est_quadrant_array_index (quadrants, (size_t) il);

      parent_to_child (q, t, qid, ctype_int, F, quad_to_local,
                       quad_to_orientations,
                       quad_to_orientations_orig,
                       (int8_t *) node_dim->array, child_offsets,
                       (p4est_locidx_t *) child_to_id->array,
                       p4est->connectivity, custom_numbering);
    }
  }
  if (overlap) {
    for (t = 0; t < p4est->connectivity->num_trees; t++) {
      p4est_locidx_t      il, istart = ghost->tree_offsets[t];
      p4est_locidx_t      iend = ghost->tree_offsets[t + 1];

      for (il = istart; il < iend; il++) {
        p4est_quadrant_t   *q =
          p4est_quadrant_array_index (&ghost->ghosts, (size_t) il);

        parent_to_child (q, t, il + Klocal, ctype_int, F, quad_to_local,
                         quad_to_orientations,
                         quad_to_orientations_orig,
                         (int8_t *) node_dim->array, child_offsets,
                         (p4est_locidx_t *) child_to_id->array,
                         p4est->connectivity, custom_numbering);
      }
    }
    P4EST_FREE (F);
  }
#ifdef P4EST_ENABLE_DEBUG
  for (il = 0; il < K * no; il++) {
    P4EST_ASSERT (quad_to_orientations[il] >= 0);
  }
#endif

  /* now that we have the list of local nodes with children,
   * we have to assign each local node a plex index: plex indices are
   * stratified, and include the cells */
  {
    int                 c;
    p4est_locidx_t      dim_counts[4] = { 0 };
    p4est_locidx_t      dim_offsets[5];
    p4est_locidx_t      dim_cone_offsets[5];
    p4est_locidx_t     *plex_to_local;
    p4est_locidx_t     *local_to_plex;
    p4est_locidx_t      Nplex = num_global_plus_children + K;
    int                *plex_to_proc;
    p4est_locidx_t      point_count;
    p4est_gloidx_t     *lnode_global_offset;
    double             *coords;

    /* count the number of points of each dimension */
    dim_counts[0] = K;
    for (il = 0; il < num_global_plus_children; il++) {
      int                 dim =
        *((int8_t *) sc_array_index (node_dim, (size_t) il));

      P4EST_ASSERT (0 <= dim && dim <= P4EST_DIM - 1);
      dim_counts[P4EST_DIM - dim]++;
    }
    /* get the offsets and the cone offsets of the dimensions */
    /* fill out_points_per_dim */
    sc_array_resize (out_points_per_dim, ctype_int + 1);
    dim_offsets[0] = 0;
    dim_cone_offsets[0] = 0;
    for (c = 0; c <= ctype_int; c++) {
      int                 dim = P4EST_DIM - c;
      p4est_locidx_t     *ppd =
        (p4est_locidx_t *) sc_array_index (out_points_per_dim, dim);

      *ppd = dim_counts[c];
      dim_offsets[c + 1] = dim_offsets[c] + dim_counts[c];
      dim_cone_offsets[c + 1] = dim_cone_offsets[c] + 2 * dim * dim_counts[c];
    }
    P4EST_ASSERT (dim_offsets[ctype_int + 1] == Nplex);
    for (c = 0; c <= ctype_int; c++) {
      int                 dim = P4EST_DIM - c;

      dim_offsets[c + 1] = dim_offsets[c] + dim_counts[c];
      dim_cone_offsets[c + 1] = dim_cone_offsets[c] + 2 * dim * dim_counts[c];
    }
    /* fill out_cone_sizes */
    sc_array_resize (out_cone_sizes, Nplex);
    for (c = 0; c <= ctype_int; c++) {
      p4est_locidx_t      pstart, pend;
      int                 dim = P4EST_DIM - c;

      pstart = dim_offsets[c];
      pend = dim_offsets[c + 1];
      for (il = pstart; il < pend; il++) {
        p4est_locidx_t     *size =
          (p4est_locidx_t *) sc_array_index (out_cone_sizes, (size_t) il);

        *size = 2 * dim;
      }
    }
    /* construct the local_to_plex, plex_to_local, and plex_to_proc arrays */
    plex_to_local = P4EST_ALLOC (p4est_locidx_t, Nplex);
    local_to_plex = P4EST_ALLOC (p4est_locidx_t, Nplex);
    plex_to_proc = P4EST_ALLOC (int, Nplex);
    sc_array_resize (out_cones, dim_cone_offsets[ctype_int + 1]);
    cones = (p4est_locidx_t *) out_cones->array;
    sc_array_resize (out_cone_orientations, dim_cone_offsets[ctype_int + 1]);
    orientations = (p4est_locidx_t *) out_cone_orientations->array;
    sc_array_resize (out_vertex_coords,
                     dim_offsets[ctype_int + 1] - dim_offsets[ctype_int]);
    coords = (double *) out_vertex_coords->array;
#ifdef P4EST_ENABLE_DEBUG
    memset (plex_to_local, -1, Nplex * sizeof (p4est_locidx_t));
    memset (local_to_plex, -1, Nplex * sizeof (p4est_locidx_t));
    memset (plex_to_proc, -1, Nplex * sizeof (int));
    memset (out_cones->array, -1,
            out_cones->elem_count * out_cones->elem_size);
    memset (out_cone_orientations->array, -1,
            out_cone_orientations->elem_count *
            out_cone_orientations->elem_size);
#endif
    point_count = 0;
    /* figure out the locations of ghosts within the base */
    Gpre =
      overlap ? (local_first ? 0
                 : (ghost->proc_offsets[mpirank] -
                    ghost->proc_offsets[0])) : 0;
    *first_local_quad = Gpre;
    Gpost = overlap ?
      (local_first ? G : (ghost->proc_offsets[mpisize] -
                          ghost->proc_offsets[mpirank + 1])) : 0;
    for (qid = 0, p = 0; qid < Gpre; qid++) {
      while (p < mpisize && ghost->proc_offsets[p + 1] <= qid) {
        p++;
      }
      P4EST_ASSERT (ghost->proc_offsets[p] <= qid &&
                    ghost->proc_offsets[p + 1] > qid);
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
      while (p < mpisize && ghost->proc_offsets[p + 1] <= qid + Gpre) {
        p++;
      }
      P4EST_ASSERT (ghost->proc_offsets[p] <= qid + Gpre &&
                    ghost->proc_offsets[p + 1] > qid + Gpre);
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
      lnode_global_offset[p + 1] =
        lnode_global_offset[p] + lnodes->global_owned_count[p];
    }
    for (il = 0, p = 0; il < num_global_plus_children; il++) {
      int                 dim =
        *((int8_t *) sc_array_index (node_dim, (size_t) il));
      int                 offset;
      p4est_locidx_t      pid, nid;
      p4est_gloidx_t      gid = -1;

      P4EST_ASSERT (0 <= dim && dim <= P4EST_DIM - 1);
      if (il < num_global) {
        gid = *((p4est_gloidx_t *) sc_array_index (all_global, il));
        while (p < mpisize && lnode_global_offset[p + 1] <= gid) {
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
    }
    for (il = 0; il < num_global; il++) {
      p4est_locidx_t      pid, nid;
      int8_t              isp;
      int                 p;

      nid = il + K;
      pid = local_to_plex[nid];
      isp = *((int8_t *) sc_array_index (is_parent, il));
      if (isp) {
        p4est_locidx_t      cstart = child_offsets[il];
        p4est_locidx_t      cend = child_offsets[il + 1], c;

        p = plex_to_proc[pid];
        for (c = cstart; c < cend; c++) {
          plex_to_proc[local_to_plex[c + K]] = p;
        }
      }
    }
    P4EST_FREE (lnode_global_offset);
#ifdef P4EST_ENABLE_DEBUG
    for (il = 0; il < Nplex; il++) {
      P4EST_ASSERT (plex_to_local[il] >= 0);
      P4EST_ASSERT (local_to_plex[il] >= 0);
      P4EST_ASSERT (plex_to_proc[il] >= 0);
    }
#endif
    /* now that we have the plex indices, we can compute out_children,
     * out_parents, out_childids */
    for (il = 0; il < Nplex; il++) {
      p4est_locidx_t      nid = plex_to_local[il] - K;
      p4est_locidx_t      parent;

      if (nid < 0) {
        continue;
      }
      parent = *((p4est_locidx_t *) sc_array_index (child_to_parent, nid));
      if (parent >= 0) {
        p4est_locidx_t     *outc =
          (p4est_locidx_t *) sc_array_push (out_children);
        p4est_locidx_t     *outp =
          (p4est_locidx_t *) sc_array_push (out_parents);
        p4est_locidx_t     *outid =
          (p4est_locidx_t *) sc_array_push (out_childids);
        int8_t              pdim =
          *((int8_t *) sc_array_index (node_dim, (size_t) parent));
        p4est_locidx_t      id;

        *outc = il;
        *outp = local_to_plex[parent + K];
        id = *((p4est_locidx_t *) sc_array_index (child_to_id, nid));
        *outid = p4est_to_plex_child_id[P4EST_DIM - 1 - pdim][id];
      }
    }

    /* compute cones and orientations */
    for (qid = 0; qid < K; qid++) {
      il = plex_to_local[qid];
      for (c = 0; c < ctype_int; c++) {
        int                 v, vstart, vend;
        int                 dim = P4EST_DIM - 1 - c;

        vstart = dim_limits[c];
        vend = dim_limits[c + 1];
        if (!c) {
          /* cell to face cones */
          p4est_locidx_t      pid;
          p4est_locidx_t      cone_off;

          pid = qid;
          cone_off = P4EST_FACES * pid;
          for (v = vstart; v < vend; v++) {
            int                 pos;

            pos = p4est_to_plex_position[0][v - vstart];
            P4EST_ASSERT (cones[cone_off + pos] == -1);

            cones[cone_off + pos] =
              local_to_plex[quad_to_local[il * V + v] + K];
            orientations[cone_off + pos] =
              p4est_to_plex_face_orientation[v][quad_to_orientations
                                                [il * no + v]];
          }
        }
        else if (c == P4EST_DIM - 1) {
          /* x to vertex cones */
          for (v = vstart; v < vend; v++) {
            int                 corner = v - vstart;
            int                 j;
            p4est_locidx_t      vid;
            int                 iter, limit =
              (quad_to_local_orig != NULL) ? 2 : 1;

            for (iter = 0; iter < limit; iter++) {
              p4est_locidx_t     *qtl =
                (iter == 0) ? quad_to_local : quad_to_local_orig;
              int8_t             *qto =
                (iter ==
                 0) ? quad_to_orientations : quad_to_orientations_orig;

              vid = local_to_plex[qtl[il * V + v] + K];
              for (j = 0; j < P4EST_DIM; j++) {
                p4est_locidx_t      pid;
                p4est_locidx_t      cone_off;
                int                 k, pos, o;

#ifndef P4_TO_P8
                k = p4est_corner_faces[corner][j];
                o = qto[il * no + k];
                pos = p4est_corner_face_corners[corner][k];
                P4EST_ASSERT (pos >= 0);
                pid = local_to_plex[qtl[il * V + k] + K];
#else
                k = p8est_corner_edges[corner][j];
                o = qto[il * no + P4EST_FACES + k];
                if (p8est_edge_corners[k][0] == corner) {
                  pos = 0;
                }
                else {
                  pos = 1;
                }
                pid = local_to_plex[qtl[il * V + P4EST_FACES + k] + K];
#endif
                if (o < 0) {
                  continue;
                }
                if (o) {
                  pos = (pos ^ 1);
                }
                P4EST_ASSERT (pid >= dim_offsets[c]
                              && pid < dim_offsets[c + 1]);
                cone_off =
                  2 * (dim + 1) * (pid - dim_offsets[c]) +
                  dim_cone_offsets[c];
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
        }
#ifdef P4_TO_P8
        else {
          /* compute face to edge cones */
          for (v = vstart; v < vend; v++) {
            int                 edge = v - vstart;
            int                 j, o;
            p4est_locidx_t      vid;
            int                 iter, limit =
              (quad_to_local_orig != NULL) ? 2 : 1;

            for (iter = 0; iter < limit; iter++) {
              p4est_locidx_t     *qtl =
                (iter == 0) ? quad_to_local : quad_to_local_orig;
              int8_t             *qto =
                (iter ==
                 0) ? quad_to_orientations : quad_to_orientations_orig;

              vid = local_to_plex[qtl[il * V + v] + K];
              o = qto[il * no + P4EST_FACES + edge];
              if (o < 0) {
                continue;
              }
              for (j = 0; j < 2; j++) {
                p4est_locidx_t      pid;
                p4est_locidx_t      cone_off;
                int                 k, f, pos, fo, cid[2], l, minc, maxc;
                int                 edgeor, faceor;

                k = p8est_edge_faces[edge][j];
                fo = qto[il * no + k];
                if (fo < 0) {
                  continue;
                }
                pid = local_to_plex[qtl[il * V + k] + K];
                for (l = 0; l < 2; l++) {
                  int                 cor, face_ex, face_in;

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
                P4EST_ASSERT (pid >= dim_offsets[c]
                              && pid < dim_offsets[c + 1]);
                cone_off =
                  2 * (dim + 1) * (pid - dim_offsets[c]) +
                  dim_cone_offsets[c];
                cone_off += pos;
                P4EST_ASSERT (cones[cone_off] == -1 ||
                              cones[cone_off] == vid);
                cones[cone_off] = vid;
                faceor = cid[0] < cid[1] ? 0 : 1;
                edgeor = p4est_to_plex_edge_orientation[f][faceor];
                P4EST_ASSERT (orientations[cone_off] == -1 ||
                              orientations[cone_off] == edgeor);
                orientations[cone_off] = edgeor;
              }
            }
          }
        }
#endif
      }
    }
    /* correct for floating children */
    for (il = 0; il < num_global; il++) {
      p4est_locidx_t      cstart, cend, ppid;
      p4est_locidx_t     *pcones, poff;
#ifdef P4_TO_P8
      p4est_locidx_t     *ornts;
#endif
      int8_t              pdim;

      cstart = child_offsets[il];
      cend = child_offsets[il + 1];
      if (cstart == cend) {
        continue;
      }
      pdim = *((int8_t *) sc_array_index (node_dim, (size_t) il));
      ppid = local_to_plex[K + il];
      poff =
        dim_cone_offsets[P4EST_DIM - pdim] + 2 * pdim * (ppid -
                                                         dim_offsets[P4EST_DIM
                                                                     - pdim]);
      pcones = &cones[poff];
#ifdef P4_TO_P8
      ornts = &orientations[poff];
#endif
      for (c = cstart; c < cend - 1; c++) {
        p4est_locidx_t      cpid;
        int8_t              cdim;
        p4est_locidx_t     *ccones, coff;
        p4est_locidx_t     *cornts;

        cpid = local_to_plex[K + c];
        cdim = *((int8_t *) sc_array_index (node_dim, (size_t) c));
        coff =
          dim_cone_offsets[P4EST_DIM - cdim] + 2 * cdim * (cpid -
                                                           dim_offsets
                                                           [P4EST_DIM -
                                                            cdim]);
        ccones = &cones[coff];
        cornts = &orientations[coff];
        if (cdim == 1) {
          int                 side = (c - cstart) & 1;

          cornts[0] = cornts[1] = 0;
          if (pdim == 1) {
            ccones[1 - side] = local_to_plex[K + cend - 1];
            ccones[side] = pcones[side];
          }
#ifdef P4_TO_P8
          else {
            p4est_locidx_t      epid, lepid, eend;

            epid = pcones[p4est_to_plex_position[1][c - (cstart + 4)]];
            lepid = plex_to_local[epid] - K;

            eend = child_offsets[lepid + 1];
            P4EST_ASSERT (eend > child_offsets[lepid]);
            ccones[custom_numbering ? (1 - side) : 1] =
              local_to_plex[K + cend - 1];
            ccones[custom_numbering ? side : 0] = local_to_plex[K + eend - 1];
          }
#endif
        }
#ifdef P4_TO_P8
        else {
          int                 j;
          p4est_locidx_t      cone_to_child[4][4] =
            { {-1, 6, 4, -4}, {-1, -2, 5, 6}, {4, 7, -3, -4}, {5, -2, -3,
                                                               7}
          };
          int                 cone_to_side[4][4] =
            { {0, -1, -1, 1}, {1, 0, -1, -1}, {-1, -1, 1, 0}, {-1, 1, 0,
                                                               -1}
          };

          for (j = 0; j < 4; j++) {
            p4est_locidx_t      nchild = cone_to_child[c - cstart][j];
            if (nchild >= 0) {
              p4est_locidx_t      prevchild =
                cone_to_child[c - cstart][(j + 3) % 4];

              ccones[j] = local_to_plex[K + cstart + nchild];
              cornts[j] = custom_numbering ? ((j < 2) ? 0 : -2)
                : ((prevchild < 0) ? 0 : -2);
            }
            else {
              int                 epid, lepid, estart;
              int                 side = cone_to_side[c - cstart][j];

              epid = pcones[-(nchild + 1)];
              lepid = plex_to_local[epid] - K;

              estart = child_offsets[lepid];
              P4EST_ASSERT (child_offsets[lepid + 1] > estart);
              cornts[j] = ornts[j];
              if (!ornts[-(nchild + 1)]) {
                ccones[j] = local_to_plex[K + estart + side];
              }
              else {
                ccones[j] = local_to_plex[K + estart + (1 - side)];
              }
            }
          }
        }
#endif
      }

    }
#ifdef P4EST_ENABLE_DEBUG
    {
      size_t              zz, count = out_cones->elem_count;

      for (zz = 0; zz < count; zz++) {
        p4est_locidx_t      cone =
          *((p4est_locidx_t *) sc_array_index (out_cones, zz));

        P4EST_ASSERT (cone >= 0);
      }
    }
#endif
    /* compute coordinates */
    for (p4est_locidx_t v = 0;
         v < dim_offsets[ctype_int + 1] - dim_offsets[ctype_int]; v++) {
      p4est_locidx_t      lid = plex_to_local[v + dim_offsets[ctype_int]] - K;
      p4est_locidx_t      parent;
      const double       *n_coords;
      double             *v_coords = &coords[v * 3];

      P4EST_ASSERT (lid >= 0);

      parent = *((p4est_locidx_t *) sc_array_index (child_to_parent, lid));
      if (parent >= 0) {
        n_coords = (double *) sc_array_index (node_coords, parent);
      }
      else {
        n_coords = (double *) sc_array_index (node_coords, lid);
      }
      for (int d = 0; d < 3; d++) {
        v_coords[d] = n_coords[d];
      }
    }

    sc_array_destroy (node_coords);

    /* cleanup */
    P4EST_FREE (quad_to_local);
    P4EST_FREE (quad_to_orientations);
    if (quad_to_local_orig) {
      P4EST_FREE (quad_to_local_orig);
      P4EST_FREE (quad_to_orientations_orig);
    }

    {
      sc_array_t         *quad_to_plex;
      sc_array_t         *lnodes_to_plex;

      /* communicate local_to_plex for quads to compute leaves and remotes */
      if (overlap) {
        p4est_locidx_t    **mirror_data;
        p4est_locidx_t     *send_data;
        sc_array_t         *quad_plex;

        mirror_data = P4EST_ALLOC (p4est_locidx_t *, num_mirrors);
        quad_plex = sc_array_new_size (sizeof (p4est_locidx_t), K);
        for (il = 0; il < Klocal; il++) {
          p4est_locidx_t     *ql =
            (p4est_locidx_t *) sc_array_index (quad_plex, il);
          *ql = local_to_plex[il];
        }
        for (il = 0; il < num_mirrors; il++) {
          p4est_quadrant_t   *q;

          q = p4est_quadrant_array_index (&ghost->mirrors, il);
          qid = q->p.piggy3.local_num;
          mirror_data[il] =
            (p4est_locidx_t *) sc_array_index (quad_plex, qid);
        }
        send_data = Klocal ? (p4est_locidx_t *) sc_array_index (quad_plex,
                                                                Klocal) :
          NULL;
        p4est_ghost_exchange_custom (p4est, ghost,
                                     (size_t) sizeof (p4est_locidx_t),
                                     (void **) mirror_data, send_data);
        P4EST_FREE (mirror_data);
        for (il = 0; il < K; il++) {
          p = plex_to_proc[il];
          if (p != mpirank) {
            p4est_locidx_t      ql;
            p4est_locidx_t     *leaf =
              (p4est_locidx_t *) sc_array_push (out_leaves);
            p4est_locidx_t     *remote =
              (p4est_locidx_t *) sc_array_push (out_remotes);

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
      lnodes_to_plex =
        sc_array_new_size (sizeof (p4est_locidx_t) * (P4EST_DIM + 1),
                           lnodes->num_local_nodes);
      quad_to_plex =
        sc_array_new_size (sizeof (p4est_locidx_t) * (P4EST_DIM + 1), V * K);
      if (lnodes->owned_count) {
        ssize_t             firstidx;

        firstidx =
          sc_array_bsearch (all_global, &lnodes->global_offset,
                            p4est_gloidx_compare);
        P4EST_ASSERT (firstidx >= 0);
        for (il = 0; il < lnodes->owned_count; il++) {
          p4est_gloidx_t     *gid =
            (p4est_gloidx_t *) sc_array_index (all_global,
                                               (size_t) firstidx + il);
          p4est_locidx_t      eid = (p4est_locidx_t) gid[1];
          p4est_locidx_t      lid = lnodes->element_nodes[eid];
          p4est_locidx_t      loc, cstart, cend, j;
          p4est_locidx_t     *lp;

          P4EST_ASSERT (gid[0] == lnodes->global_offset + il);
          P4EST_ASSERT (eid / V < lnodes->num_local_elements);
          lp =
            (p4est_locidx_t *) sc_array_index (lnodes_to_plex, (size_t) lid);

          loc = firstidx + il;
          lp[0] = local_to_plex[loc + K];
          for (j = 1; j < P4EST_DIM + 1; j++) {
            lp[j] = -1;
          }
          if (loc < num_global) {
            cstart = child_offsets[loc];
            cend = child_offsets[loc + 1];
            if (cend > cstart) {
              int8_t              ndim =
                *((int8_t *) sc_array_index (node_dim, loc)), d;

              for (d = ndim; d >= 0; d--) {
                lp[P4EST_DIM - d] = local_to_plex[cstart + K];
                if (d == ndim) {
                  cstart += 2 * ndim;
                }
#ifdef P4_TO_P8
                else if (d == 1) {
                  cstart += 4;
                }
#endif
              }
            }
          }
        }
      }
      p4est_lnodes_share_owned (lnodes_to_plex, lnodes);
      for (il = 0; il < Klocal; il++) {
        for (v = 0; v < V; v++) {
          int                 j;
          p4est_locidx_t      nid = lnodes->element_nodes[il * V + v];
          p4est_locidx_t     *lp =
            (p4est_locidx_t *) sc_array_index (lnodes_to_plex, (size_t) nid);
          p4est_locidx_t     *qp =
            (p4est_locidx_t *) sc_array_index (quad_to_plex,
                                               (size_t) (il * V + v));

          for (j = 0; j < P4EST_DIM + 1; j++) {
            qp[j] = lp[j];
          }
        }
      }
      sc_array_destroy (lnodes_to_plex);
      if (overlap) {
        p4est_locidx_t    **mirror_data;
        p4est_locidx_t     *send_data;

        mirror_data = P4EST_ALLOC (p4est_locidx_t *, num_mirrors);
        for (il = 0; il < num_mirrors; il++) {
          p4est_quadrant_t   *q;

          q = p4est_quadrant_array_index (&ghost->mirrors, il);
          qid = q->p.piggy3.local_num;
          mirror_data[il] =
            (p4est_locidx_t *) sc_array_index (quad_to_plex, qid * V);
        }
        send_data = Klocal ? (p4est_locidx_t *) sc_array_index (quad_to_plex,
                                                                Klocal * V) :
          NULL;
        p4est_ghost_exchange_custom (p4est, ghost,
                                     (size_t) V * (P4EST_DIM +
                                                   1) *
                                     sizeof (p4est_locidx_t),
                                     (void **) mirror_data, send_data);
        P4EST_FREE (mirror_data);
      }
      for (il = 0; il < num_global_plus_children; il++) {
        p4est_locidx_t      localpid = il + K;

        p = plex_to_proc[localpid];
        if (p != mpirank) {
          p4est_locidx_t      lid = plex_to_local[localpid] - K;

          if (lid < num_global) {
            p4est_gloidx_t     *gid =
              (p4est_gloidx_t *) sc_array_index (all_global, lid);
            p4est_locidx_t      eid = gid[1];
            p4est_locidx_t      pid = *((p4est_locidx_t *)
                                        sc_array_index (quad_to_plex,
                                                        (size_t) eid));
            p4est_locidx_t     *leaf =
              (p4est_locidx_t *) sc_array_push (out_leaves);
            p4est_locidx_t     *remote =
              (p4est_locidx_t *) sc_array_push (out_remotes);

            *leaf = localpid;
            remote[0] = p;
            remote[1] = pid;
          }
          else {
            p4est_locidx_t      parent =
              *((p4est_locidx_t *) sc_array_index (child_to_parent, lid));
            p4est_locidx_t      id =
              *((p4est_locidx_t *) sc_array_index (child_to_id, lid));
            int8_t              dim =
              *((int8_t *) sc_array_index (node_dim, lid));

            P4EST_ASSERT (parent >= 0);
            {
              p4est_gloidx_t     *pgid =
                (p4est_gloidx_t *) sc_array_index (all_global, parent);
              p4est_locidx_t      peid = pgid[1];
              p4est_locidx_t     *ppid =
                (p4est_locidx_t *) sc_array_index (quad_to_plex,
                                                   (size_t) peid);
              p4est_locidx_t     *leaf =
                (p4est_locidx_t *) sc_array_push (out_leaves);
              p4est_locidx_t     *remote =
                (p4est_locidx_t *) sc_array_push (out_remotes);

              *leaf = localpid;
              remote[0] = p;
              remote[1] = -1;
              if (dim == 0) {
                remote[1] = ppid[P4EST_DIM - dim];
              }
              else if (dim == P4EST_DIM - 1) {
                remote[1] = ppid[P4EST_DIM - dim] + id;
              }
#ifdef P4_TO_P8
              else {
                if (ppid[1] == -1) {    /* parent is an edge */
                  remote[1] = ppid[P4EST_DIM - dim] + id;
                }
                else {          /* parent is a facet */
                  remote[1] = ppid[P4EST_DIM - dim] + (id - 4);
                }
              }
#endif
            }
          }
        }
      }
      sc_array_destroy (quad_to_plex);
    }
    sc_array_destroy (child_to_parent);
    sc_array_destroy (child_to_id);
    sc_array_destroy (node_dim);
    P4EST_FREE (plex_to_local);
    P4EST_FREE (local_to_plex);
    P4EST_FREE (plex_to_proc);
  }
  /* cleanup */
  sc_array_destroy (is_parent);
  P4EST_FREE (child_offsets);

  sc_array_destroy (all_global);
}

void
p4est_get_plex_data_ext (p4est_t * p4est,
                         p4est_ghost_t ** ghost,
                         p4est_lnodes_t ** lnodes,
                         p4est_connect_type_t ctype,
                         int overlap, p4est_locidx_t * first_local_quad,
                         sc_array_t * out_points_per_dim,
                         sc_array_t * out_cone_sizes,
                         sc_array_t * out_cones,
                         sc_array_t * out_cone_orientations,
                         sc_array_t * out_vertex_coords,
                         sc_array_t * out_children,
                         sc_array_t * out_parents,
                         sc_array_t * out_childids,
                         sc_array_t * out_leaves, sc_array_t * out_remotes,
                         int custom_numbering)
{
  int                 ctype_int = p4est_connect_type_int (ctype);
  int                 i;
  int                 created_ghost = 0;

  if (!*ghost) {
    *ghost = p4est_ghost_new (p4est, ctype);
    created_ghost = 1;
  }
  if (!*lnodes) {
    *lnodes = p4est_lnodes_new (p4est, *ghost, -ctype_int);
  }
  if (created_ghost) {
    if (overlap) {
      p4est_ghost_support_lnodes (p4est, *lnodes, *ghost);
    }
    for (i = 1; i < overlap; i++) {
      p4est_ghost_expand_by_lnodes (p4est, *lnodes, *ghost);
    }
  }
  if (ctype != P4EST_CONNECT_FULL) {
    p4est_lnodes_destroy (*lnodes);
    *lnodes = p4est_lnodes_new (p4est, *ghost, -ctype);
  }
  p4est_get_plex_data_int (p4est, *ghost, *lnodes, overlap, 0,
                           first_local_quad, out_points_per_dim,
                           out_cone_sizes, out_cones, out_cone_orientations,
                           out_vertex_coords, out_children, out_parents,
                           out_childids, out_leaves, out_remotes,
                           custom_numbering);
}

void
p4est_get_plex_data (p4est_t * p4est, p4est_connect_type_t ctype,
                     int overlap, p4est_locidx_t * first_local_quad,
                     sc_array_t * out_points_per_dim,
                     sc_array_t * out_cone_sizes,
                     sc_array_t * out_cones,
                     sc_array_t * out_cone_orientations,
                     sc_array_t * out_vertex_coords,
                     sc_array_t * out_children,
                     sc_array_t * out_parents,
                     sc_array_t * out_childids,
                     sc_array_t * out_leaves, sc_array_t * out_remotes)
{
  p4est_ghost_t      *ghost = NULL;
  p4est_lnodes_t     *lnodes = NULL;

  p4est_get_plex_data_ext (p4est, &ghost, &lnodes, ctype, overlap,
                           first_local_quad, out_points_per_dim,
                           out_cone_sizes, out_cones,
                           out_cone_orientations, out_vertex_coords,
                           out_children, out_parents, out_childids,
                           out_leaves, out_remotes, 0);

  p4est_lnodes_destroy (lnodes);
  p4est_ghost_destroy (ghost);
}
