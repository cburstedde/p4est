/*
  This file is part of p4est.
  p4est is a C library to manage a parallel collection of quadtrees and/or
  octrees.

  Copyright (C) 2009 Carsten Burstedde, Lucas Wilcox.

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

#ifndef P8EST_LNODES_H
#define P8EST_LNODES_H

#include <p8est.h>

SC_EXTERN_C_BEGIN;

/** Store a parallel numbering of Lobatto points of a given degree > 0.
 *
 * Each element has degree+1 nodes per edge
 * and vnodes = (degree+1)^3 nodes per volume.
 * num_local_elements is the number of local quadrants in the p8est.
 * local_nodes is of dimension vnodes * num_local_elements
 * and indexes into the array global_nodes layed out as follows:
 * global_nodes = [<--------------->|<-------------------->|          ]
 *                  \ owned_offset    \ owned_count
 *                 <------------------------------------------------->
 *                  \ num_indep_nodes
 * global_nodes contains the globally unique numbers for independent nodes.
 * Hanging nodes are always local and don't have a global number.
 * They index the geometrically corresponding global indep_node of a neighbor.
 *
 * Whether nodes are hanging or not is decided based on the element faces and
 * edges. This information is encoded in face_code with one int16_t per
 * element. If no faces are hanging, the value is zero, otherwise the face_code
 * is interpreted by p8est_lnodes_decode.
 *
 * Independent nodes can be shared by multiple MPI ranks.
 * The owner rank of a node is the one from the lowest numbered element
 * on the lowest numbered octree touching the node.
 * The sharers array contains items of type p8est_lnodes_rank_t
 * that hold the ranks that own or share independent local nodes.
 * It is sorted by rank.  The rank of the current process is included.
 */
typedef struct p8est_lnodes
{
  int                 degree, vnodes;
  p4est_locidx_t      num_local_elements;
  p4est_locidx_t      num_indep_nodes;
  p4est_locidx_t      owned_offset, owned_count;
  int16_t            *face_code;
  p4est_locidx_t     *local_nodes;
  p4est_gloidx_t     *global_nodes;
  sc_array_t         *sharers;
}
p8est_lnodes_t;

/** The structure stored in the sharers array.
 *
 * shared_nodes is a sorted array of p4est_locidx_t
 * that indexes into global_nodes.  The shared_nodes array has a
 * contiguous (or empty) section of nodes owned by the current rank.
 * shared_mine_offset and shared_mine_count identify this section
 * by indexing the shared_nodes array, not the global_nodes array.
 * owned_offset and owned_count define the section of local nodes
 * that is owned by this processor (the section may be empty).
 * For the current process these coincide with those in p8est_lnodes_t.
 */
typedef struct p8est_lnodes_rank
{
  int                 rank;
  sc_array_t          shared_nodes;
  p4est_locidx_t      shared_mine_offset, shared_mine_count;
  p4est_locidx_t      owned_offset, owned_count;
}
p8est_lnodes_rank_t;

/** Decode the face_code into hanging face information.
 *
 * This is mostly for demonstration purposes.  Applications probably will
 * integrate it into their own loop over the face for performance reasons.
 *
 * \param[in] face_code as in the p8est_lnodes_t structure.
 * \param[out] hanging_face: if there are hanging faces or edges,
 *             hanging_face = -1 if the face is not hanging,
 *                          = the corner of the full face that it touches:
 *                            e.g. if face = i and hanging_face[i] = 
 *                            j, then the interpolation operator corresponding
 *                            to corner j should be used for that face.
 *             note: not touched if there are no hanging faces or edges.
 * \param[out] hanging_edge: if there are hanging faces or edges,
 *             hanging_edge = -1 if the edge is not hanging or touches a
 *                               hanging face,
 *                          = 0 if the edge is the first half of the full edge,
 *                          = 1 if the edge is the second half.
 *             not: not touched if there are no hanging faces or edges;
 * \return              true if any face or edge is hanging, false otherwise.
 */
/*@unused@*/
static inline       bool
p8est_lnodes_decode (int16_t face_code, int hanging_face[6],
                     int hanging_edge[12])
{
  SC_ASSERT (face_code >= 0);

  if (face_code) {
    int                 i;
    int16_t             c = face_code & 0x0007;
    int16_t             cwork;
    int                 f;
    int                 e;
    int16_t             work = face_code >> 3;

    memset (hanging_face, -1, 6 * sizeof (int));
    memset (hanging_edge, -1, 12 * sizeof (int));

    cwork = c;
    for (i = 0; i < 3; ++i) {
      e = p8est_corner_edges[c][i];
      hanging_edge[e] = (work & 0x0001) ? (int) (cwork & 0x0001) : -1;
      cwork >>= 1;
      work >>= 1;
    }
    for (i = 0; i < 3; ++i) {
      f = p8est_corner_faces[c][i];
      hanging_face[f] =
        (work & 0x0001) ? p8est_corner_face_corners[c][f] : -1;
      work >>= 1;
    }

    return true;
  }
  else {
    return false;
  }
}

p8est_lnodes_t     *p8est_lnodes_new (p8est_t * p8est,
                                      sc_array_t * ghost_layer, int degree);
void                p8est_lnodes_destroy (p8est_lnodes_t *);

SC_EXTERN_C_END;

#endif /* !P8EST_LNODES */
