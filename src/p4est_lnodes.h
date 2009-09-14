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

#ifndef P4EST_LNODES_H
#define P4EST_LNODES_H

#include <p4est.h>
#include <p4est_ghost.h>

SC_EXTERN_C_BEGIN;

/** Store a parallel numbering of Lobatto points of a given degree > 0.
 *
 * Each element has degree+1 nodes per face
 * and vnodes = (degree+1)^2 nodes per volume.
 * num_local_elements is the number of local quadrants in the p4est.
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
 * Whether nodes are hanging or not is decided based on the element faces.
 * This information is encoded in face_code with one int8_t per element.
 * If no faces are hanging, the value is zero, otherwise the face_code is
 * interpreted by p4est_lnodes_decode.
 *
 * Independent nodes can be shared by multiple MPI ranks.
 * The owner rank of a node is the one from the lowest numbered element
 * on the lowest numbered octree *touching* the node.
 *
 * What is meant by *touching*?
 * A quadrant is said to touch all faces/corners that are incident on it,
 * and by extension all nodes that are contained in those faces/corners.
 * 
 *         X +-----------+
 *         o |           |
 *         o |           |
 * +-----+ o |     p     |
 * |  q  | o |           |
 * |     | o |           |
 * +-----+ O +-----------+
 * 
 * In this example degree = 6.  There are 5 nodes that live on the face
 * between q and p, and one at each corner of that face.  The face is incident
 * on q, so q owns the nodes on the face (provided q is from a lower tree or has
 * a lower index than p).  The lower corner is incident on q, so q owns it as
 * well.  The upper corner is not incident on q, so q cannot own it.
 *
 * The sharers array contains items of type p4est_lnodes_rank_t
 * that hold the ranks that own or share independent local nodes.
 * It is sorted by rank.  The rank of the current process is included.
 */
typedef struct p4est_lnodes
{
  int                 degree, vnodes;
  p4est_locidx_t      num_local_elements;
  p4est_locidx_t      num_indep_nodes;
  p4est_locidx_t      owned_offset, owned_count;
  int8_t             *face_code;
  p4est_locidx_t     *local_nodes;
  p4est_gloidx_t     *global_nodes;
  sc_array_t         *sharers;
}
p4est_lnodes_t;

/** The structure stored in the sharers array.
 *
 * shared_nodes is a sorted array of p4est_locidx_t
 * that indexes into global_nodes.  The shared_nodes array has a
 * contiguous (or empty) section of nodes owned by the current rank.
 * shared_mine_offset and shared_mine_count identify this section
 * by indexing the shared_nodes array, not the global_nodes array.
 * owned_offset and owned_count define the section of local nodes
 * that is owned by this processor (the section may be empty).
 * For the current process these coincide with those in p4est_lnodes_t.
 */
typedef struct p4est_lnodes_rank
{
  int                 rank;
  sc_array_t          shared_nodes;
  p4est_locidx_t      shared_mine_offset, shared_mine_count;
  p4est_locidx_t      owned_offset, owned_count;
}
p4est_lnodes_rank_t;

/** Decode the face_code into hanging face information.
 *
 * This is mostly for demonstration purposes.  Applications probably will
 * integrate it into their own loop over the face for performance reasons.
 *
 * \param[in] face_code as in the p4est_lnodes_t structure.
 * \param[out] hanging face: if there are hanging faces,
 *             hanging_face = -1 if the face is not hanging,
 *                          = 0 if the face is the first half,
 *                          = 1 if the face is the second half.
 *             note: not touched if there are no hanging faces.
 *           
 * \return              1 if any face is hanging, 0 otherwise.
 */
/*@unused@*/
static inline int
p4est_lnodes_decode (int8_t face_code, int hanging_face[6])
{
  SC_ASSERT (face_code >= 0);

  if (face_code) {
    int                 i;
    int8_t              c = face_code & 0x03;
    int8_t              cwork;
    int                 f;
    int8_t              work = face_code >> 2;

    memset (hanging_face, -1, 4 * sizeof (int));

    cwork = c;
    for (i = 0; i < 2; ++i) {
      f = p4est_corner_faces[c][i];
      hanging_face[f] = (work & 0x01) ? (int) (cwork & 0x01) : -1;
      cwork >>= 1;
      work >>= 1;
    }

    return 1;
  }
  else {
    return 0;
  }
}

p4est_lnodes_t     *p4est_lnodes_new (p4est_t * p4est,
                                      p4est_ghost_t * ghost_layer,
                                      int degree);

void                p4est_lnodes_destroy (p4est_lnodes_t * lnodes);

void                p4est_lnodes_share_owned (sc_array_t * array,
                                              p4est_lnodes_t * lnodes,
                                              p4est_t * p4est);

sc_array_t         *p4est_lnodes_share_all_begin (sc_array_t * array,
                                                  sc_array_t ** comm_array,
                                                  p4est_lnodes_t * lnodes,
                                                  p4est_t * p4est);

void                p4est_lnodes_share_all_end (sc_array_t * comm_array);

sc_array_t         *p4est_lnodes_share_all (sc_array_t * array,
                                            p4est_lnodes_t * lnodes,
                                            p4est_t * p4est);

SC_EXTERN_C_END;

#endif /* !P4EST_LNODES */
