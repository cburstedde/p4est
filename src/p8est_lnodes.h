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
#include <p8est_ghost.h>

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
 * on the lowest numbered octree *touching* the node.
 *
 * What is meant by *touching*?
 * A quadrant is said to touch all faces/edges/corners that are incident on it,
 * and by extension all nodes that are contained in those faces/edges/corners.
 *
 *            X      +-----------+
 *             x     |\           \
 *            x      | \           \
 *             . x   |  \           \
 *            x   X  |   +-----------+
 * +-----+     . .   |   |           |
 * |\     \   X   o  +   |           |
 * | +-----+   o .    \  |     p     |
 * + |  q  |      o    \ |           |
 *  \|     |     o      \|           |
 *   +-----+      O      +-----------+
 * 
 * In this example degree = 3.  There are 4 nodes that live on the face
 * between q and p, two on each edge and one at each corner of that face.
 * The face is incident on q, so q owns the nodes marked '.' on the face
 * (provided q is from a lower tree or has a lower index than p).
 * The bottom and frone edges are incident on q, so q owns its nodes marked
 * 'o' as well.
 * The front lower corner is incident on q, so q owns its node 'O' as
 * well.  The other edges and Corner are not incident on q, so q cannot own
 * their nodes, marked 'x' and 'X'.
 * 
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
 * \return             true if any face or edge is hanging, false otherwise.
 */
/*@unused@*/
static inline int
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

    return 1;
  }
  else {
    return 0;
  }
}

p8est_lnodes_t     *p8est_lnodes_new (p8est_t * p8est,
                                      p8est_ghost_t * ghost_layer,
                                      int degree);

void                p8est_lnodes_destroy (p8est_lnodes_t *);

/** p8est_lnodes_buffer_t handles the communication of data associated with
 * nodes.
 *
 * \a send_buffers is an array of arrays: one buffer for each process to which
 * the current process sends node-data.  It should not be altered between
 * a shared_*_begin and a shared_*_end call.
 *
 * \a recv_buffers is an array of arrays that is used in lnodes_share_all_*.
 * \a recv_buffers[j] corresponds with lnodes->sharers[j]: it is the same
 * length as \a lnodes->sharers[j]->shared_nodes.  At the completion of
 * lnodes_share_all or lnodes_share_all_end, recv_buffers[j] contains the
 * node-data from the process lnodes->sharers[j]->rank
 * (unless j is the current rank, in which case recv_buffers[j] is empty).
 */
typedef struct p8est_lnodes_buffer
{
  sc_array_t         *requests; /* MPI_Request */
  sc_array_t         *send_buffers;
  sc_array_t         *recv_buffers;
}
p8est_lnodes_buffer_t;

/** p8est_lnodes_share_owned_begin
 *
 * \a node_data is a user-defined array of arbitrary type, where each entry
 * is associated with the \a lnodes->global_nodes entry of matching index.
 * For every \a lnodes->global_nodes entry that is owned by a process
 * other than the current one, the value in the \a node_data array of the
 * owning process is written directly into the \a node_data array of the current
 * process.  Values of \a node_data are not guaranteed to be sent or received
 * until the \a buffer created by p8est_lnodes_share_owned_begin is passed to
 * p8est_lnodes_share_owned_end.
 *
 * To be memory neutral, the \a buffer created by
 * p8est_lnodes_share_owned_begin must be destroying with
 * p8est_lnodes_buffer_destroy (it is not destroyed by
 * p8est_lnodes_share_owned_end).
 */
p8est_lnodes_buffer_t *p8est_lnodes_share_owned_begin (sc_array_t * node_data,
                                                       p8est_lnodes_t *
                                                       lnodes,
                                                       p8est_t * p8est);

void                p8est_lnodes_share_owned_end (p8est_lnodes_buffer_t *
                                                  buffer);

/** Equivalent to calling p8est_lnodes_share_owned_end directly after
 * p8est_lnodes_share_owned_begin.  Use if there is no local work that can be
 * done to mask the communication cost.
 */
void                p8est_lnodes_share_owned (sc_array_t * node_data,
                                              p8est_lnodes_t * lnodes,
                                              p8est_t * p8est);

/** p8est_lnodes_share_all_begin
 *
 * \a node_data is a user_defined array of arbitrary type, where each entry
 * is associated with the \a lnodes->global_nodes entry of matching index.
 * For every process that shares an entry with the current one, the value in
 * the \a node_data array of that process is written into a
 * \a buffer->recv_buffers entry as described above.  The user can then perform
 * some arbitrary work that requires the data from all processes that share a
 * node (such as reduce, max, min, etc.).  When the work concludes, the
 * \a buffer should be destroyed with p8est_lnodes_buffer_destroy.
 *
 * Values of \a node_data are not guaranteed to be send, and
 * \a buffer->recv_buffer entries are not guaranteed to be received until
 * the \a buffer created by p8est_lnodes_share_all_begin is passed to
 * p8est_lnodes_share_all_end.
 */
p8est_lnodes_buffer_t *p8est_lnodes_share_all_begin (sc_array_t * node_data,
                                                     p8est_lnodes_t * lnodes,
                                                     p8est_t * p8est);

void                p8est_lnodes_share_all_end (p8est_lnodes_buffer_t *
                                                buffer);

/** Equivalend to calling p8est_lnodes_share_all_end directly after
 * p8est_lnodes_share_all_begin.  Use if there is no local work that can be
 * done to maks the communication cost.
 */
p8est_lnodes_buffer_t *p8est_lnodes_share_all (sc_array_t * node_data,
                                               p8est_lnodes_t * lnodes,
                                               p8est_t * p8est);

void                p8est_lnodes_buffer_destroy (p8est_lnodes_buffer_t *
                                                 buffer);

SC_EXTERN_C_END;

#endif /* !P8EST_LNODES */
