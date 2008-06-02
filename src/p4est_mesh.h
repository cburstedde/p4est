/*
  This file is part of p4est.
  p4est is a C library to manage a parallel collection of quadtrees and/or
  octrees.

  Copyright (C) 2007,2008 Carsten Burstedde, Lucas Wilcox.

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

#ifndef P4EST_MESH_H
#define P4EST_MESH_H

#include <p4est.h>

/** Determine unique ordering of vertices for each quadrant.
 *
 * \param [in]  p4est              The forest whose vertices will be ordered.
 * \param [in]  identify_periodic  Boolean flag to switch on periodic b.c.
 * \param [out] num_uniq_local_vertices will be filled with the total number
 *                                      of unique vertices.
 * \param [out] quadrant_to_local_vertex an array that for each vertex of each
 *                                       quadrant holds the value of the
 *                                       corresponding unique local vertex.
 *                                       The array needs to be size (total
 *                                       number of local quadrants x 4).  It
 *                                       will be filled
 *                                       [0][0]..[0][3]..
 *                                       [num_local_quads-1][0]..
 *                                       [num_local_quads-1][3]
 */
void                p4est_order_local_vertices (p4est_t * p4est,
                                                bool identify_periodic,
                                                p4est_topidx_t *
                                                num_uniq_local_vertices,
                                                p4est_topidx_t *
                                                quadrant_to_local_vertex);

/** Generate a neighbor of a quadrant for a given node.
 *
 * The neighbor numbering is given below.
 *
 * Neighbor numbering for q, node=0:
 *
 *      ------+------+
 *            |  q   |
 *      nnum=2|nnum=3|
 *            |      |
 *      ------+------+
 *            |      |
 *      nnum=0|nnum=1|
 *
 * Neighbor numbering for q, node=1:
 *
 *            +------+------
 *            |  q   |
 *            |nnum=2|num=3
 *            |      |
 *            +------+------
 *            |      |
 *            |nnum=0|nnum=1
 *
 * Neighbor numbering for q, node=2:
 *
 *            |      |
 *      nnum=2|nnum=3|
 *            |      |
 *      ------+------+
 *            |  q   |
 *      nnum=0|nnum=1|
 *            |      |
 *      ------+------+
 *
 * Neighbor numbering for q, node=3:
 *
 *            |      |
 *            |nnum=2|nnum=3
 *            |      |
 *            +------+------
 *            |  q   |
 *            |nnum=0|nnum=1
 *            |      |
 *            +------+------
 *
 * \param [in]  q             the quadrant whose possible node \a node neighbor
 *                            will be built.
 * \param [in]  node          the node of the quadrant \a q whose possible node
 *                            neighbor list will be built.  This is given in
 *                            pixel (Morton-) ordering.
 * \param [in]  nnum          the neighbor number with the ordering described
 *                            above, if nnum==0 then it is the corner neighbor.
 * \param [in]  neighor_rlev  the relative level of the neighbor compared to
 *                            the level of \a q.
 * \param [out] neighbor      the neighbor that will be filled.
 * \param [out] neighbor_node the neighbor's node which shares with \a q
 *                            the node \a node.
 */
void                p4est_possible_node_neighbor (const p4est_quadrant_t * q,
                                                  int node, int nnum,
                                                  int neighbor_rlev,
                                                  p4est_quadrant_t * neighbor,
                                                  int *neighbor_node);

#endif /* !P4EST_MESH_H */
