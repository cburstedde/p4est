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

/** This mesh structure holds complete neighborhood information.
 * cumulative_count[i]   is the sum of local quadrants in the
 *                       local trees 0..i-1. i == local_num_trees is allowed.
 * element_offsets[i]    is the offset into local_neighbors for local
 *                       element i. i == local_num_quadrants is allowed
 *                       and contains the number of stored local neighbors.
 */
typedef struct
{
  p4est_locidx_t     *cumulative_count;
  p4est_locidx_t     *element_offsets;
  sc_array_t         *local_neighbors;
}
p4est_neighborhood_t;

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
                                                p4est_locidx_t *
                                                num_uniq_local_vertices,
                                                p4est_locidx_t *
                                                quadrant_to_local_vertex);

/** Populate lists of hanging and anchored nodes.
 */
void                p4est_collect_nodes (p4est_t * p4est,
                                         sc_array_t * ghost_layer);

/** Create neighborhood information.
 */
p4est_neighborhood_t *p4est_neighborhood_new (p4est_t * p4est);

/** Destroy neighborhood information.
 */
void                p4est_neighborhood_destroy (p4est_neighborhood_t * nhood);

#endif /* !P4EST_MESH_H */
