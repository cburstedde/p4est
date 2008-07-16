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

/** Checks a p4est to see if it is balanced.
 *
 * \param [in] p4est  The p4est to be tested.
 * \return Returns true if balanced, false otherwise.
 */
bool                p4est_is_balanced (p4est_t * p4est);

/** Gets the ghost layer
 *
 * This will gather the quadrants from each neighboring proc to build
 * one layer of face and corner based ghost elements around the ones they own.
 *
 * \param [in] p4est              The forest for which the ghost layer will
 *                                be generated.
 * \param [out] ghost_layer       An array of quadrants which make up the
 *                                ghost layer around \a p4est.  Their piggy1
 *                                data member is filled with their owner's
 *                                tree and processor ids.  Quadrants will be
 *                                ordered in \c p4est_quadrant_compare_piggy
 *                                order.  These will be quadrants inside the
 *                                neighboring tree i.e., \c
 *                                p4est_quadrant_is_inside is true for the
 *                                quadrant and the neighboring tree.
 * \return                        Returns false if it fails due to violated
 *                                2:1 constraints, true otherwise.
 */
bool                p4est_build_ghost_layer (p4est_t * p4est,
                                             sc_array_t * ghost_layer);

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

#endif /* !P4EST_MESH_H */
