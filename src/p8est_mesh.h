/*
  This file is part of p4est.
  p4est is a C library to manage a parallel collection of quadtrees and/or
  octrees.

  Copyright (C) 2008 Carsten Burstedde, Lucas Wilcox.

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

#ifndef P8EST_MESH_H
#define P8EST_MESH_H

#include <p8est.h>

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
p8est_neighborhood_t;

/** Checks a p8est to see if it is balanced.
 *
 * \param [in] p8est  The p8est to be tested.
 * \return Returns true if balanced, false otherwise.
 */
bool                p8est_is_balanced (p8est_t * p8est);

/** Gets the ghost layer
 *
 * This will gather the quadrants from each neighboring proc to build one layer
 * of face, edge and corner based ghost elements around the ones they own.
 *
 * \param [in] p8est              The forest for which the ghost layer will
 *                                be generated.
 * \param [out] ghost_layer       An array of quadrants which make up the
 *                                ghost layer around \a p4est.  Their piggy1
 *                                data member is filled with their owner's
 *                                tree and processor ids.  Quadrants will be
 *                                ordered in \c p8est_quadrant_compare_piggy
 *                                order.  These will be quadrants inside the
 *                                neighboring tree i.e., \c
 *                                p4est_quadrant_is_inside is true for the
 *                                quadrant and the neighboring tree.
 * \return                        Returns false if it fails due to violated
 *                                2:1 constraints, true otherwise.
 */
bool                p8est_build_ghost_layer (p8est_t * p8est,
                                             sc_array_t * ghost_layer);

/** Populate lists of hanging and anchored nodes.
 */
void                p8est_collect_nodes (p8est_t * p8est,
                                         sc_array_t * ghost_layer);

/** Create neighborhood information.
 */
p8est_neighborhood_t *p8est_neighborhood_new (p8est_t * p8est);

/** Destroy neighborhood information.
 */
void                p8est_neighborhood_destroy (p8est_neighborhood_t * nhood);

#endif /* P8EST_MESH_H */
