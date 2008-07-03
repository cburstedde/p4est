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
 * \warning This function may abort if the forest is not balanced.
 *
 * \param [in] p8est  The p8est to be tested.
 * \return Returns true if balanced, false otherwise.
 */
bool                p8est_is_balanced (p8est_t * p8est);

p8est_neighborhood_t *p8est_neighborhood_new (p8est_t * p8est);
void                p8est_neighborhood_destroy (p8est_neighborhood_t * nhood);

#endif /* P8EST_MESH_H */
