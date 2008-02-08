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

#ifndef __P4EST_COMMUNICATION_H__
#define __P4EST_COMMUNICATION_H__

#include <p4est.h>

typedef enum
{
  P4EST_COMM_BALANCE_FIRST_COUNT = 1,
  P4EST_COMM_BALANCE_FIRST_LOAD,
  P4EST_COMM_BALANCE_SECOND_COUNT,
  P4EST_COMM_BALANCE_SECOND_LOAD,
  P4EST_COMM_PARTITION_GIVEN,
  P4EST_COMM_PARTITION_WEIGHTED_LOW,
  P4EST_COMM_PARTITION_WEIGHTED_HIGH,
}
p4est_comm_tag_t;

/** Caculate the number and partition of quadrents.
 * \param [in,out] p4est  Adds all \c p4est->local_num_quadrant counters.
 *                        This also fills in \c p4est->global_last_quad_index
 *                        with the correct values.
 */
void                p4est_comm_count_quadrants (p4est_t * p4est);

/** Distribute the global partition boundaries.
 * \param [in,out] p4est  Fills \c p4est->global_first_position.
 */
void                p4est_comm_global_partition (p4est_t * p4est);

/** Searches the owner of a quadrant via p4est->global_first_position.
 * Assumes a tree with no overlaps.
 * \param [in] guess   Initial guess for the search.
 * \return Returns the processor id of the owner.
 */
int                 p4est_comm_find_owner (p4est_t * p4est,
                                           p4est_locidx_t which_tree,
                                           const p4est_quadrant_t * q,
                                           int guess);

#endif /* !__P4EST_COMMUNICATION_H__ */
