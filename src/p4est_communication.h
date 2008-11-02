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

#ifndef P4EST_COMMUNICATION_H
#define P4EST_COMMUNICATION_H

#include <p4est.h>

SC_EXTERN_C_BEGIN;

typedef enum
{
  P4EST_COMM_BALANCE_FIRST_COUNT = 1,
  P4EST_COMM_BALANCE_FIRST_LOAD,
  P4EST_COMM_BALANCE_SECOND_COUNT,
  P4EST_COMM_BALANCE_SECOND_LOAD,
  P4EST_COMM_PARTITION_GIVEN,
  P4EST_COMM_PARTITION_WEIGHTED_LOW,
  P4EST_COMM_PARTITION_WEIGHTED_HIGH,
  P4EST_COMM_GHOST_COUNT,
  P4EST_COMM_GHOST_LOAD,
  P4EST_COMM_NODES_QUERY,
  P4EST_COMM_NODES_REPLY,
  P4EST_COMM_SAVE
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

/** Computes information about a tree being fully owned.
 * This is determined separately for the beginning and end of the tree.
 * \param [in] p4est            The p4est to work on.
 * \param [in] which_tree       The tree in question must be partially owned.
 * \param [out] full_tree[2]    Full ownership of beginning and end of tree.
 * \param [out] tree_contact[4] True if there are neighbors across the face.
 * \param [out] firstq          Smallest possible first quadrant on this core.
 * \param [out] nextq           Smallest possible first quadrant on next core.
 *                          Any of tree_contact, firstq and nextq may be NULL.
 */
void                p4est_comm_tree_info (p4est_t * p4est,
                                          p4est_locidx_t which_tree,
                                          bool full_tree[],
                                          bool tree_contact[],
                                          const p4est_quadrant_t ** firstq,
                                          const p4est_quadrant_t ** nextq);

/** Test if the 3x3 neighborhood of a quadrant is owned by this processor.
 * \param [in] p4est             The p4est to work on.
 * \param [in] which_tree        The tree index to work on.
 * \param [in] full_tree[2]      Flags as computed by p4est_comm_tree_info.
 * \param [in] tree_contact[4]   Flags as computed by p4est_comm_tree_info.
 * \param [in] q                 The quadrant to be checked.
 * \return   Returns true iff this quadrant's 3x3 neighborhood is owned.
 */
bool                p4est_comm_neighborhood_owned (p4est_t * p4est,
                                                   p4est_locidx_t which_tree,
                                                   bool full_tree[],
                                                   bool tree_contact[],
                                                   p4est_quadrant_t * q);

/** Evaluates a boolean flag among processors.
 * \param [in] p4est       The MPI communicator of this p4est will be used.
 * \param [in] flag        The boolean flag to communicate.
 * \param [in] operation   Either MPI_BAND or MPI_BOR.
 * \return     Returns the AND resp. OR of all processors' boolean flags.
 */
bool                p4est_comm_sync_flag (p4est_t * p4est,
                                          bool flag, MPI_Op operation);

SC_EXTERN_C_END;

#endif /* !P4EST_COMMUNICATION_H */
