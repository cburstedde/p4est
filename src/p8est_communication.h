/*
  This file is part of p4est.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2010 The University of Texas System
  Written by Carsten Burstedde, Lucas C. Wilcox, and Tobin Isaac

  p4est is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  p4est is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with p4est; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

#ifndef P8EST_COMMUNICATION_H
#define P8EST_COMMUNICATION_H

#include <p8est.h>

SC_EXTERN_C_BEGIN;

typedef enum
{
  P8EST_COMM_COUNT_PERTREE = 1,
  P8EST_COMM_BALANCE_FIRST_COUNT,
  P8EST_COMM_BALANCE_FIRST_LOAD,
  P8EST_COMM_BALANCE_SECOND_COUNT,
  P8EST_COMM_BALANCE_SECOND_LOAD,
  P8EST_COMM_PARTITION_GIVEN,
  P8EST_COMM_PARTITION_WEIGHTED_LOW,
  P8EST_COMM_PARTITION_WEIGHTED_HIGH,
  P8EST_COMM_PARTITION_CORRECTION,
  P8EST_COMM_GHOST_COUNT,
  P8EST_COMM_GHOST_LOAD,
  P8EST_COMM_GHOST_EXCHANGE,
  P8EST_COMM_GHOST_EXPAND_COUNT,
  P8EST_COMM_GHOST_EXPAND_LOAD,
  P8EST_COMM_GHOST_SUPPORT_COUNT,
  P8EST_COMM_GHOST_SUPPORT_LOAD,
  P8EST_COMM_GHOST_CHECKSUM,
  P8EST_COMM_NODES_QUERY,
  P8EST_COMM_NODES_REPLY,
  P8EST_COMM_SAVE,
  P8EST_COMM_LNODES_TEST,
  P8EST_COMM_LNODES_PASS,
  P8EST_COMM_LNODES_OWNED,
  P8EST_COMM_LNODES_ALL
}
p8est_comm_tag_t;

/** Caculate the number and partition of quadrents.
 * \param [in,out] p8est  Adds all \c p8est->local_num_quadrant counters and
 *                        puts cumulative sums in p4est->global_first_quadrant.
 */
void                p8est_comm_count_quadrants (p8est_t * p8est);

/** Distribute the global partition boundaries.
 * \param [in,out] p8est        Fills \c p8est->global_first_position.
 *                              p8est->first_local_tree must be set correctly.
 *                              If this processor is not empty and
 *                              first_quad is NULL, the first quadrant
 *                              of the first local tree must be set correctly.
 * \param [in] first_quad       If not NULL will be used as first quadrant.
 */
void                p8est_comm_global_partition (p8est_t * p8est,
                                                 p8est_quadrant_t *
                                                 first_quad);

/** Compute and distribute the cumulative number of quadrants per tree.
 * \param [in] p8est    This p8est needs to have correct values for
 *                      global_first_quadrant and global_first_position.
 * \param [in,out] pertree      On input, memory for num_trees + 1 numbers.
 *                              On output, the cumulative quadrant counts.
 */
void                p8est_comm_count_pertree (p8est_t * p8est,
                                              p4est_gloidx_t * pertree);

/** Tests ownershop of a quadrant via p8est->global_first_position.
 * Assumes a tree with no overlaps.
 * \param [in] rank    Rank whose ownership is tested.
 * \return true if rank is the owner.
 */
int                 p8est_comm_is_owner (p8est_t * p8est,
                                         p4est_locidx_t which_tree,
                                         const p8est_quadrant_t * q,
                                         int rank);

/** Searches the owner of a quadrant via p8est->global_first_position.
 * Assumes a tree with no overlaps.
 * \param [in] guess   Initial guess for the search.
 * \return Returns the processor id of the owner.
 */
int                 p8est_comm_find_owner (p8est_t * p8est,
                                           p4est_locidx_t which_tree,
                                           const p8est_quadrant_t * q,
                                           int guess);

/** Computes information about a tree being fully owned.
 * This is determined separately for the beginning and end of the tree.
 * \param [in] p8est            The p8est to work on.
 * \param [in] which_tree       The tree in question must be partially owned.
 * \param [out] full_tree[2]    Full ownership of beginning and end of tree.
 * \param [out] tree_contact[6] True if there are neighbors across the face.
 * \param [out] firstq          Smallest possible first quadrant on this core.
 * \param [out] nextq           Smallest possible first quadrant on next core.
 *                          Any of tree_contact, firstq and nextq may be NULL.
 */
void                p8est_comm_tree_info (p8est_t * p8est,
                                          p4est_locidx_t which_tree,
                                          int full_tree[],
                                          int tree_contact[],
                                          const p8est_quadrant_t ** firstq,
                                          const p8est_quadrant_t ** nextq);

/** Test if the 3x3 neighborhood of a quadrant is owned by this processor.
 * \param [in] p8est            The p8est to work on.
 * \param [in] which_tree       The tree index to work on.
 * \param [in] full_tree[2]     Flags as computed by p4est_comm_tree_info.
 * \param [in] tree_contact[6]  Flags as computed by p4est_comm_tree_info.
 * \param [in] q                The quadrant to be checked.
 * \return          Returns true iff this quadrant's 3x3 neighborhood is owned.
 */
int                 p8est_comm_neighborhood_owned (p8est_t * p8est,
                                                   p4est_locidx_t which_tree,
                                                   int full_tree[],
                                                   int tree_contact[],
                                                   p8est_quadrant_t * q);

/** Evaluates true/false of a flag among processors.
 * \param [in] p8est        The MPI communicator of this p8est will be used.
 * \param [in] flag         The variable to communicate.
 * \param [in] operation    Either sc_MPI_BAND or sc_MPI_BOR (not used bitwise).
 * \return          Returns the logical AND resp. OR of all processors' flags.
 */
int                 p8est_comm_sync_flag (p8est_t * p8est,
                                          int flag, sc_MPI_Op operation);

/** Compute a parallel checksum out of local checksums.
 * \param [in] p8est       The MPI information of this p8est will be used.
 * \param [in] local_crc   Locally computed adler32 checksum.
 * \param [in] local_bytes Number of bytes used for local checksum.
 * \return                 Parallel checksum on rank 0, 0 otherwise.
 */
unsigned            p8est_comm_checksum (p8est_t * p8est,
                                         unsigned local_crc,
                                         size_t local_bytes);

SC_EXTERN_C_END;

#endif /* !P8EST_COMMUNICATION_H */
