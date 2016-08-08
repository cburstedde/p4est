/*
  This file is part of p4est.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2010 The University of Texas System
  Additional copyright (C) 2011 individual authors
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

#ifndef P4EST_COMMUNICATION_H
#define P4EST_COMMUNICATION_H

#include <p4est.h>

SC_EXTERN_C_BEGIN;

/** Assign an MPI communicator to p4est; retrieve parallel environment.
 *
 * \param [in] mpicomm    A valid MPI communicator.
 *
 * \note The provided MPI communicator is not owned by p4est.
 */
void                p4est_comm_parallel_env_assign (p4est_t * p4est,
                                                    sc_MPI_Comm mpicomm);

/** Duplicate MPI communicator and replace the current one by the duplicate.
 *
 * \note The duplicated MPI communicator is owned by p4est.
 */
void                p4est_comm_parallel_env_duplicate (p4est_t * p4est);

/** Release MPI communicator if it is owned by p4est.
 */
void                p4est_comm_parallel_env_release (p4est_t * p4est);

/** Replace the current MPI communicator by the one provided as input.
 *
 * \param [in] mpicomm    A valid MPI communicator.
 *
 * \note The provided MPI communicator is not owned by p4est.
 */
void                p4est_comm_parallel_env_replace (p4est_t * p4est,
                                                     sc_MPI_Comm mpicomm);

/** Retrieve parallel environment information.
 */
void                p4est_comm_parallel_env_get_info (p4est_t * p4est);

/** Check if the MPI communicator is valid.
 *
 * \return True if communicator is not NULL communicator, false otherwise.
 */
int                 p4est_comm_parallel_env_is_null (p4est_t * p4est);

/** Reduce MPI communicator to non-empty ranks (i.e., nonzero quadrant counts).
 *
 * \param [in/out] p4est_supercomm  Object which communicator is reduced.
 *                                  Points to NULL if this p4est does not
 *                                  exists.
 *
 * \return True if p4est exists on this MPI rank after reduction.
 */
int                 p4est_comm_parallel_env_reduce (p4est_t **
                                                    p4est_supercomm);

/** Reduce MPI communicator to non-empty ranks and add a group of ranks that
 * will remain in the reduced communicator regardless whether they are empty
 * or not.
 *
 * \param [in/out] p4est_supercomm  Object which communicator is reduced.
 *                                  Points to NULL if this p4est does not
 *                                  exists.
 * \param [in] group_add         Group of ranks that will remain in
 *                               communicator.
 * \param [in] add_to_beginning  If true, ranks will be added to the beginning
 *                               of the reduced communicator, otherwise to the
 *                               end.
 * \param[out] ranks_subcomm     If not null, array of size 'subcommsize' with
 *                               subcommrank->supercommrank map.
 *
 * \return True if p4est exists on this MPI rank after reduction.
 */
int                 p4est_comm_parallel_env_reduce_ext (p4est_t **
                                                        p4est_supercomm,
                                                        sc_MPI_Group
                                                        group_add,
                                                        int add_to_beginning,
                                                        int **ranks_subcomm);

/** Caculate the number and partition of quadrents.
 * \param [in,out] p4est  Adds all \c p4est->local_num_quadrant counters and
 *                        puts cumulative sums in p4est->global_first_quadrant.
 */
void                p4est_comm_count_quadrants (p4est_t * p4est);

/** Distribute the global partition boundaries.
 * \param [in,out] p4est        Fills \c p4est->global_first_position.
 *                              p4est->first_local_tree must be set correctly.
 *                              If this processor is not empty and
 *                              first_quad is NULL, the first quadrant
 *                              of the first local tree must be set correctly.
 * \param [in] first_quad       If not NULL will be used as first quadrant.
 */
void                p4est_comm_global_partition (p4est_t * p4est,
                                                 p4est_quadrant_t *
                                                 first_quad);

/** Compute and distribute the cumulative number of quadrants per tree.
 * \param [in] p4est    This p4est needs to have correct values for
 *                      global_first_quadrant and global_first_position.
 * \param [in,out] pertree      On input, memory for num_trees + 1 numbers.
 *                              On output, the cumulative quadrant counts.
 */
void                p4est_comm_count_pertree (p4est_t * p4est,
                                              p4est_gloidx_t * pertree);

/** Query whether a processor has no quadrants.
 * \param [in] p4est    This forests' global_first_position array must be valid.
 * \param [in] p        Valid processor id.
 * \return              True if and only if processor \p is empty.
 */
int                 p4est_comm_is_empty (p4est_t * p4est, int p);

/** Tests ownershop of a quadrant via p4est->global_first_position.
 * Assumes a tree with no overlaps.
 * \param [in] rank    Rank whose ownership is tested.
 * \return true if rank is the owner.
 */
int                 p4est_comm_is_owner (p4est_t * p4est,
                                         p4est_locidx_t which_tree,
                                         const p4est_quadrant_t * q,
                                         int rank);

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
                                          int full_tree[],
                                          int tree_contact[],
                                          const p4est_quadrant_t ** firstq,
                                          const p4est_quadrant_t ** nextq);

/** Test if the 3x3 neighborhood of a quadrant is owned by this processor.
 * \param [in] p4est            The p4est to work on.
 * \param [in] which_tree       The tree index to work on.
 * \param [in] full_tree[2]     Flags as computed by p4est_comm_tree_info.
 * \param [in] tree_contact[4]  Flags as computed by p4est_comm_tree_info.
 * \param [in] q                The quadrant to be checked.
 * \return          Returns true iff this quadrant's 3x3 neighborhood is owned.
 */
int                 p4est_comm_neighborhood_owned (p4est_t * p4est,
                                                   p4est_locidx_t which_tree,
                                                   int full_tree[],
                                                   int tree_contact[],
                                                   p4est_quadrant_t * q);

/** Evaluates true/false of a flag among processors.
 * \param [in] p4est        The MPI communicator of this p4est will be used.
 * \param [in] flag         The variable to communicate.
 * \param [in] operation    Either sc_MPI_BAND or sc_MPI_BOR (not used bitwise).
 * \return          Returns the logical AND resp. OR of all processors' flags.
 */
int                 p4est_comm_sync_flag (p4est_t * p4est,
                                          int flag, sc_MPI_Op operation);

/** Compute a parallel checksum out of local checksums.
 * \param [in] p4est       The MPI information of this p4est will be used.
 * \param [in] local_crc   Locally computed adler32 checksum.
 * \param [in] local_bytes Number of bytes used for local checksum.
 * \return                 Parallel checksum on rank 0, 0 otherwise.
 */
unsigned            p4est_comm_checksum (p4est_t * p4est,
                                         unsigned local_crc,
                                         size_t local_bytes);

SC_EXTERN_C_END;

#endif /* !P4EST_COMMUNICATION_H */
