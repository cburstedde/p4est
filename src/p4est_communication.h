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

#ifndef P4EST_COMMUNICATION_H
#define P4EST_COMMUNICATION_H

#include <p4est.h>

SC_EXTERN_C_BEGIN;

void                p4est_comm_parallel_env_create (p4est_t * p4est,
                                                    sc_MPI_Comm mpicomm);

void                p4est_comm_parallel_env_free (p4est_t * p4est);

int                 p4est_comm_parallel_env_is_null (p4est_t * p4est);

void                p4est_comm_parallel_env_assign (p4est_t * p4est,
                                                    sc_MPI_Comm mpicomm);

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

/** Defines how the communicator is obtained when transfering data.
 * It is used in \ref p4est_transfer_fixed and \ref p4est_transfer_custom.
 */
typedef enum p4est_transfer_comm
{
  P4EST_TRANSFER_COMM_SRC,      /**< Use communicator from source forest. */
  P4EST_TRANSFER_COMM_DEST,     /**< Use communicator from target forest. */
  P4EST_TRANSFER_COMM_SRC_DUP,  /**< Duplicate source communicator. */
  P4EST_TRANSFER_COMM_DEST_DUP, /**< Duplicate target communicator. */
  P4EST_TRANSFER_COMM_EXTERNAL  /**< Use user-specified communicator. */
}
p4est_transfer_comm_t;

/** Transfer data associated with one forest to a partitioned one.
 * In \ref p4est_partition, each quadrant's user data is transferred.
 * If the application maintains per-quadrant data outside of p4est,
 * this function can be used to transfer it, matching the call to partition.
 * This variant of the function assumes that the quadrant data size is fixed.
 * It sends point-to-point messages only and is blocking collective.
 * There is a split collective version; see the functions
 * \ref p4est_transfer_fixed_begin and \ref p4est_transfer_fixed_end.
 * \param [in] dest         This forest defines the target partition.
 *                          \b dest must have been derived from \b src by a
 *                          call to \ref p4est_partition.
 *                          It is ok to use \ref p4est_copy, too.
 * \param [in] src          This forest defines the original partition.
 * \param [in] which_comm   This enumeration defines how the communicator is
 *                          obtained that should be used inside this function.
 *                          If it is to be derived from either \b dest or \b
 *                          src, there is no need to specify it in \b mpicomm.
 *                          When it is duped, it will be freed internally.
 * \param [in] mpicomm      If the communicator to be used is user provided,
 *                          by specifying P4EST_TRANSFER_COMM_EXTERNAL for
 *                          \b which_comm, it must be passed in here.  Then it
 *                          must have the same size and rank as the ones stored
 *                          in \b dest and \b src.
 * \param [in] tag          This tag is used in all messages.
 * \param [in,out] dest_data    User-allocated memory of size \b data_size * \b
 *                              dest->local_num_quadrants bytes.
 * \param [in] src_data         User-allocated memory of size \b data_size * \b
 *                              src->local_num_quadrants bytes.
 * \param [in] data_size        Fixed data size per quadrant.
 */
void                p4est_transfer_fixed (p4est_t * dest, p4est_t * src,
                                          p4est_transfer_comm_t which_comm,
                                          sc_MPI_Comm mpicomm, int tag,
                                          void *dest_data,
                                          const void *src_data,
                                          size_t data_size);

typedef struct p4est_transfer_context
{
  p4est_t            *dest;
  p4est_t            *src;
  p4est_transfer_comm_t which_comm;
  sc_MPI_Comm         mpicomm;
  int                 tag;
  void               *dest_data;
  const void         *src_data;
  size_t              data_size;
}
p4est_transfer_context_t;

p4est_transfer_context_t *p4est_transfer_fixed_begin (p4est_t * dest,
                                                      p4est_t * src,
                                                      p4est_transfer_comm_t
                                                      which_comm,
                                                      sc_MPI_Comm mpicomm,
                                                      int tag,
                                                      void *dest_data,
                                                      const void *src_data,
                                                      size_t data_size);

void                p4est_transfer_fixed_end (p4est_transfer_context_t * tc);

SC_EXTERN_C_END;

#endif /* !P4EST_COMMUNICATION_H */
