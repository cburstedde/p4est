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

/** \file p8est_communication.h
 *
 * Parallel messaging and support code.
 *
 * \ingroup p8est
 */

#ifndef P8EST_COMMUNICATION_H
#define P8EST_COMMUNICATION_H

#include <p8est.h>

SC_EXTERN_C_BEGIN;

/** Given target, find index p such that `gfq[p] <= target < gfq[p + 1]`.
 * \param[in] target    The value that is searched in \a gfq. \a target
 *                      has to satisfy `gfq[0] <= target < gfq[nmemb]`.
 * \param[in] gfq       The sorted array (ascending) in that the function will
 *                      search.
 * \param [in] nmemb    Number of entries in array MINUS ONE.
 * \return              Index p such that `gfq[p] <= target < gfq[p + 1]`.
 * \note                This function differs from \ref p8est_find_partition
 *                      since \ref p8est_find_partition searches for two
 *                      targets using binary search in an optimized way
 *                      but \ref p8est_bsearch_partition only performs a
 *                      single binary search.
 */
int                 p8est_bsearch_partition (p4est_gloidx_t target,
                                             const p4est_gloidx_t * gfq,
                                             int nmemb);

/** Assign an MPI communicator to p8est; retrieve parallel environment.
 *
 * \param [in] mpicomm    A valid MPI communicator.
 *
 * \note The provided MPI communicator is not owned by p8est.
 */
void                p8est_comm_parallel_env_assign (p8est_t * p8est,
                                                    sc_MPI_Comm mpicomm);

/** Duplicate MPI communicator and replace the current one by the duplicate.
 *
 * \note The duplicated MPI communicator is owned by p8est.
 */
void                p8est_comm_parallel_env_duplicate (p8est_t * p8est);

/** Release MPI communicator if it is owned by p8est.
 */
void                p8est_comm_parallel_env_release (p8est_t * p8est);

/** Replace the current MPI communicator by the one provided as input.
 *
 * \param [in] mpicomm    A valid MPI communicator.
 *
 * \note The provided MPI communicator is not owned by p8est.
 */
void                p8est_comm_parallel_env_replace (p8est_t * p8est,
                                                     sc_MPI_Comm mpicomm);

/** Retrieve parallel environment information.
 */
void                p8est_comm_parallel_env_get_info (p8est_t * p8est);

/** Check if the MPI communicator is valid.
 *
 * \return True if communicator is not NULL communicator, false otherwise.
 */
int                 p8est_comm_parallel_env_is_null (p8est_t * p8est);

/** Reduce MPI communicator to non-empty ranks (i.e., nonzero quadrant counts).
 *
 * \param [in,out] p8est_supercomm  Object which communicator is reduced.
 *                                  Points to NULL if this p8est does not
 *                                  exists.
 *
 * \return True if p8est exists on this MPI rank after reduction.
 */
int                 p8est_comm_parallel_env_reduce (p8est_t **
                                                    p8est_supercomm);

/** Reduce MPI communicator to non-empty ranks and add a group of ranks that
 * will remain in the reduced communicator regardless whether they are empty
 * or not.
 *
 * \param [in,out] p8est_supercomm  Object which communicator is reduced.
 *                                  Points to NULL if this p8est does not
 *                                  exists.
 * \param [in] group_add         Group of ranks that will remain in
 *                               communicator.
 * \param [in] add_to_beginning  If true, ranks will be added to the beginning
 *                               of the reduced communicator, otherwise to the
 *                               end.
 * \param[out] ranks_subcomm     If not null, array of size 'subcommsize' with
 *                               subcommrank->supercommrank map.
 *
 * \return True if p8est exists on this MPI rank after reduction.
 */
int                 p8est_comm_parallel_env_reduce_ext (p8est_t **
                                                        p8est_supercomm,
                                                        sc_MPI_Group
                                                        group_add,
                                                        int add_to_beginning,
                                                        int **ranks_subcomm);

/** Calculate the number and partition of quadrants.
 * \param [in,out] p8est  Adds all \c p8est->local_num_quadrant counters and
 *                        puts cumulative sums in p8est->global_first_quadrant.
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

/** Calculate the global fist quadrant array for a uniform partition.
 *
 * \param [in] global_num_quadrants   The global number of quadrants.
 * \param [in] mpisize                The number of MPI ranks.
 * \param [in,out] gfq                At least allocated mpisize + 1
 *                                    p4est_gloidx_t. This array is
 *                                    filled with the global first
 *                                    quadrant array assuming a
 *                                    uniform partition.
 */
void                p8est_comm_global_first_quadrant (p4est_gloidx_t
                                                      global_num_quadrants,
                                                      int mpisize,
                                                      p4est_gloidx_t * gfq);

/** Compute and distribute the cumulative number of quadrants per tree.
 * \param [in] p8est    This p8est needs to have correct values for
 *                      global_first_quadrant and global_first_position.
 * \param [in,out] pertree      On input, memory for num_trees + 1 numbers.
 *                              On output, the cumulative quadrant counts.
 */
void                p8est_comm_count_pertree (p8est_t * p8est,
                                              p4est_gloidx_t * pertree);

/** Query whether a processor has no quadrants.
 * \param [in] p8est    This forests' global_first_position array must be valid.
 * \param [in] p        Valid processor id.
 * \return              True if and only if processor \p is empty.
 */
int                 p8est_comm_is_empty (p8est_t *p8est, int p);

/** Query whether a processor has no quadrants.
 * \param [in] gfq          An array encoding the partition offsets in the
 *                          global quadrant array; length \a num_procs + 1.
 * \param [in] num_procs    Number of processes in the partition.
 * \param [in] p            Valid 0 <= \a p < \a num_procs.
 * \return              True if and only if processor \a p is empty.
 */
int                 p8est_comm_is_empty_gfq (const p4est_gloidx_t *gfq,
                                             int num_procs, int p);

/** Query whether a processor has no quadrants.
 * \param [in] gfp          An array encoding the partition shape.
 *                          Non-decreasing; length \a num_procs + 1.
 * \param [in] num_procs    Number of processes in the partition.
 * \param [in] p            Valid 0 <= \a p < \a num_procs.
 * \return              True if and only if processor \a p is empty.
 */
int                 p8est_comm_is_empty_gfp (const p8est_quadrant_t *gfp,
                                             int num_procs, int p);

/** Test whether a quadrant is fully contained in a rank's owned region.
 * This function may return false when \ref p8est_comm_is_owner returns true.
 * \param [in] rank    Rank whose ownership is tested.
 *                     Assumes a forest with no overlaps.
 * \return true if rank is the owner of the whole area of the quadrant.
 */
int                 p8est_comm_is_contained (p8est_t * p8est,
                                             p4est_locidx_t which_tree,
                                             const p8est_quadrant_t * q,
                                             int rank);

/** Test ownership of a quadrant via p8est->global_first_position.
 * The quadrant is considered owned if its first descendant is owned.
 * Thus, a positive result occurs even if its last descendant overlaps
 * a higher process.
 * \param [in] p8est        Valid forest.
 * \param [in] which_tree   Valid tree number wrt. the forest.
 * \param [in] q            Valid quadrant wrt. the forest.
 * \param [in] rank         Rank whose ownership is tested.
 * \return      True if rank is the owner of the first descendant.
 */
int                 p8est_comm_is_owner (p8est_t *p8est,
                                         p4est_locidx_t which_tree,
                                         const p8est_quadrant_t *q,
                                         int rank);

/** Test ownership of a quadrant via a global_first_position array.
 * This array encodes part of the partition of a valid forest object.
 * The quadrant is considered owned if its first descendant is owned.
 * Thus, a positive result occurs even if its last descendant overlaps
 * a higher process.
 * \param [in] gfp          Position array of length \a num_procs + 1.
 * \param [in] num_procs    Number of processes in this context.
 * \param [in] num_trees    Number of trees in this context.
 * \param [in] which_tree   Valid tree number wrt. the forest.
 * \param [in] q            Valid quadrant wrt. the forest.
 * \param [in] rank         Rank whose ownership is tested.
 * \return      True if rank is the owner of the first descendant.
 */
int                 p8est_comm_is_owner_gfp
  (const p8est_quadrant_t *gfp, int num_procs, p4est_topidx_t num_trees,
   p4est_locidx_t which_tree, const p8est_quadrant_t *q, int rank);

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
 * \param [out] full_tree       Full ownership of beginning and end of tree.
 * \param [out] tree_contact    True if there are neighbors across the face.
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
 * \param [in] full_tree        Flags as computed by p8est_comm_tree_info.
 * \param [in] tree_contact     Flags as computed by p8est_comm_tree_info.
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

/** Compute a parallel partition-independent checksum out of local checksums.
 * This checksum depends on the global refinement topology.
 * It does not depend on how the mesh is partitioned.
 * The result is available on all processors.
 * \param [in] p8est       The MPI information of this p8est will be used.
 * \param [in] local_crc   Locally computed adler32 checksum.
 * \param [in] local_bytes Number of bytes used for local checksum.
 * \return                 Parallel checksum on all processors.
 */
unsigned            p8est_comm_checksum (p8est_t * p8est,
                                         unsigned local_crc,
                                         size_t local_bytes);

/** Context data to allow for split begin/end data transfer. */
typedef struct p8est_transfer_context
{
  int                 variable;
  int                 num_senders;
  int                 num_receivers;
  sc_MPI_Request     *recv_req;
  sc_MPI_Request     *send_req;
}
p8est_transfer_context_t;

/** Transfer data associated with one forest partition to another.
 * In \ref p8est_partition, each quadrant's user data is transferred.
 * If the application maintains per-quadrant data outside of the p8est object,
 * this function can be used to transfer it, matching the call to partition.
 * This variant of the function assumes that the quadrant data size is fixed.
 * It sends point-to-point messages only and is blocking collective.
 * There is a split collective version; see the functions
 * \ref p8est_transfer_fixed_begin and \ref p8est_transfer_fixed_end.
 * \param [in] dest_gfq     The target partition encoded as a \b
 *                          p8est->global_first_quadrant array.  Has \b mpisize
 *                          + 1 members, must be non-decreasing and satisfy
 *                          gfq[0] == 0, gfq[mpisize] == global_num_quadrants.
 * \param [in] src_gfq      The original partition, analogous to \b dest_gfq.
 * \param [in] mpicomm      The communicator to use.
 *                          Its mpisize must match \b dest_gfq and \b src_gfq.
 * \param [in] tag          This tag is used in all messages.  The user must
 *                          guarantee that \b mpicomm and \b tag do not
 *                          conflict with other messages in transit.
 * \param [out] dest_data   User-allocated memory of size \b data_size * \b
 *                          dest->local_num_quadrants is received into.
 * \param [in] src_data     User-allocated memory of size \b data_size * \b
 *                          src->local_num_quadrants bytes is sent from.
 * \param [in] data_size    Fixed data size per quadrant.
 */
void                p8est_transfer_fixed (const p4est_gloidx_t * dest_gfq,
                                          const p4est_gloidx_t * src_gfq,
                                          sc_MPI_Comm mpicomm, int tag,
                                          void *dest_data,
                                          const void *src_data,
                                          size_t data_size);

/** Initiate a fixed-size data transfer between partitions.
 * See \ref p8est_transfer_fixed for a full description.
 * Must be matched with \ref p8est_transfer_fixed_end for completion.
 * All parameters must stay alive until the completion has been called.
 * \param [in] dest_gfq     The target partition encoded as a \b
 *                          p8est->global_first_quadrant array.  Has \b mpisize
 *                          + 1 members, must be non-decreasing and satisfy
 *                          gfq[0] == 0, gfq[mpisize] == global_num_quadrants.
 * \param [in] src_gfq      The original partition, analogous to \b dest_gfq.
 * \param [in] mpicomm      The communicator to use.
 *                          Its mpisize must match \b dest_gfq and \b src_gfq.
 * \param [in] tag          This tag is used in all messages.  The user must
 *                          guarantee that \b mpicomm and \b tag do not
 *                          conflict with other messages in transit.
 * \param [out] dest_data   User-allocated memory of size \b data_size * \b
 *                          dest->local_num_quadrants bytes is received into.
 *                          It must not be accessed before completion with
 *                          \ref p8est_transfer_fixed_end.
 * \param [in] src_data     User-allocated memory of size \b data_size * \b
 *                          src->local_num_quadrants bytes is sent from.
 *                          It must not be accessed before completion with
 *                          \ref p8est_transfer_fixed_end.
 * \param [in] data_size    Fixed data size per quadrant.
 * \return                  The context object must be passed to the matching
 *                          call to \ref p8est_transfer_fixed_end.
 */
p8est_transfer_context_t *p8est_transfer_fixed_begin (const p4est_gloidx_t *
                                                      dest_gfq,
                                                      const p4est_gloidx_t *
                                                      src_gfq,
                                                      sc_MPI_Comm mpicomm,
                                                      int tag,
                                                      void *dest_data,
                                                      const void *src_data,
                                                      size_t data_size);

/** Complete a fixed-size data transfer between partitions.
 * Waits for remaining messages to complete and frees the transfer context.
 * \param [in] tc       Context data from \ref p8est_transfer_fixed_begin.
 *                      Is deallocated before this function returns.
 */
void                p8est_transfer_fixed_end (p8est_transfer_context_t * tc);

/** Transfer variable-size quadrant data between partitions.
 * (See \ref p8est_transfer_fixed that is optimized for fixed-size data.)
 * The destination process may not know the data size for the elements it
 * receives.  In this case the sizes need to be obtained separately in advance,
 * for example by calling \ref p8est_transfer_fixed with \b src_sizes as
 * payload data, or alternatively its split begin/end versions.
 * \param [in] dest_gfq     The target partition encoded as a \b
 *                          p8est->global_first_quadrant array.  Has \b mpisize
 *                          + 1 members, must be non-decreasing and satisfy
 *                          gfq[0] == 0, gfq[mpisize] == global_num_quadrants.
 * \param [in] src_gfq      The original partition, analogous to \b dest_gfq.
 * \param [in] mpicomm      The communicator to use.
 *                          Its mpisize must match \b dest_gfq and \b src_gfq.
 * \param [in] tag          This tag is used in all messages.  The user must
 *                          guarantee that \b mpicomm and \b tag do not
 *                          conflict with other messages in transit.
 * \param [out] dest_data   User-allocated memory of
 *                          sum_{i in \b dest->local_num_quadrants} \b
 *                          dest_sizes [i] many bytes is received into.
 *                          See below about how to choose its size.
 * \param [in] dest_sizes   User-allocated memory of one integer for each
 *                          quadrant, storing the data size to receive for it.
 *                          We understand that the sizes are often not known a
 *                          priori, in which case they can be obtained by a
 *                          prior call to \ref p8est_transfer_fixed.
 *                          Optionally the split begin/end versions can be used
 *                          for added flexibility and overlapping of messages.
 *                          We use the type int to minimize the message size,
 *                          and to conform to MPI that has no type for size_t.
 * \param [in] src_data     User-allocated memory of
 *                          sum_{i in \b src->local_num_quadrants} \b
 *                          src_sizes [i] many bytes is sent from.
 * \param [in] src_sizes    User-allocated memory of one integer for each
 *                          quadrant, storing the data size to send for it.
 *                          We use the type int to minimize the message size,
 *                          and to conform to MPI that has no type for size_t.
 */
void                p8est_transfer_custom (const p4est_gloidx_t * dest_gfq,
                                           const p4est_gloidx_t * src_gfq,
                                           sc_MPI_Comm mpicomm, int tag,
                                           void *dest_data,
                                           const int *dest_sizes,
                                           const void *src_data,
                                           const int *src_sizes);

/** Initiate a variable-size data transfer between partitions.
 * See \ref p8est_transfer_custom for a full description.
 * Must be matched with \ref p8est_transfer_custom_end for completion.
 * All parameters must stay alive until the completion has been called.
 * \param [in] dest_gfq     The target partition encoded as a \b
 *                          p8est->global_first_quadrant array.  Has \b mpisize
 *                          + 1 members, must be non-decreasing and satisfy
 *                          gfq[0] == 0, gfq[mpisize] == global_num_quadrants.
 * \param [in] src_gfq      The original partition, analogous to \b dest_gfq.
 * \param [in] mpicomm      The communicator to use.
 *                          Its mpisize must match \b dest_gfq and \b src_gfq.
 * \param [in] tag          This tag is used in all messages.  The user must
 *                          guarantee that \b mpicomm and \b tag do not
 *                          conflict with other messages in transit.
 * \param [out] dest_data   User-allocated memory of
 *                          sum_{i in \b dest->local_num_quadrants} \b
 *                          dest_sizes [i] many bytes is received into.
 *                          It must not be accessed before completion with
 *                          \ref p8est_transfer_custom_end.
 *                          See below about how to choose its size.
 * \param [in] dest_sizes   User-allocated memory of one integer for each
 *                          quadrant, storing the data size to receive for it.
 *                          We understand that the sizes are often not known a
 *                          priori, in which case they can be obtained by a
 *                          prior call to \ref p8est_transfer_fixed.
 *                          Optionally the split begin/end versions can be used
 *                          for added flexibility and overlapping of messages.
 *                          We use the type int to minimize the message size,
 *                          and to conform to MPI that has no type for size_t.
 * \param [in] src_data     User-allocated memory of
 *                          sum_{i in \b src->local_num_quadrants} \b
 *                          src_sizes [i] many bytes is sent from.
 *                          It must not be accessed before completion with
 *                          \ref p8est_transfer_custom_end.
 * \param [in] src_sizes    User-allocated memory of one integer for each
 *                          quadrant, storing the data size to send for it.
 *                          We use the type int to minimize the message size,
 *                          and to conform to MPI that has no type for size_t.
 * \return                  The context object must be passed to the matching
 *                          call to \ref p8est_transfer_custom_end.
 */
p8est_transfer_context_t *p8est_transfer_custom_begin (const p4est_gloidx_t *
                                                       dest_gfq,
                                                       const p4est_gloidx_t *
                                                       src_gfq,
                                                       sc_MPI_Comm mpicomm,
                                                       int tag,
                                                       void *dest_data,
                                                       const int *dest_sizes,
                                                       const void *src_data,
                                                       const int *src_sizes);

/** Complete a variable-size data transfer between partitions.
 * Waits for remaining messages to complete and frees the transfer context.
 * \param [in] tc       Context data from \ref p8est_transfer_custom_begin.
 *                      Is deallocated before this function returns.
 */
void                p8est_transfer_custom_end (p8est_transfer_context_t * tc);

/** Transfer variable-count item data between partitions.
 * Each quadrant may have a different number of items (including 0).
 * (See \ref p8est_transfer_fixed that is optimized for fixed-count data,
 *  and \ref p8est_transfer_custom for data that is not itemized at all.)
 * The destination process may not know the item count for the elements it
 * receives.  In this case the counts need to be obtained separately in advance,
 * for example by calling \ref p8est_transfer_fixed with \b src_counts as
 * payload data, or alternatively its split begin/end versions.
 * \param [in] dest_gfq     The target partition encoded as a \b
 *                          p8est->global_first_quadrant array.  Has \b mpisize
 *                          + 1 members, must be non-decreasing and satisfy
 *                          gfq[0] == 0, gfq[mpisize] == global_num_quadrants.
 * \param [in] src_gfq      The original partition, analogous to \b dest_gfq.
 * \param [in] mpicomm      The communicator to use.
 *                          Its mpisize must match \b dest_gfq and \b src_gfq.
 * \param [in] tag          This tag is used in all messages.  The user must
 *                          guarantee that \b mpicomm and \b tag do not
 *                          conflict with other messages in transit.
 * \param [out] dest_data   User-allocated memory of
 *                          sum_{i in \b dest->local_num_quadrants} item_size
 *                          * dest_counts[i] many bytes is received into.
 *                          See below about how to choose its size.
 *                          If dest has data to transfer, must be non-NULL.
 * \param [in] dest_counts  User-allocated memory of one integer for each
 *                          quadrant, storing the item count to receive for it.
 *                          We understand that the counts are often not known a
 *                          priori, in which case they can be obtained by a
 *                          prior call to \ref p8est_transfer_fixed.
 *                          Optionally the split begin/end versions can be used
 *                          for added flexibility and overlapping of messages.
 *                          We use the type int to minimize the message size,
 *                          and to conform to MPI that has no type for size_t.
 *                          If dest has quadrants, must be non-NULL.
 * \param [in] src_data     User-allocated memory of
 *                          sum_{i in \b src->local_num_quadrants} item_size
 *                          * src_counts[i] many bytes is sent from.
 *                          If src has data to transfer, must be non-NULL.
 * \param [in] src_counts   User-allocated memory of one integer for each
 *                          quadrant, storing the item count to send for it.
 *                          We use the type int to minimize the message size,
 *                          and to conform to MPI that has no type for size_t.
 *                          If src has quadrants, must be non-NULL.
 * \param [in] item_size    Data size for each item in bytes.
 */
void                p8est_transfer_items
  (const p4est_gloidx_t * dest_gfq, const p4est_gloidx_t * src_gfq,
   sc_MPI_Comm mpicomm, int tag,
   void *dest_data, const int *dest_counts,
   const void *src_data, const int *src_counts, size_t item_size);

/** Initiate a variable-count item transfer between partitions.
 * See \ref p8est_transfer_items for a full description.
 * This functions calls asynchronous MPI send/receive and returns.
 * Must be matched with \ref p8est_transfer_items_end for completion,
 * which calls blocking MPI wait until all messages have been processed.
 * All parameters must stay alive until the completion has been called.
 */
p8est_transfer_context_t *p8est_transfer_items_begin
  (const p4est_gloidx_t * dest_gfq, const p4est_gloidx_t * src_gfq,
   sc_MPI_Comm mpicomm, int tag,
   void *dest_data, const int *dest_counts,
   const void *src_data, const int *src_counts, size_t item_size);

/** Complete a variable-count item transfer between partitions.
 * Waits for remaining messages to complete and frees the transfer context.
 * \param [in] tc       Context data from \ref p8est_transfer_items_begin.
 *                      Is deallocated before this function returns.
 */
void                p8est_transfer_items_end (p8est_transfer_context_t * tc);

/** Complete any of the transfer_begin functions.
 * The specialized transfer_end functions are recommended over this one
 * for slightly stricter error checking: \ref p8est_transfer_fixed_end,
 * \ref p8est_transfer_custom_end, and \ref p8est_transfer_items_end.
 *
 * \param [in,out] tc       A valid context from one of the begin functions.
 *                          This function waits for remaining communications
 *                          to complete and frees the transfer context.
 */
void                p8est_transfer_end (p8est_transfer_context_t * tc);

SC_EXTERN_C_END;

#endif /* !P8EST_COMMUNICATION_H */
