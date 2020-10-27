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

/** \file p4est_algorithms.h
 *
 * Routines for managing quadrants as elements of trees and subtrees.
 * In addition, some high level algorithms such as \ref p4est_partition_given.
 *
 * \ingroup p4est
 */

#ifndef P4EST_ALGORITHMS_H
#define P4EST_ALGORITHMS_H

#include <p4est.h>
#include <p4est_extended.h>

SC_EXTERN_C_BEGIN;

/** Alloc and initialize the user data of a valid quadrant.
 * \param [in]  which_tree 0-based index of this quadrant's tree.
 * \param [in,out]  quad       The quadrant to be initialized.
 * \param [in]  init_fn    User-supplied callback function to init data.
 */
void                p4est_quadrant_init_data (p4est_t * p4est,
                                              p4est_topidx_t which_tree,
                                              p4est_quadrant_t * quad,
                                              p4est_init_t init_fn);

/** Free the user data of a valid quadrant.
 * \param [in,out]  quad The quadrant whose data shall be freed.
 */
void                p4est_quadrant_free_data (p4est_t * p4est,
                                              p4est_quadrant_t * quad);

/** Computes a machine-independent checksum of a list of quadrants.
 * \param [in] quadrants       Array of quadrants.
 * \param [in,out] checkarray  Temporary array of elem_size 4.
 *                             Will be resized to quadrants->elem_count * 3.
 *                             If it is NULL, will be allocated internally.
 * \param [in] first_quadrant  Index of the quadrant to start with.
 *                             Can be between 0 and elem_count (inclusive).
 */
unsigned            p4est_quadrant_checksum (sc_array_t * quadrants,
                                             sc_array_t * checkarray,
                                             size_t first_quadrant);

/** Report whether a quadrant fits into the limits of a quadrant range.
 * The range's boundaries are determined by its first and last descendant.
 * Such descendants are for example stored in each tree of a p4est_t.
 * \param [in] fd           First descendant quadrant of a range.
 *                          Thus its level must be P4EST_QMAXLEVEL.
 *                          Must be valid, its data is ignored.
 * \param [in] ld           Last descendant quadrant of a range.
 *                          Thus it must not be smaller than \a fd
 *                          and its level must be P4EST_QMAXLEVEL.
 *                          Must be valid, its data is ignored.
 * \param [in] quadrant     Quadrant checked to be fully inside the range.
 *                          Must at least be extended, its data is ignored.
 * \return                  True if and only if quadrant is contained in range.
 * \see p4est_quadrant_is_extended
 */
int                 p4est_quadrant_in_range (const p4est_quadrant_t * fd,
                                             const p4est_quadrant_t * ld,
                                             const p4est_quadrant_t *
                                             quadrant);

/** Test if a tree is sorted in Morton ordering.
 * \return Returns true if sorted, false otherwise.
 * \note Duplicate quadrants are not allowed.
 */
int                 p4est_tree_is_sorted (p4est_tree_t * tree);

/** Test if a tree is sorted in Morton ordering and linear.
 * \return Returns true if linear, false otherwise.
 * \note Linear means that the tree has no overlaps.
 */
int                 p4est_tree_is_linear (p4est_tree_t * tree);

/** Test if a tree is sorted in Morton ordering and complete.
 * \return Returns true if complete, false otherwise.
 * \note Complete means that the tree has no holes and no overlaps.
 */
int                 p4est_tree_is_complete (p4est_tree_t * tree);

/** Check if a tree is sorted/linear except for diagonally outside corners.
 * \param [in]  check_linearity  Boolean for additional check for linearity.
 * \return Returns true if almost sorted/linear, false otherwise.
 */
int                 p4est_tree_is_almost_sorted (p4est_tree_t * tree,
                                                 int check_linearity);

/** Print the quadrants in a tree.
 * Prints one line per quadrant with x, y, level and a string.
 * The string denotes the relation to the previous quadrant and can be:
 *   Fn  for the first quadrant in the tree with child id n
 *   I   for identical quadrants
 *   R   for a quadrant that with smaller Morton index
 *   Cn  for child id n
 *   Sn  for sibling with child id n
 *   D   for a descendant
 *   Nn   for a next quadrant in the tree with no holes in between and child id n
 *   qn  for a general quadrant whose child id is n
 * \param [in] tree        Any (possibly incomplete, unsorted) tree to be printed.
 */
void                p4est_tree_print (int log_priority, p4est_tree_t * tree);

/** Locally check forest/connectivity structures for equality.
 * \param [in] p4est1    The first forest to be compared.
 * \param [in] p4est2    The second forest to be compared.
 * \param [in] compare_data     Also check if quadrant data are identical.
 * \return          Returns true if forests and their connectivities are equal.
 */
int                 p4est_is_equal (p4est_t * p4est1, p4est_t * p4est2,
                                    int compare_data);

/** Check a forest for validity and allreduce the result.
 * Some properties of a valid forest are:
 *    the quadrant counters are consistent
 *    all trees are complete
 *    all non-local trees are empty
 * This function is collective!
 * It is also relatively expensive, so its use in production should be limited.
 * \param [in] p4est    The forest to be tested.
 *                      Itself and its connectivity must be non-NULL.
 * \return              Returns true if valid, false otherwise.
 */
int                 p4est_is_valid (p4est_t * p4est);

/** Compute the overlap of a number of insulation layers with a tree.
 * Every quadrant out of the insulation layer of the quadrants in \a in
 * except the quadrant itself is checked for overlap of quadrants
 * from all trees that are smaller by at least two levels and thus
 * can cause a split. Those elements that cause a split (as determined by the
 * p4est_balance_*_test routines) create quadrants in \a out that will
 * reproduce those splits when \a in is balanced.
 * Note: Use this version if you are using less than full balance.
 *
 * \param [in] p4est    The p4est to work on.
 * \param [in] in       A piggy-sorted linear list of quadrants.
 *                      The piggy2->from_tree member must be set.
 * \param [in,out] out  A piggy-sorted subset of tree->quadrants.
 * \param [in] balance  The type of balance condition that should be enforced.
 * \param [in] borders  Array of arrays of tree border elements: if not NULL,
 *                      this will be used to fill \a out.
 * \param [in] inseeds  The seeds that \a in generates locally.
 */
void                p4est_tree_compute_overlap (p4est_t * p4est,
                                                sc_array_t * in,
                                                sc_array_t * out,
                                                p4est_connect_type_t balance,
                                                sc_array_t * borders,
                                                sc_array_t * inseeds);

/** Gets the reduced representation of the overlap that results from using
 * p4est_tree_compute_overlap_new
 * \param [in,out] out  A piggy-sorted subset of tree->quadrants.
  */
void                p4est_tree_uniqify_overlap (sc_array_t * out);

/** Removes quadrants that are outside the owned tree boundaries from a tree.
 * \param [in,out] p4est    The p4est to work on.
 * \param [in] which_tree   Index to a sorted owned tree in the p4est.
 * \return                  Returns the number of removed quadrants.
 */
size_t              p4est_tree_remove_nonowned (p4est_t * p4est,
                                                p4est_topidx_t which_tree);

/** Constructs a minimal linear octree between two octants.
 *
 * This is algorithm 2 from H. Sundar, R.S. Sampath and G. Biros
 * with the additional improvements that we do not require sorting
 * and the runtime is O(N).
 *
 * \pre \a q1 < \a q2 in the Morton ordering.
 *
 * \param [in]  p4est      Used for the memory pools and quadrant init.
 * \param [in]  q1         First input quadrant.  Data init'ed if included.
 * \param [in]  include_q1 Flag to specify whether q1 is included.
 * \param [in]  q2         Second input quadrant.  Data init'ed if included.
 * \param [in]  include_q2 Flag to specify whether q2 is included.
 * \param [out] tree       Initialized tree with zero elements.
 * \param [in]  which_tree The 0-based index of \a tree which is needed for
 *                         the \c p4est_quadrant_init_data routine.
 * \param [in]  init_fn    Callback function to initialize the user_data
 *                         which is already allocated automatically.
 */
void                p4est_complete_region (p4est_t * p4est,
                                           const p4est_quadrant_t * q1,
                                           int include_q1,
                                           const p4est_quadrant_t * q2,
                                           int include_q2,
                                           p4est_tree_t * tree,
                                           p4est_topidx_t which_tree,
                                           p4est_init_t init_fn);

/** Completes a sorted tree within a p4est. It may have exterior quadrants.
 * The completed tree will have only owned quadrants and no overlap.
 * Note that the tree's counters (\a quadrants_per_level, \a maxlevel) must be
 * correct for the quadrants in the incoming tree.
 *
 * \param [in,out] p4est      The p4est to work on.
 * \param [in]     which_tree The 0-based index of the subtree to complete.
 * \param [in]     init_fn    Callback function to initialize the user_data
 *                            which is already allocated automatically.
 */
void                p4est_complete_subtree (p4est_t * p4est,
                                            p4est_topidx_t which_tree,
                                            p4est_init_t init_fn);

/** Balances a sorted tree within a p4est. It may have exterior quadrants.
 * The completed tree will have only owned quadrants and no overlap.
 * \param [in,out] p4est      The p4est to work on.
 * \param [in]     btype      The balance type (face or corner).
 * \param [in]     which_tree The 0-based index of the subtree to balance.
 * \param [in]     init_fn    Callback function to initialize the user_data
 *                            which is already allocated automatically.
 */
void                p4est_balance_subtree (p4est_t * p4est,
                                           p4est_connect_type_t btype,
                                           p4est_topidx_t which_tree,
                                           p4est_init_t init_fn);

void                p4est_balance_border (p4est_t * p4est,
                                          p4est_connect_type_t btype,
                                          p4est_topidx_t which_tree,
                                          p4est_init_t init_fn,
                                          p4est_replace_t replace_fn,
                                          sc_array_t * borders);

/** Remove overlaps from a sorted list of quadrants.
 *
 * This is algorithm 8 from H. Sundar, R.S. Sampath and G. Biros
 * with the additional improvement that it works in-place.
 *
 * \param [in]     p4est used for the memory pool and quadrant free.
 * \param [in,out] tree   A sorted tree to be linearized in-place.
 * \return                Returns the number of removed quadrants.
 */
size_t              p4est_linearize_tree (p4est_t * p4est,
                                          p4est_tree_t * tree);

/** Compute correction of partition for a process.
 *
 * The correction denotes how many quadrants the process with id \a rank takes
 * from (if correction positive) or gives to (if correction negative) the
 * previous process with id \a rank-1 in order to assign a family of quadrants
 * to one process.
 * The process with the highest number of quadrants of a family gets all
 * quadrants belonging to this family from other processes. If this applies to
 * several processes, then the process with the lowest id gets the quadrants.
 * A process can give more quadrants than it owns, if it passes quadrants from
 * other processes.
 *
 * \param [in] partition       first global quadrant index for each process (+1)
 * \param [in] num_procs       number of processes
 * \param [in] rank            process id for which correction is computed
 * \param [in] min_quadrant_id minimal global quadrant index of family
 * \param [in] max_quadrant_id maximal global quadrant index of family
 * \return                     correction for process \a rank
 */
p4est_locidx_t      p4est_partition_correction (p4est_gloidx_t *
                                                partition,
                                                int num_procs,
                                                int rank,
                                                p4est_gloidx_t
                                                min_quadrant_id,
                                                p4est_gloidx_t
                                                max_quadrant_id);

/** Correct partition counters to allow one level of coarsening.
 * No quadrants are actually shipped, just the desired number is updated.
 * This function guarantees that empty processors remain empty.
 * This function is collective and acts as a synchronization point.
 *
 * \param [in] p4est                     forest whose partition is corrected
 * \param [in,out] num_quadrants_in_proc partition that will be corrected
 * \return                               total absolute number of moved
 *                                       quadrants.  In practice, at most
 *                                       a small number per processor.
 */
p4est_gloidx_t      p4est_partition_for_coarsening (p4est_t * p4est,
                                                    p4est_locidx_t *
                                                    num_quadrants_in_proc);

/** Find next non-empty process.
 *
 * Finds the next process id >= \a rank which is not empty according to
 * \a num_quadrants_in_proc.
 *
 * \param [in] rank                  process id where search starts
 * \param [in] num_proc              number of processes
 * \param [in] num_quadrants_in_proc number of quadrants for each process
 * \return                           process id of a non empty process
 */
int                 p4est_next_nonempty_process (int rank,
                                                 int num_procs,
                                                 p4est_locidx_t *
                                                 num_quadrants_in_proc);

/** Partition \a p4est given the number of quadrants per proc.
 *
 * Given the desired number of quadrants per proc \a num_quadrants_in_proc
 * the forest \a p4est is partitioned.
 *
 * \param [in,out] p4est the forest that is partitioned.
 * \param [in]     num_quadrants_in_proc  an integer array of the number of
 *                                        quadrants desired per processor.
 * \return  Returns the global count of shipped quadrants.
 */
p4est_gloidx_t      p4est_partition_given (p4est_t * p4est,
                                           const p4est_locidx_t *
                                           num_quadrants_in_proc);

SC_EXTERN_C_END;

#endif /* !P4EST_ALGORITHMS_H */
