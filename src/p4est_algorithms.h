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

#ifndef P4EST_ALGORITHMS_H
#define P4EST_ALGORITHMS_H

#include <p4est.h>

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

/** Test if a tree is sorted in Morton ordering.
 * \return Returns true if sorted, false otherwise.
 * \note Duplicate quadrants are not allowed.
 */
bool                p4est_tree_is_sorted (p4est_tree_t * tree);

/** Test if a tree is sorted in Morton ordering and linear.
 * \return Returns true if linear, false otherwise.
 * \note Linear means that the tree has no overlaps.
 */
bool                p4est_tree_is_linear (p4est_tree_t * tree);

/** Test if a tree is sorted in Morton ordering and complete.
 * \return Returns true if complete, false otherwise.
 * \note Complete means that the tree has no holes and no overlaps.
 */
bool                p4est_tree_is_complete (p4est_tree_t * tree);

/** Check if a tree is sorted/linear except for diagonally outside corners.
 * \param [in]  check_linearity  Boolean for additional check for linearity.
 * \return Returns true if almost sorted/linear, false otherwise.
 */
bool                p4est_tree_is_almost_sorted (p4est_tree_t * tree,
                                                 bool check_linearity);

/** Print the quadrants in a tree.
 * Prints one line per quadrant with x, y, level and a string.
 * The string denotes the relation to the previous quadrant and can be:
 *   Fn  for the first quadrant in the tree with child id n
 *   I   for identical quadrants
 *   R   for a quadrant that with smaller Morton index
 *   Cn  for child id n
 *   Sn  for sibling with child id n
 *   D   for a descendent
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
bool                p4est_is_equal (p4est_t * p4est1, p4est_t * p4est2,
                                    bool compare_data);

/** Check a forest for validity and allreduce the result.
 * Some properties of a valid forest are:
 *    the quadrant counters are consistent
 *    all trees are complete
 *    all non-local trees are empty
 * \param [in] p4est    The forest to be tested.
 * \return              Returns true if valid, false otherwise.
 */
bool                p4est_is_valid (p4est_t * p4est);

/** Find the lowest position tq in a quadrant array such that tq >= q.
 * \return  Returns the id of the matching quadrant
 *                  or -1 if not found or the array is empty.
 */
ssize_t             p4est_find_lower_bound (sc_array_t * array,
                                            const p4est_quadrant_t * q,
                                            size_t guess);

/** Given a sorted \a array of quadrants that have a common ancestor at level
 * \a level, compute the \a indices of the first quadrant in each of the common
 * ancestor's children at level \a level + 1;
 * \param [in] array     The sorted array of quadrants of level > \a level.
 * \param [in] level     The level at which there is a common ancestor.
 * \param [in,out] indices     The indices of the first quadrant in each of
 *                             the ancestors's children, plus an additional
 *                             index on the end.  The quadrants of \a array
 *                             that are descendents of child i have indices
 *                             between indices[i] and indices[i + 1] - 1.  If
 *                             indices[i] = indices[i+1], this indicates that
 *                             no quadrant in the array is contained in
 *                             child i.
 */
void                p4est_split_array (sc_array_t * array, int level,
                                       size_t indices[]);

/** Find the highest position tq in a quadrant array such that tq <= q.
 * \return  Returns the id of the matching quadrant
 *                  or -1 if not found or the array is empty.
 */
ssize_t             p4est_find_higher_bound (sc_array_t * array,
                                             const p4est_quadrant_t * q,
                                             size_t guess);

/** Compute the overlap of a number of insulation layers with a tree.
 * Every quadrant out of the insulation layer of the quadrants in \a in
 * except the quadrant itself is checked for overlap of quadrants
 * from all trees that are smaller by at least two levels and thus
 * can cause a split. Those elements that overlap are put into \a out.
 * \param [in] p4est    The p4est to work on.
 * \param [in] in       A piggy-sorted linear list of quadrants.
 * \param [in,out] out  A piggy-sorted subset of tree->quadrants.
 */
void                p4est_tree_compute_overlap (p4est_t * p4est,
                                                sc_array_t * in,
                                                sc_array_t * out);

/** Removes duplicate quadrant from the output array of compute_overlap.
 * \param [in] skip     A piggy-sorted list of quadrants to be skipped.
 * \param [in,out] out  A piggy-sorted subset of tree->quadrants.
  */
void                p4est_tree_uniqify_overlap (sc_array_t * skip,
                                                sc_array_t * out);

/** Removes quadrants that are outside the owned tree boundaries from a tree.
 * \param [in,out] p4est    The p4est to work on.
 * \param [in] which_tree   Index to a sorted owned tree in the p4est.
 * \return                  Returns the number of removed quadrants.
 */
size_t              p4est_tree_remove_nonowned (p4est_t * p4est,
                                                p4est_topidx_t which_tree);

/** Constructs a minimal linear octree between two octants.
 *
 * This is alogorithm 2 from H. Sundar, R.S. Sampath and G. Biros
 * with the additional improvements that we do not require sorting
 * and the runtime is O(N).
 *
 * \pre \a q1 < \a q2 in the Morton ordering.
 *
 * \param [in]  p4est      Used for the memory pools and quadrant init.
 * \param [in]  q1         First input quadrant.  User data will not change.
 * \param [in]  include_q1 Flag is set to true if q1 is included.
 * \param [in]  q2         Second input quadrant.  User data will not change.
 * \param [in]  include_q2 Flag is set to true if q2 is included.
 * \param [out] tree       Initialized tree with zero elements.
 * \param [in]  which_tree The 0-based index of \a tree which is needed for
 *                         the \c p4est_quadrant_init_data routine.
 * \param [in]  init_fn    Callback function to initialize the user_data
 *                         which is already allocated automatically.
 */
void                p4est_complete_region (p4est_t * p4est,
                                           const p4est_quadrant_t * q1,
                                           bool include_q1,
                                           const p4est_quadrant_t * q2,
                                           bool include_q2,
                                           p4est_tree_t * tree,
                                           p4est_topidx_t which_tree,
                                           p4est_init_t init_fn);

/** Completes a sorted tree within a p4est. It may have exterior quadrants.
 * The completed tree will have only owned quadrants and no overlap.
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
                                           p4est_balance_type_t btype,
                                           p4est_topidx_t which_tree,
                                           p4est_init_t init_fn);

/** Remove overlaps from a sorted list of quadrants.
 *
 * This is alogorithm 8 from H. Sundar, R.S. Sampath and G. Biros
 * with the additional improvement that it works in-place.
 *
 * \param [in]     p4est used for the memory pool and quadrant free.
 * \param [in,out] tree   A sorted tree to be linearized in-place.
 * \return                Returns the number of removed quadrants.
 */
size_t              p4est_linearize_tree (p4est_t * p4est,
                                          p4est_tree_t * tree);

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
