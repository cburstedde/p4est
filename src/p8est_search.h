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

#ifndef P8EST_SEARCH_H
#define P8EST_SEARCH_H

/** \file p8est_search.h
 * Search through quadrants, the local part of a forest, or the partition.
 *
 * This file provides several helper functions and recursive algorithms.
 * \ingroup p8est
 */

#include <p8est.h>

SC_EXTERN_C_BEGIN;

/** Find the lowest position tq in a quadrant array such that tq >= q.
 * \return  Returns the id of the matching quadrant
 *                  or -1 if not found or the array is empty.
 */
ssize_t             p8est_find_lower_bound (sc_array_t * array,
                                            const p8est_quadrant_t * q,
                                            size_t guess);

/** Find the highest position tq in a quadrant array such that tq <= q.
 * \return  Returns the id of the matching quadrant
 *                  or -1 if not found or the array is empty.
 */
ssize_t             p8est_find_higher_bound (sc_array_t * array,
                                             const p8est_quadrant_t * q,
                                             size_t guess);

/** Given a sorted \b array of quadrants that have a common ancestor at level
 * \b level, compute the \b indices of the first quadrant in each of the common
 * ancestor's children at level \b level + 1.
 * \param [in] array     The sorted array of quadrants of level > \b level.
 * \param [in] level     The level at which there is a common ancestor.
 * \param [in,out] indices     The indices of the first quadrant in each of
 *                             the ancestors's children, plus an additional
 *                             index on the end.  The quadrants of \b array
 *                             that are descendants of child i have indices
 *                             between indices[i] and indices[i + 1] - 1.  If
 *                             indices[i] = indices[i+1], this indicates that
 *                             no quadrant in the array is contained in
 *                             child i.
 */
void                p8est_split_array (sc_array_t * array, int level,
                                       size_t indices[]);

/** Find the boundary points touched by a range of quadrants.
 *
 * Given two smallest quadrants, \b lq and \b uq, that mark the first and the
 * last quadrant in a range of quadrants, determine which portions of the tree
 * boundary the range touches.
 * \param [in] lq        The smallest quadrant at the start of the range: if
 *                       NULL, the tree's first quadrant is taken to be the
 *                       start of the range.
 * \param [in] uq        The smallest quadrant at the end of the range: if
 *                       NULL, the tree's last quadrant is taken to be the
 *                       end of the range.
 * \param [in] level     The level of the containing quadrant whose boundaries
 *                       are tested: 0 if we want to test the boundaries of the
 *                       whole tree.
 * \param [in,out] faces       An array of size 6 that is filled: faces[i] is
 *                             true if the range touches that face.
 * \param [in,out] edges       An array of size 12 that is filled: edges[i] is
 *                             true if the range touches that edge.
 * \param [in,out] corners     An array of size 8 that is filled: corners[i] is
 *                             true if the range touches that corner.
 *                             \b faces, \b edges or \b corners may be NULL.
 * \return  Returns an int32_t encoded with the same information in \b faces,
 *          \b edges and \b corners: the first (least) six bits represent the
 *          six faces, the next twelve bits represent the twelve edges, the
 *          next eight bits represent the eight corners.
 */
int32_t             p8est_find_range_boundaries (p8est_quadrant_t * lq,
                                                 p8est_quadrant_t * uq,
                                                 int level, int faces[],
                                                 int edges[], int corners[]);

/** Callback function to query the match of a "point" with a quadrant.
 *
 * This function can be called in two roles:  Per-quadrant, in which case the
 * parameter \b point is NULL, or per-point, possibly many times per quadrant.
 *
 * \param [in] p8est        The forest to be queried.
 * \param [in] which_tree   The tree id under consideration.
 * \param [in] quadrant     The quadrant under consideration.
 *                          This quadrant may be coarser than the quadrants
 *                          that are contained in the forest (an ancestor), in
 *                          which case it is a temporary variable and not part
 *                          of the forest storage.  Otherwise, it is a leaf and
 *                          points directly into the forest storage.
 * \param [in] local_num    If the quadrant is not a leaf, this is -1.  Otherwise
 *                          it is the (non-negative) index of the quadrant
 *                          relative to the processor-local quadrant storage.
 * \param [in] point        Representation of a "point"; user-defined.
 *                          If \b point is NULL, the callback may be used to
 *                          prepare quadrant-related search meta data.
 * \return                  If \b point is NULL, true if the search confined to
 *                          \b quadrant should be executed, false to skip it.
 *                          Else, true if point may be contained in the
 *                          quadrant and false otherwise; the return value has
 *                          no effect on a leaf.
 */
typedef int         (*p8est_search_query_t) (p8est_t * p8est,
                                             p4est_topidx_t which_tree,
                                             p8est_quadrant_t * quadrant,
                                             p4est_locidx_t local_num,
                                             void *point);

/** Search through the local part of a forest.
 * The search is especially efficient if multiple targets, called "points"
 * below, are searched for simultaneously.
 *
 * The search runs over all local quadrants and proceeds recursively top-down.
 * For each tree, it may start at the root of that tree, or further down at the
 * root of the subtree that contains all of the tree's local quadrants.
 * Likewise, some intermediate levels in the recursion may be skipped.
 * Its outer loop is thus a depth-first, processor-local forest traversal.
 * Each quadrant in that loop either is a leaf, or a (direct or indirect)
 * strict ancestor of a leaf.  On entering a new quadrant, a user-provided
 * quadrant-callback is executed.
 *
 * As a convenience, the user may provide anonymous "points" that are tracked
 * down the forest.  This way one search call may be used for multiple targets.
 * The set of points that potentially matches a given quadrant diminishes from
 * the root down to the leaves:  For each quadrant, an inner loop over the
 * potentially matching points executes a point-callback for each candidate
 * that determines whether the point may be a match.  If not, it is discarded
 * immediately, otherwise it is passed to the next finer level.
 * The callback is allowed to return true for the same point and more than one
 * quadrant; in this case more than one matching quadrant may be identified.
 * The callback is also allowed to return false for all children of a quadrant
 * that it returned true for earlier.
 * The points can really be anything, p4est does not perform any
 * interpretation, just passes the pointer along to the callback function.
 *
 * \param [in] p8est        The forest to be searched.
 * \param [in] search_quadrant_fn   Executed once for each quadrant that is
 *                          entered.  This quadrant is always local, if not
 *                          itself than at least one child of it.  If the
 *                          callback returns false, this quadrant and its
 *                          descendants are excluded from the search.
 *                          Its \b point argument is always NULL.
 *                          May be NULL in which case it is ignored.
 * \param [in] search_point_fn      If \b points is not NULL, must be not NULL.
 *                          Must return true for any possible matching point.
 *                          If \b points is NULL, this callback is ignored.
 * \param [in] points       User-defined array of "points".
 *                          If NULL, only the \b search_quadrant_fn callback
 *                          is executed.  If that is NULL, this function noops.
 *                          If not NULL, the \b search_point_fn is called on
 *                          its members during the search.
 */
void                p8est_search (p8est_t * p8est,
                                  p8est_search_query_t search_quadrant_fn,
                                  p8est_search_query_t search_point_fn,
                                  sc_array_t * points);

/** Callback function for the traversal recursion.
 * \param [in] p8est        The forest to traverse.
 *                          Its local quadrants are never accessed.
 * \param [in] which_tree   The tree number under consideration.
 * \param [in] quadrant     This quadrant is not from local forest storage,
 *                          and its user data is undefined.  It represents
 *                          the branch of the forest in the top-down recursion.
 * \param [in] pfirst       The lowest processor that owns part of \b quadrant.
 *                          Guaranteed to be non-empty.
 * \param [in] plast        The highest processor that owns part of \b quadrant.
 *                          Guaranteed to be non-empty.  If this is equal to
 *                          \b pfirst, then the recursion will stop for
 *                          quadrant's branch after this function returns.
 * \return                  If false, the recursion at quadrant is terminated.
 *                          If true, it continues if \b pfirst < \b plast.
 */
typedef int         (*p8est_traverse_query_t) (p8est_t * p8est,
                                               p4est_topidx_t which_tree,
                                               p8est_quadrant_t * quadrant,
                                               int pfirst, int plast);

/** Traverse the global partition top-down.
 * We proceed top-down through the partition, identically on all processors
 * except for the results of a user-provided callback.  The recursion will only
 * go down branches that are split between multiple processors.  The callback
 * function can be used to stop a branch recursion even for split branches.
 * \note Traversing the whole processor partition will likely by inefficient,
 *       so sensible use of the callback function is advised.
 * \param [in] p8est        The forest to traverse.
 *                          Its local quadrants are never accessed.
 * \param [in] traverse_fn  This function controls the recursion,
 *                          which only continues deeper if this
 *                          callback returns true for a branch quadrant.
 */
void                p8est_traverse (p8est_t * p8est,
                                    p8est_traverse_query_t traverse_fn);

SC_EXTERN_C_END;

#endif /* !P8EST_SEARCH_H */
