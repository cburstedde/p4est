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

#ifndef P4EST_SEARCH_H
#define P4EST_SEARCH_H

/** \file p4est_search.h
 * Search through quadrants, the local part of a forest, or the partition.
 *
 * This file provides several helper functions and recursive algorithms.
 * \ingroup p4est
 */

#include <p4est.h>

SC_EXTERN_C_BEGIN;

/** Find the lowest position tq in a quadrant array such that tq >= q.
 * \return  Returns the id of the matching quadrant
 *                  or -1 if array < q or the array is empty.
 */
ssize_t             p4est_find_lower_bound (sc_array_t * array,
                                            const p4est_quadrant_t * q,
                                            size_t guess);

/** Find the highest position tq in a quadrant array such that tq <= q.
 * \return  Returns the id of the matching quadrant
 *                  or -1 if array > q or the array is empty.
 */
ssize_t             p4est_find_higher_bound (sc_array_t * array,
                                             const p4est_quadrant_t * q,
                                             size_t guess);

/** Split an array of quadrants by the children of an ancestor.
 *
 * Given a sorted \b array of quadrants that have a common ancestor at level
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
void                p4est_split_array (sc_array_t * array, int level,
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
 * \param [in,out] faces       An array of size 4 that is filled: faces[i] is
 *                             true if the range touches that face.
 * \param [in,out] corners     An array of size 4 that is filled: corners[i] is
 *                             true if the range touches that corner.
 *                             \b faces or \b corners may be NULL.
 * \return  Returns an int32_t encoded with the same information in \b faces
 *          and \b corners: the first (least) four bits represent the four
 *          faces, the next four bits represent the four corners.
 */
int32_t             p4est_find_range_boundaries (p4est_quadrant_t * lq,
                                                 p4est_quadrant_t * uq,
                                                 int level, int faces[],
                                                 int corners[]);

/** Callback function to query the match of a "point" with a quadrant.
 *
 * This function can be called in two roles:  Per-quadrant, in which case the
 * parameter \b point is NULL, or per-point, possibly many times per quadrant.
 *
 * \param [in] p4est        The forest to be queried.
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
typedef int         (*p4est_search_query_t) (p4est_t * p4est,
                                             p4est_topidx_t which_tree,
                                             p4est_quadrant_t * quadrant,
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
 * \param [in] p4est        The forest to be searched.
 * \param [in] search_quadrant_fn   Executed once for each quadrant that is
 *                          entered.  This quadrant is always local, if not
 *                          itself then at least one child of it.  If the
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
void                p4est_search (p4est_t * p4est,
                                  p4est_search_query_t search_quadrant_fn,
                                  p4est_search_query_t search_point_fn,
                                  sc_array_t * points);

SC_EXTERN_C_END;

#endif /* !P4EST_SEARCH_H */
