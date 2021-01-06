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
 * This file provides several helper functions and a couple highlevel
 * recursive search algorithms.  These can be used to search for a collection
 * of user-defined "points" through the forest.  There are three flavors of
 * the main search algorithm:
 *
 * 1. \ref p8est_search_local
 *
 *    This function examines the processor-local part of the refinement tree.
 *    It proceeds top-down along all subtrees that have at least one local leaf.
 *    Non-local subtrees are ignored in an optimized way.  Use this function to
 *    compare points against the local branches and leaves.
 *
 * 2. \ref p8est_search_partition
 *
 *    This function examines the parallel partition that is known to all
 *    processors without knowing about actual leaves on remote processors.
 *    Use this to find the processors relevant for a collection of points,
 *    which can then be used to send each point to its assigned processor.
 *
 * 3. \ref p8est_search_all
 *
 *    This function combines the first two into one algorithm.  Note that when
 *    the parallel partition is not of interest, \ref p8est_search_local is
 *    recommended instead since it employs optimizations that are only possible
 *    on the local processor.
 *
 * \ingroup p8est
 */

#include <p8est.h>

SC_EXTERN_C_BEGIN;

/** Binary search in partition array.
 * Given two targets \a my_begin and \a my_end, find offsets such that
 * `search_in[begin] >= my_begin`, `my_end <= search_in[end]`.
 * If more than one index satisfies the conditions, then the minimal index is the
 * result. If there is no index that satisfies the conditions, then \a begin
 * and \a end are tried to set equal such that `search_in[begin] >= my_end`.
 * If \a my_begin is less or equal than the smallest value of \a search_in
 * \a begin is set to 0 and if \a my_end is bigger or equal than the largest
 * value of \a search_in \a end is set to \a num_procs - 1.
 * If none of the above conditions is satisfied, the output is not well defined.
 * We require `my_begin <= my_begin'.
 * \param [in] num_procs    Number of processes to get the length of \a search_in.
 * \param [in] search_in    The sorted array (ascending) in that the function will search.
 *                          If `k` indexes search_in, then `0 <= k < num_procs`.
 * \param [in] my_begin     The first target that defines the start of the search window.
 * \param [in] my_end       The second target that defines the end of the search window.
 * \param [in,out] begin    The first offset such that `search_in[begin] >= my_begin`.
 * \param [in,out] end      The second offset such that `my_end <= search_in[end]`.
 */
void                p8est_find_partition (const int num_procs,
                                          p4est_gloidx_t * search_in,
                                          p4est_gloidx_t my_begin,
                                          p4est_gloidx_t my_end,
                                          p4est_gloidx_t * begin,
                                          p4est_gloidx_t * end);

/** Find the lowest position tq in a quadrant array such that tq >= q.
 * \return  Returns the id of the matching quadrant
 *                  or -1 if array < q or the array is empty.
 */
ssize_t             p8est_find_lower_bound (sc_array_t * array,
                                            const p8est_quadrant_t * q,
                                            size_t guess);

/** Find the highest position tq in a quadrant array such that tq <= q.
 * \return  Returns the id of the matching quadrant
 *                  or -1 if array > q or the array is empty.
 */
ssize_t             p8est_find_higher_bound (sc_array_t * array,
                                             const p8est_quadrant_t * q,
                                             size_t guess);

/** Search a local quadrant by its cumulative number in the forest.
 *
 * We perform a binary search over the processor-local trees,
 * which means that it is advisable NOT to use this function if possible,
 * and to try to maintain O(1) tree context information in the calling code.
 *
 * \param [in]  p8est           Forest to be worked with.
 * \param [in]  cumulative_id   Cumulative index over all trees of quadrant.
 * \param [in,out] which_tree   If not NULL, the input value can be -1
 *                              or an initial guess for the quadrant's tree.
 *                              An initial guess must be the index of a
 *                              nonempty local tree.
 *                              Output is the tree of returned quadrant.
 * \param [out] quadrant_id     If not NULL, the number of quadrant in tree.
 * \return                      The identified quadrant.
 */
p8est_quadrant_t   *p8est_find_quadrant_cumulative (p8est_t * p8est,
                                                    p4est_locidx_t
                                                    cumulative_id,
                                                    p4est_topidx_t *
                                                    which_tree,
                                                    p4est_locidx_t *
                                                    quadrant_id);

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
 * \param [in] p4est        The forest to be queried.
 * \param [in] which_tree   The tree id under consideration.
 * \param [in] quadrant     The quadrant under consideration.
 *                          This quadrant may be coarser than the quadrants
 *                          that are contained in the forest (an ancestor), in
 *                          which case it is a temporary variable and not part
 *                          of the forest storage.  Otherwise, it is a leaf and
 *                          points directly into the forest storage.
 * \param [in] local_num    If the quadrant is not a leaf, this is < 0.
 *                          Otherwise it is the (non-negative) index of the
 *                          quadrant relative to the processor-local storage.
 * \param [in] point        Representation of a "point"; user-defined.
 *                          If \b point is NULL, the callback may be used to
 *                          prepare quadrant-related search meta data.
 * \return                  If \b point is NULL, true if the search confined to
 *                          \b quadrant should be executed, false to skip it.
 *                          Else, true if point may be contained in the
 *                          quadrant and false otherwise; the return value has
 *                          no effect on a leaf.
 */
typedef int         (*p8est_search_local_t) (p8est_t * p4est,
                                             p4est_topidx_t which_tree,
                                             p8est_quadrant_t * quadrant,
                                             p4est_locidx_t local_num,
                                             void *point);

/** This typedef is provided for backwards compatibility. */
typedef p8est_search_local_t p8est_search_query_t;

/** Search through the local part of a forest.
 * The search is especially efficient if multiple targets, called "points"
 * below, are searched for simultaneously.
 *
 * The search runs over all local quadrants and proceeds recursively top-down.
 * For each tree, it may start at the root of that tree, or further down at the
 * root of the subtree that contains all of the tree's local quadrants.
 * Likewise, some intermediate levels in the recursion may be skipped if the
 * processor-local part is contained in a single deeper subtree.
 * The outer loop is thus a depth-first, processor-local forest traversal.
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
 * in the current branch, otherwise it is passed to the next deeper level.
 * The callback is allowed to return true for the same point and more than one
 * quadrant; in this case more than one matching quadrant may be identified.
 * The callback is also allowed to return false for all children of a quadrant
 * that it returned true for earlier.  If the point callback returns false for
 * all points relevant to a quadrant, the recursion stops.
 * The points can really be anything, p4est does not perform any
 * interpretation, just passes the pointer along to the callback function.
 *
 * If points are present and the first quadrant callback returned true, we
 * execute it a second time after calling the point callback for all current
 * points.  This can be used to gather and postprocess information about the
 * points more easily.  If it returns false, the recursion stops.
 *
 * If the points are a NULL array, they are ignored and the recursion proceeds
 * by querying the per-quadrant callback.  If the points are not NULL but an
 * empty array, the recursion will stop immediately!
 *
 * \param [in] p4est        The forest to be searched.
 * \param [in] call_post    If true, call quadrant callback both pre and post.
 * \param [in] quadrant_fn  Executed once when a quadrant is entered, and once
 *                          when it is left (the second time only if points are
 *                          present and the first call returned true).
 *                          This quadrant is always local, if not completely
 *                          than at least one descendant of it.  If the
 *                          callback returns false, this quadrant and its
 *                          descendants are excluded from the search recursion.
 *                          Its \b point argument is always NULL.
 *                          Callback may be NULL in which case it is ignored.
 * \param [in] point_fn     If \b points is not NULL, must be not NULL.
 *                          Shall return true for any possible matching point.
 *                          If \b points is NULL, this callback is ignored.
 * \param [in] points       User-defined array of "points".
 *                          If NULL, only the \b quadrant_fn callback
 *                          is executed.  If that is NULL, this function noops.
 *                          If not NULL, the \b point_fn is called on
 *                          its members during the search.
 */
void                p8est_search_local (p8est_t * p4est, int call_post,
                                        p8est_search_local_t quadrant_fn,
                                        p8est_search_local_t point_fn,
                                        sc_array_t * points);

/** This function is provided for backwards compatibility.
 * We call \ref p8est_search_local with call_post = 0.
 */
void                p8est_search (p8est_t * p4est,
                                  p8est_search_query_t quadrant_fn,
                                  p8est_search_query_t point_fn,
                                  sc_array_t * points);

/** Callback function for the partition recursion.
 * \param [in] p4est        The forest to traverse.
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
 *                          \b quadrant's branch after this function returns.
 * \param [in,out] point    Pointer to a user-defined point object.
 *                          If called per-quadrant, this is NULL.
 * \return                  If false, the recursion at quadrant is terminated.
 *                          If true, it continues if \b pfirst < \b plast.
 */
typedef int         (*p8est_search_partition_t) (p8est_t * p4est,
                                                 p4est_topidx_t which_tree,
                                                 p8est_quadrant_t * quadrant,
                                                 int pfirst, int plast,
                                                 void *point);

/** Traverse the global partition top-down.
 * We proceed top-down through the partition, identically on all processors
 * except for the results of two user-provided callbacks.  The recursion will only
 * go down branches that are split between multiple processors.  The callback
 * functions can be used to stop a branch recursion even for split branches.
 * This function offers the option to search for arbitrary user-defined points
 * analogously to \ref p4est_search_local.
 * \note Traversing the whole processor partition will be at least O(P),
 *       so sensible use of the callback function is advised to cut it short.
 * \param [in] p4est        The forest to traverse.
 *                          Its local quadrants are never accessed.
 * \param [in] call_post    If true, call quadrant callback both pre and post.
 * \param [in] quadrant_fn  This function controls the recursion,
 *                          which only continues deeper if this
 *                          callback returns true for a branch quadrant.
 *                          It is allowed to set this to NULL.
 * \param [in] point_fn     This function decides per-point whether it is
 *                          followed down the recursion.
 *                          Must be non-NULL if \b points are not NULL.
 * \param [in] points       User-provided array of \b points that are
 *                          passed to the callback \b point_fn.
 *                          See \ref p8est_search_local for details.
 */
void                p8est_search_partition (p8est_t * p4est, int call_post,
                                            p8est_search_partition_t
                                            quadrant_fn,
                                            p8est_search_partition_t point_fn,
                                            sc_array_t * points);

/** Callback function for the top-down search through the whole forest.
 * \param [in] p4est        The forest to search.
 *                          We recurse through the trees one after another.
 * \param [in] which_tree   The current tree number.
 * \param [in] quadrant     The current quadrant in the recursion.
 *                          This quadrant is either a non-leaf tree branch
 *                          or a leaf.  If the quadrant is contained in the
 *                          local partition, we know which, otherwise we don't.
 *                          Let us first consider the situation when \b
 *                          quadrant is local, which is indicated by both \b
 *                          pfirst and \b plast being equal to \b
 *                          p4est->mpirank.  Then the parameter \b local_num is
 *                          negative for non-leaves and the number of the
 *                          quadrant as a leaf in local storage otherwise.
 *                          Only if the quadrant is a local leaf, it points to
 *                          the actual local storage and can be used to access
 *                          user data etc., and the recursion terminates.
 *                          The other possibility is that \b pfirst < \b plast,
 *                          in which case we proceed with the recursion,
 *                          or both are equal to the same remote rank, in
 *                          which case the recursion terminates.  Either way,
 *                          the quadrant is not from local forest storage.
 *
 * \param [in] pfirst       The lowest processor that owns part of \b quadrant.
 *                          Guaranteed to be non-empty.
 * \param [in] plast        The highest processor that owns part of \b quadrant.
 *                          Guaranteed to be non-empty.
 * \param [in] local_num    If \b quadrant is a local leaf, this number is the
 *                          index of the leaf in local quadrant storage.
 *                          Else, this is a negative value.
 *
 * \param [in,out] point    User-defined representation of a point.  This
 *                          parameter distinguishes two uses of the callback.
 *                          For each quadrant, the callback is first called
 *                          with a NULL point, and if this callback returns
 *                          true, once for each point tracked in this branch.
 *                          The return value for a point determines whether
 *                          it shall be tracked further down the branch or not,
 *                          and has no effect on a local leaf.
 *                          The call with a NULL point is intended to prepare
 *                          quadrant-related search meta data that is common to
 *                          all points, and/or to efficiently terminate the
 *                          recursion for all points in the branch in one call.
 *
 * \return                  If false, the recursion at \b quadrant terminates.
 *                          If true, it continues if \b pfirst < \b plast or
 *                          if they are both equal to \b p4est->mpirank and
 *                          the recursion has not reached a leaf yet.
 */
typedef int         (*p8est_search_all_t) (p8est_t * p8est,
                                           p4est_topidx_t which_tree,
                                           p8est_quadrant_t * quadrant,
                                           int pfirst, int plast,
                                           p4est_locidx_t local_num,
                                           void *point);

/** Perform a top-down search on the whole forest.
 *
 * This function combines the functionality of \ref p4est_search_local and \ref
 * p4est_search_partition; their documentation applies for the most part.
 *
 * The recursion proceeds from the root quadrant of each tree until
 * (a) we encounter a remote quadrant that covers only one processor, or
 * (b) we encounter a local leaf quadrant.
 * In other words, we proceed with the recursion into a quadrant's children if
 * (a) the quadrant is split between two or more processors, no matter whether
 * one of them is the calling processor or not, or (b) if the quadrant is on
 * the local processor but we have not reached a leaf yet.
 *
 * The search can track one or more points, which are abstract placeholders.
 * They are matched against the quadrants traversed using a callback function.
 * The result of the callback function can be used to stop a recursion early.
 * The user determines how a point is interpreted, we only pass it around.
 *
 * Note that in the remote case (a), we may terminate the recursion even if
 * the quadrant is not a leaf, which we have no means of knowing.  Still,
 * this case is sufficient to determine the processor ownership of a point.
 *
 * \note
 * This is a very powerful function that can become slow if not used carefully.
 *
 * \note
 * As with the two other search functions in this file, calling it once with
 * many points is generally much faster than calling it once for each point.
 * Using multiple points also allows for a per-quadrant termination of the
 * recursion in addition to a more costly per-point termination.
 *
 * \note
 * This function works fine when used for the special cases that either the
 * partition or the local quadrants are not of interest.  However, in the case
 * of querying only local information we expect that \ref p4est_search_local
 * will be faster since it employs specific local optimizations.
 *
 * \param [in] p4est        The forest to be searched.
 * \param [in] call_post    If true, call quadrant callback both pre and post.
 * \param [in] quadrant_fn  Executed once for each quadrant that is entered.
 *                          If the callback returns false, this quadrant and
 *                          its descendants are excluded from the search, and
 *                          the points in this branch are not queried further.
 *                          Its \b point argument is always NULL.
 *                          Callback may be NULL in which case it is ignored.
 * \param [in] point_fn     Executed once for each point that is relevant for a
 *                          quadrant of the search.  If it returns true, the
 *                          point is tracked further down that branch, else it
 *                          is discarded from the queries for the children.
 *                          If \b points is not NULL, this callback must be not
 *                          NULL.  If \b points is NULL, it is not called.
 * \param [in] points       User-defined array of points.  We do not interpret
 *                          a point, just pass it into the callbacks.
 *                          If NULL, only the \b quadrant_fn callback
 *                          is executed.  If that is NULL, the whole function
 *                          noops.  If not NULL, the \b point_fn is
 *                          called on its members during the search.
 */
void                p8est_search_all (p8est_t * p4est, int call_post,
                                      p8est_search_all_t quadrant_fn,
                                      p8est_search_all_t point_fn,
                                      sc_array_t * points);

SC_EXTERN_C_END;

#endif /* !P8EST_SEARCH_H */
