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

#include <sc.h>

#ifndef P8EST_MULTI_OVERSET_H
#define P8EST_MULTI_OVERSET_H

SC_EXTERN_C_BEGIN;

/** Return true if \a p8est is an artifical p8est.
 *
 * This function can be used in \ref p8est_intersect_t callbacks to
 * distinguish temporary, artificial p8est instances that were created for a
 * local partition search (and only contain some meta information) from a real
 * p8est used for a local search.
 */
int p8est_is_meta (p8est_t *p8est);

/** Callback function for checking if a point intersects a quadrant.
 *
 * This function can be passed to \ref p8est_multi_overset to search for all
 * pairs of a quadrant from the background p8est and a query point from a
 * overset mesh, which is contained in the quadrant.
 * It will be called both in a \ref p8est_search_partition and a
 * \ref p8est_search_local.
 * Use \ref p8est_is_meta, to determine which is the case.
 * \param [in] p8est            The p8est we search in.
 *                              This is either a valid p8est used in a local
 *                              search or a temporary, artificial p8est created
 *                              for a partition search. In the latter case most
 *                              entries will be initialized to -1. Its member
 *                              mpirank will be set, whenever we are on a leaf
 *                              of the partition search tree.
 * \param [in] which_tree       The current tree number.
 * \param [in] quadrant         The quadrant under consideration.
 *                              When on a leaf of a valid p8est, this is a
 *                              quadrant from the local forest storage of the
 *                              valid p8est.
 *                              Otherwise, this is a temporary quadrant that is
 *                              not from local forest storage. It is a valid
 *                              quadrant with user data set to NULL. It
 *                              represents an ancestor of the leaves of the
 *                              forest in the top-down recursion.
 * \param [in] local_num        If local_num is -1, we are on a temporary quadrant.
 *                              Otherwise, this is a valid local_num indexing
 *                              into local quadrant storage of \a p8est.
 * \param [in] point            An abstract user-defined point.
 * \param [in] user             Arbitrary user data passed on from the function
 *                              calling this callback.
 * \return                      True, if there is a possible intersection of the
 *                              quadrant or one of its ancestors and the point.
 *                              Return false if there is definitely no
 *                              intersection.
 *                              If we are on a leaf of the local search
 *                              (local_num is non-negative) or the partition
 *                              search (p8est_is_meta and the mpirank of
 *                              \a p8est is non-negative.) this callback should
 *                              do an exact test for intersection.
 *                              Else, the return value may be a false positive,
 *                              we'll be fine.
 */
typedef int (*p8est_intersect_t) (p8est_t *p8est,
                                  p4est_topidx_t which_tree,
                                  p8est_quadrant_t *quadrant,
                                  p4est_locidx_t local_num,
                                  void *point,
                                  void *user);

/** Callback function for interpolation of a point intersecting a quadrant.
 * Functions of this type can be passed to \ref p8est_multi_overset to compute
 * interpolation data for query points based on the background p8est.
 * \param [in] p8est            A valid p8est.
 * \param [in] which_tree       The number of the tree containing \a quadrant.
 * \param [in] quadrant         A leaf quadrant intersecting \a point.
 * \param [in] local_num        The index of the leaf in local quadrant storage.
 * \param [in] point            A pointer to a user-defined query point.
 * \param [out] intpl_data      A pointer to a user-defined interpolation data
 *                              structure. When passing functions of this type
 *                              to \ref p8est_multi_overset, the interpolation
 *                              values should be computed based on \a point and
 *                              the data stored in \a quadrant and \a p8est.
 */
typedef void (*p8est_interpolate_point_t) (p8est_t *p8est,
                                           p4est_topidx_t which_tree,
                                           p8est_quadrant_t *quadrant,
                                           p4est_locidx_t local_num,
                                           void *point,
                                           void *intpl_data);

/** Execute multi-mesh overset algorithm.
 *
 * Exchange interpolation data between a p8est background mesh \a bgp8est and
 * several overset meshes, each residing on its own contiguous subblock of ranks
 * of \a glocomm. Every overset mesh process comes with a local array of
 * query points \a qpoints of form (x,y,z,v) provided by the user, which are
 * searched in \a bgp8est.
 * For each query point the function determines if its a receptor point or a
 * donor point based on its associated volume \a v.
 * A point is a receptor point, if
 *  - it lies on the overset boundary ( \a v is -1)
 *  - the points volume \a v is greater than the volume of the quadrant of
 *    \a bgp8est it is located in and the quadrant does not contain another
 *    query point which lies on the wall boundary ( \a v is -2).
 * The function computes interpolation data for every receptor point and sends
 * it back to the corresponding overset process.
 * \param [in] glocomm          Global communicator over all meshes.
 * \param [in] headcomm         If global rank is first for a mesh, a
 *                              communicator over all first mesh ranks.
 * \param [in] rolecomm         Separate communicator over each mesh.
 * \param [in] myrole           Index of mesh: 0 for background mesh,
 *                              starting from 1 for overset meshes.
 * \param [in] num_meshes       Number of meshes including background.
 * \param [in] mesh_offsets     Array of ascending global ranks,
 *                              one for the first of each mesh, and then
 *                              one more for the end (exclusive of the last).
 * \param [in] bgp8est          For \a myrole zero, the background forest.
 *                              NULL otherwise.
 * \param [in] qpoints          Query points: 4-double-tuples (x, y, z, v).
 *                              \b v is -1 for the overset boundary, -2 for
 *                              the wall boundary, and a non-negative
 *                              representative volume of the point otherwise.
 *                              The points are provided process-local over
 *                              all overset meshes present on this process.
 *                              NULL for \a myrole zero.
 * \param [in] intsc_fn         Callback function used for searching points in
 *                              \a bgp8est (see also \ref p8est_intersect_t).
 *                              This function receives a quadrant and a point
 *                              from \a qpoints. It returns true, if the point
 *                              is contained in the quadrant or one of its
 *                              ancestors, and false otherwise.
 *                              The function will be called both during a
 *                              local search in \a bgp8est and during a
 *                              partition search in \a bgp8est on a overset
 *                              process, where only limited information about
 *                              the forest is available. In the latter case,
 *                              a temporary, artifical quadrant and p8est will
 *                              be passed to the callback.
 *                              It needs to be provided on all processes.
 * \param [out] intpl_data      An array of interpolation data. All elements
 *                              are instances of a user-defined interpolation
 *                              data type. The data is computed by
 *                              calling \a intpl_fn for each receptor point and
 *                              the quadrant containing it.
 *                              The interpolation data is only returned for
 *                              receptor points, not for donor points, for which
 *                              it would be uninitialized. In combination with
 *                              \a intpl_indices the array forms a sparse vector
 *                              representation for the interpolation data of
 *                              \a qpoints.
 *                              NULL for \a myrole zero.
 * \param [out] intpl_indices   An array of size_t entries indexing into
 *                              \a qpoints. The array has as many entries as
 *                              \a intpl_data and relates the interpolation data
 *                              elements to their associated query point.
 *                              NULL for \a myrole zero.
 * \param [in] intpl_fn         Callback function used for interpolating points
 *                              (see also \ref p8est_interpolate_point_t).
 *                              This function receives a leaf quadrant of
 *                              \a bgp8est and a point in \a qpoints for which
 *                              it is supposed to compute interpolation data.
 *                              The function is only called, if the point is
 *                              contained in the quadrant according to
 *                              \a intsc_fn and it was assigned the role of
 *                              receptor point based on its volume \a v.
 *                              The function will receive a pointer to the
 *                              actual \a bgp8est and the corresponding quadrant.
 *                              Thereby, it has access to all user-pointers
 *                              contained in these.
 *                              NULL for \a myrole unequal to zero.
 * \param [in] user             Context data passed to callback.
 */
void                 p8est_multi_overset
  (sc_MPI_Comm glocomm, sc_MPI_Comm headcomm, sc_MPI_Comm rolecomm,
   int myrole, int num_meshes, const int *mesh_offsets,
   p8est_t *bgp8est, sc_array_t *qpoints, p8est_intersect_t intsc_fn,
   sc_array_t *intpl_data, sc_array_t *intpl_indices,
   p8est_interpolate_point_t intpl_fn, void *user);

SC_EXTERN_C_END;

#endif /* !P8EST_MULTI_OVERSET_H */
