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

#ifndef P4EST_BUILD_H
#define P4EST_BUILD_H

#include <p4est.h>

SC_EXTERN_C_BEGIN;

/** \file p4est_build.h
 * Create a new p4est object by adding individual quadrants in order.
 * This can for example be driven by running \ref p4est_search_local.
 * This allows to create a heavily coarsened forest in one pass.
 * It is also legal to add more highly refined quadrants.
 *
 * The only rules are to respect the original partition boundary
 * and to add non-overlapping quadrants in Morton order.
 */

/** Context object for building a new p4est from individual quadrants.
 */
typedef struct p4est_build p4est_build_t;

/** Allocate a context for building a new forest.
 * \param [in] from         This forest is used as a template for creation.
 * \param [in] data_size    Data size of the created forest, may be zero.
 * \param [in] init_fn      This functions is called for created quadrants,
 *                          added manually by \ref p4est_build_add
 *                          or by the internal completion of the subtrees.
 *                          It may be overridden for added quadrants by
 *                          \ref p4est_build_init_add.
 *                          NULL leaves the quadrant data uninitialized.
 * \param [in] user_pointer Registered into the newly built forest.
 * \return                  A context that needs to be processed further.
 */
p4est_build_t      *p4est_build_new (p4est_t * from,
                                     size_t data_size,
                                     p4est_init_t init_fn,
                                     void *user_pointer);

/** Set a dedicated initialization callback for manually added quadrants.
 * \param [in,out] build    The building context at any stage.
 * \param [in] add_init_fn  Henceforth used for quadrants added by
 *                          \ref p4est_build_add.
 *                          NULL leaves the quadrant data uninitialized.
 */
void                p4est_build_init_add (p4est_build_t * build,
                                          p4est_init_t add_init_fn);

/** This function is usable from a \ref p4est_search_local_t callback.
 * It can also be used outside of a search context using proper care.
 *
 * It may be called multiple times in order of trees and then quadrants.
 * The quadrant added in each call must fit entirely into the current tree.
 * This means that inner nodes of the tree may not be legal to pass in here.
 * It is safest to call this function only on leaves of the original tree.
 * However, other calls are possible if subsequent quadrants do not overlap.
 *
 * It is legal to call this function twice with the same quadrant.
 * In this case the second call does nothing.
 *
 * \param [in,out] build    The building context must be passed through.
 * \param [in] which_tree   The tree number is passed from the search callback.
 * \param [in] quadrant     The quadrant is passed from the search callback.
 * \return                  True if the quadrant was added, false if it was
 *                          identical to the previous one and thus not added.
 */
int                 p4est_build_add (p4est_build_t * build,
                                     p4est_topidx_t which_tree,
                                     p4est_quadrant_t * quadrant);

/** Finalize the construction of the new forest after adding quadrants.
 * \param [in,out] build    The building context will be deallocated inside.
 * \return                  A valid forest object.
 *                          Its revision number is set to zero.
 */
p4est_t            *p4est_build_complete (p4est_build_t * build);

SC_EXTERN_C_END;

#endif /* ! P4EST_BUILD_H */
