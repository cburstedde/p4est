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

#ifndef P4EST_SEARCH_BUILD_H
#define P4EST_SEARCH_BUILD_H

#include <p4est_search.h>

/** \file p4est_search_build.h
 * Create a new p4est object by running \ref p4est_search_local.
 * This allows to create a heavily coarsened forest in one pass.
 */

/** Context object for building a new p4est from search callbacks.
 */
typedef struct p4est_search_build p4est_search_build_t;

/** Allocate a context for building a new forest.
 * \param [in] from         This forest is used as a template for creation.
 * \param [in] data_size    Data size of the created forest, may be zero.
 * \param [in] init_fn      This functions is called for created quadrants.
 *                          NULL leaves the quadrant data uninitialized.
 *                          This function can be overridden on a per-quadrant
 *                          basis in \ref p4est_search_build_local.
 * \param [in] user_pointer Registered into the newly built forest.
 * \return                  A context that needs to be processed further.
 */
p4est_search_build_t *p4est_search_build_new (p4est_t * from,
                                              size_t data_size,
                                              p4est_init_t init_fn,
                                              void *user_pointer);

/** This function is usable from a \ref p4est_search_local_t callback.
 * \param [in,out] build    The building context must be passed through.
 * \param [in] which_tree   The tree number is passed from the search callback.
 * \param [in] quadrant     The quadrant is passed from the search callback.
 * \param [in] local_num    The quadrant number is passed from search callback.
 * \param [in] init_quadrant    If NULL, use value from \ref
 *                          p4est_search_build_new.  Otherwise, this function
 *                          is passed on and may be used to initialize this
 *                          quadrant's data.
 *
 * \return                  TODO: figure out what to do on return.
 */
int                 p4est_search_build_add (p4est_search_build_t * build,
                                            p4est_topidx_t which_tree,
                                            p4est_quadrant_t * quadrant,
                                            p4est_locidx_t local_num,
                                            p4est_init_t init_quadrant);

/** Finalize the construction of the new forest after the search.
 * \param [in,out] build    The building context will be deallocated inside.
 * \return                  A valid forest object.
 */
p4est_t            *p4est_search_build_complete (p4est_search_build_t *
                                                 build);

#endif /* ! P4EST_SEARCH_BUILD_H */
