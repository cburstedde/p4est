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

/** \file p4est_dune.h
 *
 * Provide various functions to interface dune to a p4est mesh.
 *
 * For example, populate element corner and face number tables, where the
 * range for each codimension begins with 0 and is contiguous.
 *
 * We also provide a face-only variant of the iterator that visits all small
 * quadrants on one side of a hanging face from separate callbacks in order.
 *
 * \ingroup p4est
 */

#ifndef P4EST_DUNE_H
#define P4EST_DUNE_H

#include <p4est_iterate.h>

SC_EXTERN_C_BEGIN;

/** Volume and face iterator over a 2:1 balanced forest.
 * See \ref p4est_dune_iterate for a detailed reference.
 */
void                p4est_dune_iterate_balanced (p4est_t *p4est,
                                                 p4est_ghost_t *ghost_layer,
                                                 void *user_data,
                                                 p4est_iter_volume_t
                                                 iter_volume,
                                                 p4est_iter_face_t iter_face);

/** Execute user supplied callbacks at every local volume and face.
 *
 * Execute the user-supplied callback functions at every volume and face in
 * the local forest.  When multiple small quadrants are on the other side of
 * the face of a large quadrant, the callback is executed for every small
 * face neighbor in turn, with the same large quadrant on the other side.
 * The \a user_data pointer is not touched by the iterator, only passed to
 * each of the callbacks.  Any of the callbacks may be NULL.
 * The callback functions are interspersed with each other, i.e. some face
 * callbacks will occur between volume callbacks:
 *
 * 1. Volume callbacks occur in the sorted Morton-index order.
 * 2. A face callback is not executed until after the volume callbacks have
 *    been executed for the quadrants that share it.
 * 3. The face callbacks for the small elements at the same hanging face
 *    are executed right after each other in their Morton-index order.
 * 4. Callbacks are not executed at faces that only involve ghost
 *    quadrants, i.e. that are not adjacent in the local section of the
 *    forest.
 *
 * The convention for the context information passed to the callbacks
 * is the same as for \ref p4est_iterate with the following exception:
 * For a hanging face, we call the callback anew for every small neighbor
 * quadrant with the same large neighbor quadrant.  The hanging side data
 * (see \ref p4est_iter_face_side_hanging_t) is populated at index [0] each
 * time, and we store the number of the small quadrant within the order of
 * its face siblings in the variable is.hanging.quadid[1], values in [0, 2).
 *
 * \param[in] p4est          The forest to iterate over.
 *                           It is not required to be 2:1 balanced.
 * \param[in] ghost_layer    Valid ghost structure or NULL.
 * \param[in,out] user_data  optional context to supply to each callback.
 * \param[in] iter_volume    callback function for every quadrant interior.
 * \param[in] iter_face      callback function for every face between.
 */
void                p4est_dune_iterate (p4est_t *p4est,
                                        p4est_ghost_t *ghost_layer,
                                        void *user_data,
                                        p4est_iter_volume_t iter_volume,
                                        p4est_iter_face_t iter_face);

SC_EXTERN_C_END;

#endif /* P4EST_DUNE_H */
