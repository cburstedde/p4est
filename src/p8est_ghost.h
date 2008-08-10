/*
  This file is part of p4est.
  p4est is a C library to manage a parallel collection of quadtrees and/or
  octrees.

  Copyright (C) 2008 Carsten Burstedde, Lucas Wilcox.

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

#ifndef P8EST_GHOST_H
#define P8EST_GHOST_H

#include <p8est.h>

/** Gets the processor id of a quadrant's owner.
 * The quadrant can lie outside of a tree across faces (and only faces).
 *
 * \param [in] p8est  The forest in which to search for a quadrant.
 * \param [in] treeid The tree id for which the quadrant belongs.
 * \param [in] face   Supply a face direction if known, or -1 otherwise.
 * \param [in] q      The quadrant that is being searched for.
 *
 * \return Processor id of the owner
 *                or -1 if the quadrant lies outside of the mesh.
 *
 * \warning Does not work for tree edge or corner neighbors.
 */
int                 p8est_quadrant_find_owner (p8est_t * p8est,
                                               p4est_topidx_t treeid,
                                               int face,
                                               const p8est_quadrant_t * q);

/** Builds the ghost layer.
 *
 * This will gather the quadrants from each neighboring proc to build one layer
 * of face, edge and corner based ghost elements around the ones they own.
 *
 * \param [in] p8est              The forest for which the ghost layer will
 *                                be generated.
 * \param [in] include_diagonals  Include neighbors across edges and corners.
 * \param [in,out] ghost_layer    An array of quadrants which make up the
 *                                ghost layer around \a p4est.  Their piggy1
 *                                data member is filled with their owner's
 *                                tree and processor ids.  Quadrants will be
 *                                ordered in \c p8est_quadrant_compare_piggy
 *                                order.  These will be quadrants inside the
 *                                neighboring tree i.e., \c
 *                                p4est_quadrant_is_inside is true for the
 *                                quadrant and the neighboring tree.
 * \param [out] ghost_owner       If not NULL, this array is allocated and
 *                                filled with one rank per ghost.
 * \return                        Returns false if it fails due to violated
 *                                2:1 constraints, true otherwise.
 */
bool                p8est_build_ghost_layer (p8est_t * p8est,
                                             bool include_diagonals,
                                             sc_array_t * ghost_layer,
                                             int **ghost_owner);

/** Checks if quadrant exists in the local forest or the ghost layer.
 *
 * For quadrants across tree corners/edges it checks if the quadrant exists
 * in any of the corner/edge neighbors.
 *
 * \param [in]  p4est        The forest in which to search for \a q
 * \param [in]  ghost_layer  The ghost layer in which to search for \a q
 * \param [in]  treeid       The tree id for which \a q belongs.
 * \param [in]  q            The quadrant that is being searched for.
 * \param [out] exists_arr   Filled for tree corner/edge cases.  An entry
 *                           for each corner/edge neighbor is set to 1 if
 *                           it exists in the local forest or ghost_layer.
 *
 * \return true if the quadrant exists in the local forest or in the
 *              ghost_layer, and false if doesn't exist in either.
 */
bool                p8est_quadrant_exists (p8est_t * p8est,
                                           sc_array_t * ghost_layer,
                                           p4est_topidx_t treeid,
                                           const p8est_quadrant_t * q,
                                           sc_array_t * exists_arr);

/** Checks a p8est to see if it is balanced.
 *
 * This function builds the ghost layer and discards it when done.
 *
 * \param [in] p8est  The p8est to be tested.
 * \return Returns true if balanced, false otherwise.
 */
bool                p8est_is_balanced (p8est_t * p8est);

#endif /* !P8EST_GHOST_H */
