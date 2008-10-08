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

#ifndef P4EST_GHOST_H
#define P4EST_GHOST_H

#include <p4est.h>

SC_EXTERN_C_BEGIN;

/** Gets the processor id of a quadrant's owner.
 * The quadrant can lie outside of a tree across faces (and only faces).
 *
 * \param [in] p4est  The forest in which to search for a quadrant.
 * \param [in] treeid The tree id for which the quadrant belongs.
 * \param [in] face   Supply a face direction if known, or -1 otherwise.
 * \param [in] q      The quadrant that is being searched for.
 *
 * \return Processor id of the owner
 *                or -1 if the quadrant lies outside of the mesh.
 *
 * \warning Does not work for tree edge or corner neighbors.
 */
int                 p4est_quadrant_find_owner (p4est_t * p4est,
                                               p4est_topidx_t treeid,
                                               int face,
                                               const p4est_quadrant_t * q);

/** Builds the ghost layer.
 *
 * This will gather the quadrants from each neighboring proc to build
 * one layer of face and corner based ghost elements around the ones they own.
 *
 * \param [in] p4est              The forest for which the ghost layer will
 *                                be generated.
 * \param [in] include_diagonals  Include neighbors across corners.
 * \param [in,out] ghost_layer    An array of quadrants which make up the
 *                                ghost layer around \a p4est.  Their piggy1
 *                                data member is filled with their owner's
 *                                tree and processor ids.  Quadrants will be
 *                                ordered in \c p4est_quadrant_compare_piggy
 *                                order.  These will be quadrants inside the
 *                                neighboring tree i.e., \c
 *                                p4est_quadrant_is_inside is true for the
 *                                quadrant and the neighboring tree.
 * \param [out] ghost_owner       If not NULL, this array is allocated and
 *                                filled with one rank per ghost.
 * \return                        Returns false if it fails due to violated
 *                                2:1 constraints, true otherwise.
 */
bool                p4est_build_ghost_layer (p4est_t * p4est,
                                             bool include_diagonals,
                                             sc_array_t * ghost_layer,
                                             int **ghost_owner);

/** Checks if quadrant exists in the local forest or the ghost layer.
 *
 * For quadrants across tree boundaries it checks if the quadrant exists
 * across any face, but not across corners.
 *
 * \param [in]  p4est        The forest in which to search for \a q.
 * \param [in]  ghost_layer  The ghost layer in which to search for \a q.
 * \param [in]  treeid       The tree id for which \a q belongs.
 * \param [in]  q            The quadrant that is being searched for.
 * \param [in,out] face      On input, face id across which \a q was created.
 *                           On output, the neighbor's face number.
 * \param [out] hang         If not NULL, signals that q is bigger than
 *                           the quadrant it came from.  The child id
 *                           of that originating quadrant is passed in hang.
 *                           On output, hang holds the hanging face number
 *                           of \a q that is in contact with its originator.
 * \param [out] owner_rank   Filled with the rank of the owner
 *                           if it is found and undefined otherwise.
 *
 * \return      Returns the local number of \a q if the quadrant exists
 *              in the local forest or in the ghost_layer.  Otherwise,
 *              returns -2 for a domain boundary and -1 if not found.
 */
p4est_locidx_t      p4est_face_quadrant_exists (p4est_t * p4est,
                                                sc_array_t * ghost_layer,
                                                p4est_topidx_t treeid,
                                                const p4est_quadrant_t * q,
                                                int *face, int *hang,
                                                int *owner_rank);

/** Checks if quadrant exists in the local forest or the ghost layer.
 *
 * For quadrants across tree corners it checks if the quadrant exists
 * in any of the corner neighbors.
 *
 * \param [in]  p4est        The forest in which to search for \a q
 * \param [in]  ghost_layer  The ghost layer in which to search for \a q
 * \param [in]  treeid       The tree id for which \a q belongs.
 * \param [in]  q            The quadrant that is being searched for.
 * \param [out] exists_arr   Filled for tree corner cases.  An entry
 *                           for each corner neighbor is set to 1 if
 *                           it exists in the local forest or ghost_layer.
 *
 * \return true if the quadrant exists in the local forest or in the
 *              ghost_layer, and false if doesn't exist in either.
 */
bool                p4est_quadrant_exists (p4est_t * p4est,
                                           sc_array_t * ghost_layer,
                                           p4est_topidx_t treeid,
                                           const p4est_quadrant_t * q,
                                           sc_array_t * exists_arr);

/** Checks a p4est to see if it is balanced.
 *
 * This function builds the ghost layer and discards it when done.
 *
 * \param [in] p4est  The p4est to be tested.
 * \return Returns true if balanced, false otherwise.
 */
bool                p4est_is_balanced (p4est_t * p4est);

SC_EXTERN_C_END;

#endif /* !P4EST_GHOST_H */
