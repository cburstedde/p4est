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

#ifndef P8EST_BALANCE_H
#define P8EST_BALANCE_H

#include <p8est.h>

SC_EXTERN_C_BEGIN;

/** Determines if quadrant \a q causes quadrant \a p to split under the given
 * \a balance condition.
 *
 * \param [in] q         Test quadrant.
 * \param [in] p         Trial quadrant.
 * \param [in] balance   Balance condition.
 * \param [out] seeds    optional array: if \a seeds is not NULL, then it will
 *                       be resized and filled with descendants of \a p such
 *                       that the coarsest balanced subtree rooted at \a p
 *                       that contains all of \a seeds is also the coarset
 *                       subtree rooted at \a p that is entirely balanced with
 *                       \a q.
 * \return               True if \a q causes \a p to split.
 */
int                 p8est_balance_seeds (p8est_quadrant_t * q,
                                         p8est_quadrant_t * p,
                                         p8est_connect_type_t balance,
                                         sc_array_t * seeds);

/** Same as p8est_balance_seeds, optimized for the case when it is already
 * known that \a q is outside of a certain \a face of \a p.
 */
int                 p8est_balance_seeds_face (p8est_quadrant_t * q,
                                              p8est_quadrant_t * p,
                                              int face, p8est_connect_type_t
                                              balance, sc_array_t * seeds);

/** Same as p8est_balance_seeds, optimized for the case when it is already
 * known that \a q is outside of a certain \a edge of \a p.
 */
int                 p8est_balance_seeds_edge (p8est_quadrant_t * q,
                                              p8est_quadrant_t * p,
                                              int face, p8est_connect_type_t
                                              balance, sc_array_t * seeds);

/** Same as p8est_balance_seeds, optimized for the case when it is already
 * known that \a q is outside of a certain \a corner of \a p.
 */
int                 p8est_balance_seeds_corner (p8est_quadrant_t * q,
                                                p8est_quadrant_t * p,
                                                int face, p8est_connect_type_t
                                                balance, sc_array_t * seeds);
SC_EXTERN_C_END;

#endif
