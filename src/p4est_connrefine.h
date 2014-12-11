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

/** \file p4est_connrefine.h
 *
 * Refine a connectivity.
 *
 * \ingroup p4est
 */

#ifndef P4EST_CONNREFINE_H
#define P4EST_CONNREFINE_H

#include <p4est_connectivity.h>

/** Uniformly refine a connectivity.
 * This is useful if you would like to uniformly refine by something other
 * than a power of 2.
 *
 * \param [in] conn         a valid connectivity
 * \param [in] num_per_edge the number of new trees in each direction
 *
 * \return a refined connectivity.
 */
p4est_connectivity_t *p4est_connectivity_refine (p4est_connectivity_t * conn,
                                                 int num_per_edge);
#endif
