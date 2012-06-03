/*
  This file is part of p4est.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2010 The University of Texas System
  Copyright (C) 2012 Carsten Burstedde
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

#ifndef P8EST_IO_H
#define P8EST_IO_H

#include <p8est.h>

/** Extract all quadrant x y z level data on the local processor.
 * Optionally extracts the quadrant data as well into a separate array.
 * \param [in] p8est    The forest is not modified.
 * \param [in,out] data If not NULL, pointer to a pointer that will be set
 *                      to a newly allocated array with per-quadrant data.
 *                      Must be NULL if p4est->data_size == 0.
 * \return              An array of type p8est_qcoord_t that contains
 *                      x y z level for each quadrant on this processor.
 *                      The tree information is not extracted.
 */
sc_array_t        *p4est_deflate_quadrants (p8est_t * p8est,
                                              sc_array_t **data);

#endif /* !P8EST_IO_H */
