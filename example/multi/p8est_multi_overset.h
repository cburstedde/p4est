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

/** Execute multi-mesh overset algorithm.
 * \param [in] glorank          Rank within global communicator.
 * \param [in] myrole           Index of mesh: 0 for background mesh,
 *                              starting from 1 for overset meshes.
 * \param [in] num_meshes       Number of meshes including background.
 * \param [in] points           An allocated sc_array of query points
 *                              of user-defined type. The data is never
 *                              touched by p4est but only passed to
 *                              user-defined callbacks.
 * \param [in] mesh_offsets     Array of ascending global ranks,
 *                              one for the first of each mesh, and then
 *                              one more for the end (exclusive of the last).
 */
void                p8est_multi_overset
  (int glorank, int myrole, int num_meshes, sc_array_t * points,
   const int *mesh_offsets);

#endif /* !P8EST_MULTI_OVERSET_H */
