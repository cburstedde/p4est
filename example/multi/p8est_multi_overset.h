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

/** Execute multi-mesh overset algorithm.
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
 * \param [in] bgp4est          For \a myrole zero, the background forest.
 *                              NULL otherwise.
 * \param [in] qpoints          Query points: 4-double-tuples (x, y, z, v).
 *                              \b a is -1 for the overset boundary, -2 for
 *                              the wall boundary, and a non-negative
 *                              representative volume of the point otherwise.
 *                              The points are gathered process-local over
 *                              all overset meshes present on this process.
 *                              NULL for \a myrole zero.
 */
void                 p8est_multi_overset
  (sc_MPI_Comm glocomm, sc_MPI_Comm headcomm, sc_MPI_Comm rolecomm,
   int myrole, int num_meshes, const int *mesh_offsets,
   p4est_t *bgp4est, sc_array_t *qpoints);

SC_EXTERN_C_END;

#endif /* !P8EST_MULTI_OVERSET_H */
