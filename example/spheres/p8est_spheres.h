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

/** \file p8est_spheres.h
 * Functions to create a random distribution of spheres in the unit cube.
 * These spheres can be intersected with quadrants.
 * This file ignores vertex coordinates in a \ref p8est_connectivity_t,
 * and neither does transformations using \ref p8est_geometry_t.
 */

#ifndef P8EST_SPHERES_H
#define P8EST_SPHERES_H

#include <p8est.h>

SC_EXTERN_C_BEGIN;

/** The cube diagonal in 3D is sqrt (3.). */
#define P8EST_CUBE_DIAG                 (1.7320508075688772)

typedef struct p8est_sphere
{
  double              center[P8EST_DIM];
  double              radius;
}
p8est_sphere_t;

typedef struct p8est_spheres
{
  p8est_t            *p4est;
}
p8est_spheres_t;

/** Use quadrant coordinates to represent it with a sphere type.
 * \param [in] quadrant     Valid quadrant.
 * \param [out] sph         On output, contains quadrant dimensions.
 */
void                p8est_quadrant_sphere_box
  (const p8est_quadrant_t * quadrant, p8est_sphere_t * sph);

/** Optimistically check for intersection of sphere surface and cube volume.
 * \param [in] box  We use the sphere data structure to define a box
 *                  by its center and half width, the latter stored as radius.
 * \param [in] sph  Sphere whose surface is checked against the box volume.
 * \param [in] t    Factor wrt. the sphere's radius for thickness of surface.
 * \return          True if there might be an intersection, false if not.
 */
int                 p8est_sphere_match_approx
  (const p8est_sphere_t * box, const p8est_sphere_t * sph, double t);

/** Exactly check for intersection of sphere surface and cube volume.
 * \param [in] box  We use the sphere data structure to define a box
 *                  by its center and half width, the latter stored as radius.
 * \param [in] sph  Sphere whose surface is checked against the box volume.
 * \param [in] t    Factor wrt. the sphere's radius for thickness of surface.
 * \return          True if and only if sphere surface intersects box volume.
 */
int                 p8est_sphere_match_exact
  (const p8est_sphere_t * box, const p8est_sphere_t * sph, double t);

SC_EXTERN_C_END;

#endif /* !P8EST_SPHERES_H */
