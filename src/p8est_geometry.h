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

#ifndef P8EST_GEOMETRY_H
#define P8EST_GEOMETRY_H

#include <p8est.h>

SC_EXTERN_C_BEGIN;

typedef struct p8est_geometry p8est_geometry_t;

/** Forward transformation from vertex frame to physical space.
 * The vertex space "abc" is defined per octree and spanned by the vertices
 * at its corners; see p8est_connectivity.h.
 * The physical space "xyz" is user-defined, currently used for VTK output.
 */
typedef void        (*p8est_geometry_X_t) (p8est_geometry_t * geom,
                                           p4est_topidx_t which_tree,
                                           const double abc[3],
                                           double xyz[3]);

/** This structure can be created by the user,
 * p4est will never change its contents.
 */
struct p8est_geometry
{
  const char         *name;     /**< User's choice is arbitrary. */
  void               *user;     /**< User's choice is arbitrary. */
  p8est_geometry_X_t  X;        /**< Coordinate transformation. */
};

/** Create a geometry structure for the identity transformation.
 * This function is just for demonstration since a NULL geometry works too.
 * \return          Geometry structure which must be freed with P4EST_FREE.
 */
p8est_geometry_t   *p8est_geometry_new_identity (void);

/** Create a geometry structure for the spherical shell of 24 trees.
 * This is suitable for forests obtained with p8est_connectivity_new_shell.
 * \param [in] R2   The outer radius of the shell.
 * \param [in] R1   The inner radius of the shell.
 * \return          Geometry structure which must be freed with P4EST_FREE.
 */
p8est_geometry_t   *p8est_geometry_new_shell (double R2, double R1);

/** Create a geometry structure for the solid sphere of 13 trees.
 * This is suitable for forests obtained with p8est_connectivity_new_sphere.
 * \param [in] R2   The outer radius of the sphere.
 * \param [in] R1   The outer radius of the inner shell.
 * \param [in] R0   The inner radius of the inner shell.
 * \return          Geometry structure which must be freed with P4EST_FREE.
 */
p8est_geometry_t   *p8est_geometry_new_sphere (double R2, double R1,
                                               double R0);

SC_EXTERN_C_END;

#endif /* !P8EST_GEOMETRY_H */
