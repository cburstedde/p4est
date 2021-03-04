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

/** \file p4est_geometry.h transforms from vertex frame to physical space.
 *
 * \ingroup p4est
 */

#ifndef P4EST_GEOMETRY_H
#define P4EST_GEOMETRY_H

#include <p4est_connectivity.h>

SC_EXTERN_C_BEGIN;

/** This object encapsulates a custom geometry transformation. */
typedef struct p4est_geometry p4est_geometry_t;

/** Forward transformation from the reference unit square to physical space.
 * Note that the two-dimensional connectivities have 3D vertex coordinates
 * that can be used in the transformation if so desired.
 * The physical space "xyz" is user-defined, currently used for VTK output.
 */
typedef void        (*p4est_geometry_X_t) (p4est_geometry_t * geom,
                                           p4est_topidx_t which_tree,
                                           const double abc[3],
                                           double xyz[3]);

/** Destructor prototype for a user-allocated \a p4est_geometry_t.
 * It is invoked by p4est_geometry_destroy.  If the user chooses to
 * reserve the structure statically, simply don't call p4est_geometry_destroy.
 */
typedef void        (*p4est_geometry_destroy_t) (p4est_geometry_t * geom);

/** This structure can be filled or allocated by the user.
 * p4est will never change its contents.
 */
struct p4est_geometry
{
  const char         *name;     /**< User's choice is arbitrary. */
  void               *user;     /**< User's choice is arbitrary. */
  p4est_geometry_X_t  X;        /**< Coordinate transformation. */
  p4est_geometry_destroy_t destroy;     /**< Destructor called by
                                             p4est_geometry_destroy.  If
                                             NULL, P4EST_FREE is called. */
};

/** Can be used to conveniently destroy a geometry structure.
 * The user is free not to call this function at all if they handle the
 * memory of the p4est_geometry_t in their own way.
 */
void                p4est_geometry_destroy (p4est_geometry_t * geom);

/** Create a geometry structure based on the vertices in a connectivity.
 * The transformation is constructed using bilinear interpolation.
 * \param [in] conn A p4est_connectivity_t with valid vertices.  We do NOT
 *                  take ownership and expect this structure to stay alive.
 * \return          Geometry structure; use with p4est_geometry_destroy.
 */
p4est_geometry_t   *p4est_geometry_new_connectivity (p4est_connectivity_t *
                                                     conn);

/** Create a geometry for mapping the 3d sphere using 2d connectivity icosahedron.
 *
 */
p4est_geometry_t   *p4est_geometry_new_icosahedron (p4est_connectivity_t *
                                                    conn, double R);

/** Create a geometry for mapping 2d shell.
 *  This a direct adaptation of geometric shell in 3d.
 */
p4est_geometry_t   *p4est_geometry_new_shell2d (p4est_connectivity_t * conn,
                                                double R2, double R1);

/**
 * disk2d geometry associated to disk2d connectivity.
 *
 * \param[in] R0 radius of the inner circle
 * \param[in] R1 radius of the outer circle (external border)
 *
 *
 * This geometry is meant to be used with the disk2d connectivity,
 * \ref p4est_connectivity_new_disk2d which is a 5-tree connectivty
 * to map the disk.
 */
p4est_geometry_t   *p4est_geometry_new_disk2d (p4est_connectivity_t * conn,
                                               double R0, double R1);

SC_EXTERN_C_END;

#endif /* !P4EST_GEOMETRY_H */
