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

/** \file p8est_geometry.h transforms from vertex frame to physical space.
 *
 * \ingroup p8est
 */

#ifndef P8EST_GEOMETRY_H
#define P8EST_GEOMETRY_H

#include <p8est_connectivity.h>

SC_EXTERN_C_BEGIN;

/** This object encapsulates a custom geometry transformation. */
typedef struct p8est_geometry p8est_geometry_t;

/** Forward transformation from the reference unit square to physical space.
 * The physical space "xyz" is user-defined, currently used for VTK output.
 */
typedef void        (*p8est_geometry_X_t) (p8est_geometry_t * geom,
                                           p4est_topidx_t which_tree,
                                           const double abc[3],
                                           double xyz[3]);

/** Destructor prototype for a user-allocated \a p8est_geometry_t.
 * It is invoked by p8est_geometry_destroy.  If the user chooses to
 * reserve the structure statically, simply don't call p4est_geometry_destroy.
 */
typedef void        (*p8est_geometry_destroy_t) (p8est_geometry_t * geom);

/** This structure can be created by the user,
 * p4est will never change its contents.
 */
struct p8est_geometry
{
  const char         *name;     /**< User's choice is arbitrary. */
  void               *user;     /**< User's choice is arbitrary. */
  p8est_geometry_X_t  X;        /**< Coordinate transformation. */
  p8est_geometry_destroy_t destroy;     /**< Destructor called by
                                             p8est_geometry_destroy.  If
                                             NULL, P4EST_FREE is called. */
};

/** Can be used to conveniently destroy a geometry structure.
 * The user is free not to call this function at all if they handle the
 * memory of the p8est_geometry_t in their own way.
 */
void                p8est_geometry_destroy (p8est_geometry_t * geom);

/** Create a geometry structure based on the vertices in a connectivity.
 * The transformation is constructed using trilinear interpolation.
 * \param [in] conn A p8est_connectivity_t with valid vertices.  We do NOT
 *                  take ownership and expect this structure to stay alive.
 * \return          Geometry structure; use with p4est_geometry_destroy.
 */
p8est_geometry_t   *p8est_geometry_new_connectivity (p8est_connectivity_t *
                                                     conn);

/** Create a geometry structure for the spherical shell of 24 trees.
 * \param [in] conn Result of p8est_connectivity_new_shell or equivalent.
 *                  We do NOT take ownership and expect it to stay alive.
 * \param [in] R2   The outer radius of the shell.
 * \param [in] R1   The inner radius of the shell.
 * \return          Geometry structure; use with p4est_geometry_destroy.
 */
p8est_geometry_t   *p8est_geometry_new_shell (p8est_connectivity_t * conn,
                                              double R2, double R1);

/** Create a geometry structure for the solid sphere of 13 trees.
 * \param [in] conn Result of p8est_connectivity_new_sphere or equivalent.
 *                  We do NOT take ownership and expect it to stay alive.
 * \param [in] R2   The outer radius of the sphere.
 * \param [in] R1   The outer radius of the inner shell.
 * \param [in] R0   The inner radius of the inner shell.
 * \return          Geometry structure; use with p4est_geometry_destroy.
 */
p8est_geometry_t   *p8est_geometry_new_sphere (p8est_connectivity_t * conn,
                                               double R2, double R1,
                                               double R0);

/** Create a geometry structure for the torus.
 *
 * This geometry maps a revolution torus, obtained using
 * \ref p8est_connectivity_new_torus
 *
 * The torus is divided into into segments arround the revolution axis,
 * each segments is made of 5 trees; so here we provided the geometric
 * transformation in a piecewise manner for each tree of the connectivity.
 *
 * \param [in] conn Result of p8est_connectivity_new_torus or equivalent.
 *                  We do NOT take ownership and expect it to stay alive.
 *
 * \param [in] R0   The inner radius of the 2d disk slice.
 * \param [in] R1   The outer radius of the 2d disk slice.
 * \param [in] R2   The outer radius of the torus.
 * \return          Geometry structure; use with p4est_geometry_destroy.
 */
p8est_geometry_t   *p8est_geometry_new_torus (p8est_connectivity_t * conn,
                                              double R0, double R1,
                                              double R2);

SC_EXTERN_C_END;

#endif /* !P8EST_GEOMETRY_H */
