/*
  This file is part of p4est.
  p4est is a C library to manage a parallel collection of quadtrees and/or
  octrees.

  Copyright (C) 2007-2009 Carsten Burstedde, Lucas Wilcox.

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

#ifndef P8EST_GEOMETRY_H
#define P8EST_GEOMETRY_H

#include <p8est.h>

/* Naming convention for different coordinate systems:
 * xyz              physical domain
 * abc              p4est vertex domain
 * rst              reference domain [-1,1]^3
 */

typedef struct p8est_geometry p8est_geometry_t;

typedef void        (*p8est_geometry_X_t) (p8est_geometry_t * geom,
                                           p4est_topidx_t which_tree,
                                           const double abc[3],
                                           double xyz[3]);
typedef double      (*p8est_geometry_D_t) (p8est_geometry_t * geom,
                                           p4est_topidx_t which_tree,
                                           const double abc[3]);
typedef double      (*p8est_geometry_J_t) (p8est_geometry_t * geom,
                                           p4est_topidx_t which_tree,
                                           const double abc[3],
                                           double J[3][3]);

struct p8est_geometry
{
  p8est_geometry_X_t  X;
  p8est_geometry_D_t  D;
  p8est_geometry_J_t  J, Jit;   /* both return the determinant of J */
};

/** Compute the inverse transpose Jacobian by calling
 * the geom->J function and transpose inverting the result.
 * \return          The determinant of the Jacobian J (not of Jit).
 */
double              p8est_geometry_Jit (p8est_geometry_t * geom,
                                        p4est_topidx_t which_tree,
                                        const double abc[3],
                                        double Jit[3][3]);

/** The identity transformation.
 */
void                p8est_geometry_identity_X (p8est_geometry_t * geom,
                                               p4est_topidx_t which_tree,
                                               const double abc[3],
                                               double xyz[3]);

/** The Jacobi determinant of the identity transformation.
 * \return          The determinant of the Jacobian.
 */
double              p8est_geometry_identity_D (p8est_geometry_t * geom,
                                               p4est_topidx_t which_tree,
                                               const double abc[3]);

/** The Jacobian matrix of the identity transformation.
 * \return          The determinant of the Jacobian J.
 */
double              p8est_geometry_identity_J (p8est_geometry_t * geom,
                                               p4est_topidx_t which_tree,
                                               const double abc[3],
                                               double J[3][3]);

/** Create a geometry structure for the identity transformation.
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

#endif /* !P8EST_GEOMETRY_H */
