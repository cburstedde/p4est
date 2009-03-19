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

#include <p8est_geometry.h>

double
p8est_geometry_Jit (p8est_geometry_t * geom,
                    p4est_topidx_t which_tree,
                    const double xyz[3], double Jit[3][3])
{
  double              J[3][3];
  double              detJ, idetJ;

  idetJ = 1. / (detJ = geom->J (geom, which_tree, xyz, J));

  Jit[0][0] = (J[1][1] * J[2][2] - J[1][2] * J[2][1]) * idetJ;
  Jit[0][1] = (J[1][2] * J[2][0] - J[1][0] * J[2][2]) * idetJ;
  Jit[0][2] = (J[1][0] * J[2][1] - J[1][1] * J[2][0]) * idetJ;

  Jit[1][0] = (J[0][2] * J[2][1] - J[0][1] * J[2][2]) * idetJ;
  Jit[1][1] = (J[0][0] * J[2][2] - J[0][2] * J[2][0]) * idetJ;
  Jit[1][2] = (J[0][1] * J[2][0] - J[0][0] * J[2][1]) * idetJ;

  Jit[2][0] = (J[0][1] * J[1][2] - J[1][1] * J[0][2]) * idetJ;
  Jit[2][1] = (J[0][2] * J[1][0] - J[1][2] * J[0][0]) * idetJ;
  Jit[2][2] = (J[0][0] * J[1][1] - J[1][0] * J[0][1]) * idetJ;

  return detJ;
}

void
p8est_geometry_identity_X (p8est_geometry_t * geom,
                           p4est_topidx_t which_tree,
                           const double xyz[3], double XYZ[3])
{
  memcpy (XYZ, xyz, 3 * sizeof (double));
}

double
p8est_geometry_identity_D (p8est_geometry_t * geom,
                           p4est_topidx_t which_tree, const double xyz[3])
{
  return 1.;
}

double
p8est_geometry_identity_J (p8est_geometry_t * geom,
                           p4est_topidx_t which_tree,
                           const double xyz[3], double J[3][3])
{
  J[0][0] = J[1][1] = J[2][2] = 1.;
  J[0][1] = J[1][2] = J[2][0] = 0.;
  J[1][0] = J[2][1] = J[0][2] = 0.;

  return 1.;
}

p8est_geometry_t   *
p8est_geometry_new_identity (void)
{
  p8est_geometry_t   *geom;

  geom = P4EST_ALLOC (p8est_geometry_t, 1);

  geom->X = p8est_geometry_identity_X;
  geom->D = p8est_geometry_identity_D;
  geom->J = geom->Jit = p8est_geometry_identity_J;      /* identical here */

  return geom;
}
