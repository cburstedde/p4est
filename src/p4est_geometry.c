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

/**
 * \file p4est_geometry.c
 * We provide the identity transformation for reference.
 * Please implement p4est_geometry_t as you see fit.
 */

#include <p4est_geometry.h>

static void
p4est_geometry_identity_X (p4est_geometry_t * geom,
                           p4est_topidx_t which_tree,
                           const double abc[3], double xyz[3])
{
  memcpy (xyz, abc, 3 * sizeof (double));
}

p4est_geometry_t   *
p4est_geometry_new_identity (void)
{
  p4est_geometry_t   *geom;

  geom = P4EST_ALLOC_ZERO (p4est_geometry_t, 1);

  geom->name = "p4est_identity";
  geom->X = p4est_geometry_identity_X;

  return geom;
}
