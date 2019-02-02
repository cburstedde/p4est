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

#ifndef P4EST_TO_P8EST_SPHERES_H
#define P4EST_TO_P8EST_SPHERES_H

#ifdef P4EST_H
#error "The include files p4est.h and p4est_to_p8est_spheres.h cannot be combined"
#endif

#include <p4est_to_p8est.h>

#define P4EST_CUBE_DIAG                 P8EST_CUBE_DIAG

#define p4est_sphere_t                  p8est_sphere_t
#define p4est_spheres_t                 p8est_spheres_t

#define p4est_quadrant_sphere_box       p8est_quadrant_sphere_box
#define p4est_sphere_match_approx       p8est_sphere_match_approx
#define p4est_sphere_match_exact        p8est_sphere_match_exact

#endif /* !P4EST_TO_P8EST_SPHERES_H */
