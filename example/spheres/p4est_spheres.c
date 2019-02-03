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

#ifndef P4_TO_P8
#include <p4est_bits.h>
#include "p4est_spheres.h"
#else
#include <p8est_bits.h>
#include "p8est_spheres.h"
#endif

#define P4EST_CUBE_DIAG_INV     (1. / P4EST_CUBE_DIAG)
#define P4EST_ROOT_LEN_INV      (1. / P4EST_ROOT_LEN)

void
p4est_quadrant_sphere_box (const p4est_quadrant_t * quadrant,
                           p4est_sphere_t * sph)
{
  double              r;

  P4EST_ASSERT (p4est_quadrant_is_valid (quadrant));
  P4EST_ASSERT (sph != NULL);

  sph->radius = r =
    .5 * P4EST_QUADRANT_LEN (quadrant->level) * P4EST_ROOT_LEN_INV;
  sph->center[0] = P4EST_ROOT_LEN_INV * quadrant->x + r;
  sph->center[1] = P4EST_ROOT_LEN_INV * quadrant->y + r;
#ifdef P4_TO_P8
  sph->center[2] = P4EST_ROOT_LEN_INV * quadrant->z + r;
#endif
}

/* Optimistic check.  We may return 1 without being sure. */
int
p4est_sphere_match_approx (const p4est_sphere_t * box,
                           const p4est_sphere_t * sph, double t)
{
  const double        br = box->radius;
  const double        sr = sph->radius;
  double              fdist[P4EST_DIM];
  double              fmax;
  double              cd, dd;

  P4EST_ASSERT (t >= 0.);

  fdist[0] = fabs (box->center[0] - sph->center[0]);
  fdist[1] = fabs (box->center[1] - sph->center[1]);
  fmax = SC_MAX (fdist[0], fdist[1]);
#ifdef P4_TO_P8
  fdist[2] = fabs (box->center[2] - sph->center[2]);
  fmax = SC_MAX (fmax, fdist[2]);
#endif

  /* sphere and quadrant are necessarily disjoint */
  cd = br + (1. + t) * sr;
  if (fmax > cd) {
    /* this result is correct */
    return 0;
  }

  /* quadrant contained in sphere's in-cube: no surface intersection */
  dd = P4EST_CUBE_DIAG_INV * (1. - t) * sr - br;
  if (fmax < dd) {
    /* if dd is negative, then this case does not arise */
    return 0;
  }

  /* we are returning 1 below, no need for this test */
#if 0
  /* sphere contained in quadrant is a definite match */
  dd = br - (1. + t) * sr;
  if (fmax < dd) {
    /* if dd is negative, then this case does not arise */
    return 1;
  }
#endif

  /* otherwise, this match may be true, but we are not certain */
  return 1;
}

/** Exactly check for intersection of sphere surface and cube volume */
int
p4est_sphere_match_exact (const p4est_sphere_t * box,
                          const p4est_sphere_t * sph, double t)
{
  const double        br = box->radius;
  const double        sr = sph->radius;
  int                 i;
  int                 outsi[P4EST_DIM];
  double              fdist[P4EST_DIM];
  double              m, ssmin, ssmax;
  double              rmin2, rmax2;

  P4EST_ASSERT (t >= 0.);

  /* squared inner and outer radius */
  rmin2 = (1. - t) * sr;
  rmin2 = SC_SQR (rmin2);
  rmax2 = (1. + t) * sr;
  rmax2 = SC_SQR (rmax2);

  /* figure out position of sphere's center relative to quadrant */
  for (i = 0; i < P4EST_DIM; ++i) {
    m = sph->center[i] - box->center[i];
    outsi[i] = (fdist[i] = fabs (m)) > br;
  }

  /* find minimum distance in all directions aligned with quadrant */
  /* find maximum distance to any corner in the same loop */
  ssmin = ssmax = 0.;
  for (i = 0; i < P4EST_DIM; ++i) {
    if (outsi[i]) {
      m = fdist[i] - br;
      P4EST_ASSERT (m >= 0.);
      ssmin += m * m;
    }
    m = fdist[i] + br;
    ssmax += m * m;
  }
  return ssmin <= rmax2 && rmin2 <= ssmax;
}
