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

#ifndef P4_TO_P8
#include <p4est_geometry.h>
#else
#include <p8est_geometry.h>
#endif

void
p4est_geometry_destroy (p4est_geometry_t * geom)
{
  if (geom->destroy != NULL) {
    geom->destroy (geom);
  }
  else {
    P4EST_FREE (geom);
  }
}

static void
p4est_geometry_connectivity_X (p4est_geometry_t * geom,
                               p4est_topidx_t which_tree,
                               const double abc[3], double xyz[3])
{
  p4est_connectivity_t *connectivity = (p4est_connectivity_t *) geom->user;
  const p4est_topidx_t *tree_to_vertex = connectivity->tree_to_vertex;
  const double       *v = connectivity->vertices;
  double              eta_x, eta_y, eta_z = 0.;
  int                 j, k;
  p4est_topidx_t      vt[P4EST_CHILDREN];

  /* retrieve corners of the tree */
  for (k = 0; k < P4EST_CHILDREN; ++k) {
    vt[k] = tree_to_vertex[which_tree * P4EST_CHILDREN + k];
  }

  /* these are reference coordinates in [0, 1]**d */
  eta_x = abc[0];
  eta_y = abc[1];
  eta_z = abc[2];

  /* bi/trilinear transformation */
  for (j = 0; j < 3; ++j) {
    /* *INDENT-OFF* */
    xyz[j] =
           ((1. - eta_z) * ((1. - eta_y) * ((1. - eta_x) * v[3 * vt[0] + j] +
                                                  eta_x  * v[3 * vt[1] + j]) +
                                  eta_y  * ((1. - eta_x) * v[3 * vt[2] + j] +
                                                  eta_x  * v[3 * vt[3] + j]))
#ifdef P4_TO_P8
            +     eta_z  * ((1. - eta_y) * ((1. - eta_x) * v[3 * vt[4] + j] +
                                                  eta_x  * v[3 * vt[5] + j]) +
                                  eta_y  * ((1. - eta_x) * v[3 * vt[6] + j] +
                                                  eta_x  * v[3 * vt[7] + j]))
#endif
           );
    /* *INDENT-ON* */
  }
}

p4est_geometry_t   *
p4est_geometry_new_connectivity (p4est_connectivity_t * conn)
{
  p4est_geometry_t   *geom;

  P4EST_ASSERT (conn->vertices != NULL);

  geom = P4EST_ALLOC_ZERO (p4est_geometry_t, 1);

  geom->name = P4EST_STRING "_connectivity";
  geom->user = conn;
  geom->X = p4est_geometry_connectivity_X;

  return geom;
}
