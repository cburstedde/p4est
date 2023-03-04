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

/** \file p4est_trimesh.h
 *
 * Generate a conforming triangle mesh from a 2:1 balanced p4est.
 * This mesh is represented by augmenting the \ref p4est_lnodes structure.
 */

#ifndef P4EST_TRIMESH_H
#define P4EST_TRIMESH_H

#include <p4est_lnodes.h>

SC_EXTERN_C_BEGIN;

typedef struct p4est_trimesh
{
  p4est_lnodes_t     *lnodes;
}
p4est_trimesh_t;

p4est_trimesh_t    *p4est_trimesh_new (p4est_t * p4est,
                                       p4est_ghost_t * ghost, int with_edge);

void                p4est_trimesh_destroy (p4est_trimesh_t * trimesh);

SC_EXTERN_C_END;

#endif /* !P4EST_TRIMESH_H */
