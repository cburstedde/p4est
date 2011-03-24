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

#ifndef P4EST_MESH_H
#define P4EST_MESH_H

#include <p8est.h>

SC_EXTERN_C_BEGIN;

/** This structure contains complete mesh information on the forest.
 * It stores the locally relevant neighborhood, that is, all locally owned
 * quadrants and one layer of adjacent ghost quadrants and their owners.
 *
 * All vertices of the locally owned quadrants are stored in xyz triples.
 * For each local quadrant, its eight corners point into the vertices array.
 * For each ghost quadrant, its owner rank is stored in ghost_to_proc.
 *
 * The quad_to_quad list stores one value for each local quadrant's face.
 * The quad_to_face list has equally many entries which are either:
 * 1. A value of v = 0..23 which indicates one same-size neighbor.
 *    This value is decoded as v = r * 6 + nf, where nf = 0..5 is the
 *    neigbbor's connecting face number and r = 0..3 is the relative
 *    orientation of the neighbor's face, see p8est_connectivity.h.
 * 2. A value of v = 24..119 which indicates a double-size neighbor.
 *    This value is decoded as v = 24 + h * 24 + r * 6 + nf, where
 *    r and nf are as above and h = 0..3 is the number of the subface.
 * 3. A value of v = -24..1 indicates four half-size neighbors.
 *    In this case the corresponding quad_to_quad index points into the
 *    quad_to_half table which stores four quadrant numbers per index,
 *    and the orientation of the smaller faces follow from 24 + v.
 */
typedef struct
{
  p4est_locidx_t      local_num_vertices;
  p4est_locidx_t      local_num_quadrants;
  p4est_locidx_t      ghost_num_quadrants;

  double             *vertices;
  p4est_locidx_t     *quad_to_vertex;   /* 8 indices for each local quad */
  int                *ghost_to_proc;    /* 1 integer for each ghost quad */

  p4est_locidx_t     *quad_to_quad;     /* 1 index for each of the 6 faces */
  int8_t             *quad_to_face;     /* encodes orientation/2:1 status */
  p4est_locidx_t     *quad_to_half;     /* stores half-size neigbors */
}
p8est_mesh_t;

SC_EXTERN_C_END;

#endif /* !P8EST_MESH_H */
