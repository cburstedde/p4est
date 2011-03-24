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

#include <p4est_ghost.h>

SC_EXTERN_C_BEGIN;

typedef struct
{
  p4est_locidx_t      local_num_vertices;
  p4est_locidx_t      local_num_quadrants;
  p4est_locidx_t      ghost_num_quadrants;

  double             *vertices;
  p4est_locidx_t     *quad_to_vertex;   /* 4 indices for each local quad */
  int                *ghost_to_proc;    /* 1 integer for each ghost quad */
  p4est_locidx_t     *ghost_to_index;   /* 1 remote index for each ghost */

  p4est_locidx_t     *quad_to_quad;     /* 1 index for each of the 4 faces */
  int8_t             *quad_to_face;     /* encodes orientation/2:1 status */
  sc_array_t         *quad_to_half;     /* stores half-size neigbors */
}
p4est_mesh_t;

p4est_mesh_t       *p4est_mesh_new (p4est_t * p4est, p4est_ghost_t * ghost,
                                    p4est_balance_type_t btype);
void                p4est_mesh_destroy (p4est_mesh_t * mesh);

SC_EXTERN_C_END;

#endif /* !P4EST_MESH_H */
