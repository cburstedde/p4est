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

/** This structure contains complete mesh information on the forest.
 * It stores the locally relevant neighborhood, that is, all locally owned
 * quadrants and one layer of adjacent ghost quadrants and their owners.
 *
 * All vertices of the locally owned quadrants are stored in xyz triples.
 * For each local quadrant, its four corners point into the vertices array.
 * Some vertices are duplicated and no effort is made to uniquify them.
 * For each ghost quadrant, its owner rank is stored in ghost_to_proc,
 * and its number in it's owners range of local quadrants in ghost_to_index.
 *
 * The quad_to_quad list stores one value for each local quadrant's face.
 * This value is in 0..local_num_quadrants-1 for local quadrants, or in
 * local_num_quadrants + (0..ghost_num_quadrants-1) for ghost quadrants.
 * The quad_to_face list has equally many entries which are either:
 * 1. A value of v = 0..7 indicates one same-size neighbor.
 *    This value is decoded as v = r * 4 + nf, where nf = 0..3 is the
 *    neigbbor's connecting face number and r = 0..1 is the relative
 *    orientation of the neighbor's face, see p4est_connectivity.h.
 * 2. A value of v = 8..23 indicates a double-size neighbor.
 *    This value is decoded as v = 8 + h * 8 + r * 4 + nf, where
 *    r and nf are as above and h = 0..1 is the number of the subface.
 * 3. A value of v = -8..-1 indicates two half-size neighbors.
 *    In this case the corresponding quad_to_quad index points into the
 *    quad_to_half array which stores two quadrant numbers per index,
 *    and the orientation of the smaller faces follows from 8 + v.
 *    The entries of quad_to_half encode between local and ghost quadrant
 *    in the same way as the quad_to_quad values described above.
 * A quadrant on the boundary of the forest sees itself and its face number.
 */
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

/** Calculate the memory usage of the mesh structure.
 * \param [in] mesh     Mesh structure.
 * \return              Memory used in bytes.
 */
size_t              p4est_mesh_memory_used (p4est_mesh_t * mesh);

/** Create a p4est_mesh structure.
 * The vertex information will be filled if p4est->connectivity contains
 * vertices.  Currently only face neighborhood information is stored.
 * \param [in] p4est    A forest that is fully 2:1 balanced.
 * \param [in] ghost    The ghost layer created from the provided p4est.
 * \param [in] btype    Currently ignored, only face neighbors are stored.
 * \return              A fully allocated mesh structure.
 */
p4est_mesh_t       *p4est_mesh_new (p4est_t * p4est, p4est_ghost_t * ghost,
                                    p4est_connect_type_t btype);

/** Destroy a p4est_mesh structure.
 * \param [in] mesh     Mesh structure previously created by p4est_mesh_new.
 */
void                p4est_mesh_destroy (p4est_mesh_t * mesh);

SC_EXTERN_C_END;

#endif /* !P4EST_MESH_H */
