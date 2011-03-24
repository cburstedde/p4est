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

#include <p4est_mesh.h>

p4est_mesh_t       *
p4est_mesh_new (p4est_t * p4est, p4est_ghost_t * ghost,
                p4est_balance_type_t btype)
{
  int                 rank;
  p4est_locidx_t      lq, ng;
  p4est_locidx_t      jl;
  p4est_mesh_t       *mesh;
  p4est_quadrant_t   *quad;

  mesh = P4EST_ALLOC (p4est_mesh_t, 1);

  mesh->local_num_vertices = 0;
  lq = mesh->local_num_quadrants = p4est->local_num_quadrants;
  ng = mesh->ghost_num_quadrants = (p4est_locidx_t) ghost->ghosts.elem_count;

  mesh->vertices = NULL;
  mesh->quad_to_vertex = P4EST_ALLOC (p4est_locidx_t, P4EST_CHILDREN * lq);
  mesh->ghost_to_proc = P4EST_ALLOC (int, ng);
  mesh->ghost_to_index = P4EST_ALLOC (p4est_locidx_t, ng);
  mesh->quad_to_quad = P4EST_ALLOC (p4est_locidx_t, P4EST_FACES * lq);
  mesh->quad_to_face = P4EST_ALLOC (int8_t, P4EST_FACES * lq);
  mesh->quad_to_half = NULL;

  /* Populate ghost information */
  rank = 0;
  for (jl = 0; jl < ng; ++jl) {
    while (ghost->proc_offsets[rank + 1] <= jl) {
      ++rank;
      P4EST_ASSERT (rank < p4est->mpisize);
    }
    quad = p4est_quadrant_array_index (&ghost->ghosts, (size_t) jl);
    mesh->ghost_to_proc[jl] = rank;
    mesh->ghost_to_index[jl] = quad->p.piggy3.local_num;
  }

  return mesh;
}

void
p4est_mesh_destroy (p4est_mesh_t * mesh)
{
  P4EST_FREE (mesh->vertices);
  P4EST_FREE (mesh->quad_to_vertex);
  P4EST_FREE (mesh->ghost_to_proc);
  P4EST_FREE (mesh->ghost_to_index);
  P4EST_FREE (mesh->quad_to_quad);
  P4EST_FREE (mesh->quad_to_face);
  P4EST_FREE (mesh->quad_to_half);
  P4EST_FREE (mesh);
}
