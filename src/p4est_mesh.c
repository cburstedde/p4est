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

#ifndef P4_TO_P8
#include <p4est_iterate.h>
#include <p4est_mesh.h>
#include <p4est_nodes.h>
#else
#include <p8est_iterate.h>
#include <p8est_mesh.h>
#include <p8est_nodes.h>
#endif

static void
mesh_iter_face (p4est_iter_face_info_t * info, void *user_data)
{
  int                 h;
  int                 swapsides;
  p4est_mesh_t       *mesh = (p4est_mesh_t *) user_data;
  p4est_locidx_t      jl, jl2, jls[P4EST_HALF];
  p4est_locidx_t      in_qtoq, halfindex;
  p4est_locidx_t     *halfentries;
  p4est_tree_t       *tree;
  p4est_iter_face_side_t *side, *side2, *tempside;

  if (info->sides.elem_count == 1) {
    /* this face is on an outside boundary of the forest */
    P4EST_ASSERT (info->orientation == 0);
    P4EST_ASSERT (info->tree_boundary);
    side = (p4est_iter_face_side_t *) sc_array_index (&info->sides, 0);
    P4EST_ASSERT (0 <= side->treeid &&
                  side->treeid < info->p4est->connectivity->num_trees);
    P4EST_ASSERT (0 <= side->face && side->face < P4EST_FACES);
    P4EST_ASSERT (!side->is_hanging && !side->is.full.is_ghost);
    tree = p4est_tree_array_index (info->p4est->trees, side->treeid);
    jl = side->is.full.quadid + tree->quadrants_offset;
    P4EST_ASSERT (0 <= jl && jl < mesh->local_num_quadrants);
    in_qtoq = P4EST_FACES * jl + side->face;
    mesh->quad_to_quad[in_qtoq] = jl;   /* put in myself and my own face */
    mesh->quad_to_face[in_qtoq] = side->face;
  }
  else {
    /* this face is between two quadrants */
    P4EST_ASSERT (info->orientation == 0 || info->tree_boundary);
    P4EST_ASSERT (info->sides.elem_count == 2);
    side = (p4est_iter_face_side_t *) sc_array_index (&info->sides, 0);
    side2 = (p4est_iter_face_side_t *) sc_array_index (&info->sides, 1);
    P4EST_ASSERT (info->tree_boundary || side->treeid == side2->treeid);
    P4EST_ASSERT (!side->is_hanging || !side2->is_hanging);
    if (!side->is_hanging && !side2->is_hanging) {
      /* same-size face neighbors */
      P4EST_ASSERT (!side->is.full.is_ghost || !side2->is.full.is_ghost);

      /* determine both quadrant numbers */
      if (!side->is.full.is_ghost) {
        tree = p4est_tree_array_index (info->p4est->trees, side->treeid);
        jl = side->is.full.quadid + tree->quadrants_offset;
        P4EST_ASSERT (0 <= jl && jl < mesh->local_num_quadrants);
      }
      else {
        jl = mesh->local_num_quadrants + side->is.full.quadid;
      }
      if (!side2->is.full.is_ghost) {
        tree = p4est_tree_array_index (info->p4est->trees, side2->treeid);
        jl2 = side2->is.full.quadid + tree->quadrants_offset;
        P4EST_ASSERT (0 <= jl2 && jl2 < mesh->local_num_quadrants);
      }
      else {
        jl2 = mesh->local_num_quadrants + side2->is.full.quadid;
      }

      /* encode quadrant neighborhood */
      if (!side->is.full.is_ghost) {
        in_qtoq = P4EST_FACES * jl + side->face;
        P4EST_ASSERT (mesh->quad_to_quad[in_qtoq] == -1);
        P4EST_ASSERT (mesh->quad_to_face[in_qtoq] == -25);
        mesh->quad_to_quad[in_qtoq] = jl2;
        mesh->quad_to_face[in_qtoq] =
          P4EST_FACES * info->orientation + side2->face;
      }
      if (!side2->is.full.is_ghost) {
        in_qtoq = P4EST_FACES * jl2 + side2->face;
        P4EST_ASSERT (mesh->quad_to_quad[in_qtoq] == -1);
        P4EST_ASSERT (mesh->quad_to_face[in_qtoq] == -25);
        mesh->quad_to_quad[in_qtoq] = jl;
        mesh->quad_to_face[in_qtoq] =
          P4EST_FACES * info->orientation + side->face;
      }
    }
    else {
      /* one of the faces is hanging, rename so it's always side2 */
      swapsides = side->is_hanging;
      if (swapsides) {
        tempside = side;
        side = side2;
        side2 = tempside;
      }
      P4EST_ASSERT (!side->is_hanging && side2->is_hanging);

      /* determine quadrant number for non-hanging large face */
      if (!side->is.full.is_ghost) {
        tree = p4est_tree_array_index (info->p4est->trees, side->treeid);
        jl = side->is.full.quadid + tree->quadrants_offset;
        P4EST_ASSERT (0 <= jl && jl < mesh->local_num_quadrants);
      }
      else {
        jl = mesh->local_num_quadrants + side->is.full.quadid;
      }

      /* determine quadrant numbers for all hanging faces */
      for (h = 0; h < P4EST_HALF; ++h) {
        if (!side2->is.hanging.is_ghost[h]) {
          tree = p4est_tree_array_index (info->p4est->trees, side2->treeid);
          jls[h] = side2->is.hanging.quadid[h] + tree->quadrants_offset;
          P4EST_ASSERT (0 <= jls[h] && jls[h] < mesh->local_num_quadrants);
        }
        else {
          jls[h] = mesh->local_num_quadrants + side2->is.hanging.quadid[h];
        }
      }

      /* encode quadrant neighborhood */
      if (!side->is.full.is_ghost) {
        in_qtoq = P4EST_FACES * jl + side->face;
        P4EST_ASSERT (mesh->quad_to_quad[in_qtoq] == -1);
        P4EST_ASSERT (mesh->quad_to_face[in_qtoq] == -25);
        halfindex = (p4est_locidx_t) mesh->quad_to_half->elem_count;
        mesh->quad_to_quad[in_qtoq] = halfindex;
        mesh->quad_to_face[in_qtoq] =
          P4EST_FACES * (info->orientation - P4EST_HALF) + side2->face;
        halfentries = (p4est_locidx_t *) sc_array_push (mesh->quad_to_half);
        for (h = 0; h < P4EST_HALF; ++h) {
          halfentries[h] = jls[h];
        }
      }
      for (h = 0; h < P4EST_HALF; ++h) {
        if (!side2->is.hanging.is_ghost[h]) {
          in_qtoq = P4EST_FACES * jls[h] + side2->face;
          P4EST_ASSERT (mesh->quad_to_quad[in_qtoq] == -1);
          P4EST_ASSERT (mesh->quad_to_face[in_qtoq] == -25);
          mesh->quad_to_quad[in_qtoq] = jl;
          mesh->quad_to_face[in_qtoq] =
            P4EST_FACES * (info->orientation + (h + 1) * P4EST_HALF) +
            side->face;
        }
      }
    }
  }
}

size_t
p4est_mesh_memory_used (p4est_mesh_t * mesh)
{
  size_t              lqz, ngz;

  lqz = (size_t) mesh->local_num_quadrants;
  ngz = (size_t) mesh->ghost_num_quadrants;

  return sizeof (p4est_mesh_t) +
    ((P4EST_CHILDREN + P4EST_FACES) * lqz + ngz) * sizeof (p4est_locidx_t) +
    ngz * sizeof (int) + (P4EST_FACES * lqz) * sizeof (int8_t) +
    3 * (size_t) mesh->local_num_vertices * sizeof (double) +
    sc_array_memory_used (mesh->quad_to_half, 1);
}

p4est_mesh_t       *
p4est_mesh_new (p4est_t * p4est, p4est_ghost_t * ghost,
                p4est_balance_type_t btype)
{
  int                 rank;
  size_t              ncz, zz;
  p4est_locidx_t      lq, ng;
  p4est_locidx_t      jl;
  p4est_mesh_t       *mesh;
  p4est_nodes_t      *nodes;
  p4est_quadrant_t   *quad;

  P4EST_ASSERT (p4est_is_balanced (p4est, P4EST_BALANCE_FULL));

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
  mesh->quad_to_half = sc_array_new (P4EST_HALF * sizeof (p4est_locidx_t));

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

  /* Fill arrays with default values */
  memset (mesh->quad_to_quad, -1, P4EST_FACES * lq * sizeof (p4est_locidx_t));
  memset (mesh->quad_to_face, -25, P4EST_FACES * lq * sizeof (int8_t));

  /* Call the forest iterator to collect face connectivity */
  p4est_iterate (p4est, ghost, mesh, NULL, mesh_iter_face,
#ifdef P4_TO_P8
                 NULL,
#endif
                 NULL);

  /* Construct non-unique vertex information if vertices are present */
  if (p4est->connectivity->num_vertices > 0) {
    P4EST_ASSERT (p4est->connectivity->vertices != NULL);
    nodes = p4est_nodes_new (p4est, NULL);      /* fast local version */
    P4EST_ASSERT (nodes->num_local_quadrants == mesh->local_num_quadrants);
    P4EST_ASSERT (nodes->num_owned_shared == 0);
    P4EST_ASSERT (nodes->face_hangings.elem_count == 0);
#ifdef P4_TO_P8
    P4EST_ASSERT (nodes->edge_hangings.elem_count == 0);
#endif
    P4EST_ASSERT (nodes->shared_indeps.elem_count == 0);
    memcpy (mesh->quad_to_vertex, nodes->local_nodes,
            P4EST_CHILDREN * lq * sizeof (p4est_locidx_t));
    ncz = (p4est_locidx_t) nodes->indep_nodes.elem_count;
    mesh->vertices = P4EST_ALLOC (double, 3 * ncz);
    mesh->local_num_vertices = (p4est_locidx_t) ncz;
    for (zz = 0; zz < ncz; ++zz) {
      quad = p4est_quadrant_array_index (&nodes->indep_nodes, zz);
      p4est_qcoord_to_vertex (p4est->connectivity,
                              quad->p.which_tree, quad->x, quad->y,
#ifdef P4_TO_P8
                              quad->z,
#endif
                              mesh->vertices + 3 * zz);
    }
    p4est_nodes_destroy (nodes);
  }
  else {
    memset (mesh->quad_to_vertex, -1,
            P4EST_CHILDREN * lq * sizeof (p4est_locidx_t));
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
  sc_array_destroy (mesh->quad_to_half);
  P4EST_FREE (mesh);
}
