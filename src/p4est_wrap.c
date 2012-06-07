/*
  This file is part of p4est.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2012 Carsten Burstedde

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
#include <p4est_extended.h>
#include <p4est_wrap.h>
#else
#include <p8est_bits.h>
#include <p8est_extended.h>
#include <p8est_wrap.h>
#endif

p4est_wrap_t       *
p4est_wrap_new (int initial_level)
{
  p4est_wrap_t       *pp;

  pp = SC_ALLOC (p4est_wrap_t, 1);
  pp->p4est_dim = P4EST_DIM;
  pp->p4est_half = P4EST_HALF;
  pp->p4est_faces = P4EST_FACES;
  pp->p4est_children = P4EST_CHILDREN;
#ifndef P4_TO_P8
  pp->conn = p4est_connectivity_new_unitsquare ();
#else
  pp->conn = p8est_connectivity_new_unitcube ();
#endif
  pp->p4est = p4est_new_ext (MPI_COMM_WORLD, pp->conn,
                             0, initial_level, 1, 0, NULL, NULL);
  pp->ghost = p4est_ghost_new (pp->p4est, P4EST_CONNECT_FULL);
  pp->mesh = p4est_mesh_new (pp->p4est, pp->ghost, P4EST_CONNECT_FULL);

  return pp;
}

void
p4est_wrap_destroy (p4est_wrap_t * pp)
{
  p4est_mesh_destroy (pp->mesh);
  p4est_ghost_destroy (pp->ghost);
  p4est_destroy (pp->p4est);
  p4est_connectivity_destroy (pp->conn);

  SC_FREE (pp);
}

static p4est_wrap_leaf_t *
p4est_wrap_leaf_info (p4est_wrap_leaf_t * leaf)
{
#ifdef P4EST_DEBUG
  int                 nface;
  p4est_mesh_t       *mesh = leaf->pp->mesh;
#endif
  p4est_quadrant_t    corner;

  leaf->total_quad = leaf->tree->quadrants_offset + leaf->which_quad;
  leaf->quad = p4est_quadrant_array_index (&leaf->tree->quadrants,
                                           leaf->which_quad);

  leaf->level = (int) leaf->quad->level;
  p4est_qcoord_to_vertex (leaf->pp->conn, leaf->which_tree,
                          leaf->quad->x, leaf->quad->y,
#ifdef P4_TO_P8
                          leaf->quad->z,
#endif
                          leaf->lowerleft);
  p4est_quadrant_corner_node (leaf->quad, P4EST_CHILDREN - 1, &corner);
  p4est_qcoord_to_vertex (leaf->pp->conn, leaf->which_tree,
                          corner.x, corner.y,
#ifdef P4_TO_P8
                          corner.z,
#endif
                          leaf->upperright);

#ifdef P4EST_DEBUG
  printf ("C: Leaf level %d tree %d tree_leaf %d local_leaf %d\n",
          leaf->level, leaf->which_tree, leaf->which_quad, leaf->total_quad);
  for (nface = 0; nface < P4EST_FACES; ++nface) {
    printf ("C: Leaf face %d leaf %d\n", nface,
            mesh->quad_to_quad[P4EST_FACES * leaf->total_quad + nface]);
  }
#endif

  return leaf;
}

p4est_wrap_leaf_t  *
p4est_wrap_leaf_first (p4est_wrap_t * pp)
{
  p4est_wrap_leaf_t  *leaf;
  p4est_t            *p4est = pp->p4est;

  if (p4est->local_num_quadrants == 0) {
    return NULL;
  }

  leaf = SC_ALLOC (p4est_wrap_leaf_t, 1);
  leaf->pp = pp;
  leaf->which_tree = p4est->first_local_tree;
  leaf->tree = p4est_tree_array_index (p4est->trees, leaf->which_tree);
  P4EST_ASSERT (leaf->tree->quadrants.elem_size > 0);
  leaf->which_quad = 0;

  return p4est_wrap_leaf_info (leaf);
}

p4est_wrap_leaf_t  *
p4est_wrap_leaf_next (p4est_wrap_leaf_t * leaf)
{
  p4est_t            *p4est = leaf->pp->p4est;

  P4EST_ASSERT (leaf != NULL);

  if ((size_t) leaf->which_quad + 1 == leaf->tree->quadrants.elem_count) {
    ++leaf->which_tree;
    if (leaf->which_tree > p4est->last_local_tree) {
      SC_FREE (leaf);
      return NULL;
    }
    leaf->tree = p4est_tree_array_index (p4est->trees, leaf->which_tree);
    P4EST_ASSERT (leaf->tree->quadrants.elem_size > 0);
    leaf->which_quad = 0;
  }
  else {
    ++leaf->which_quad;
  }

  return p4est_wrap_leaf_info (leaf);
}
