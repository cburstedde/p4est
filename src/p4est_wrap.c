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

static void
init_callback (p4est_t * p4est, p4est_topidx_t which_tree,
               p4est_quadrant_t * q)
{
  q->p.user_int = 0;
}

static int
refine_callback (p4est_t * p4est, p4est_topidx_t which_tree,
                 p4est_quadrant_t *q)
{
  return q->p.user_int == P4EST_WRAP_REFINE;
}

static int
coarsen_callback (p4est_t * p4est, p4est_topidx_t which_tree,
                  p4est_quadrant_t *q[])
{
  int                 k;

  for (k = 0; k < P4EST_CHILDREN; ++k) {
    if (q[k]->p.user_int != P4EST_WRAP_COARSEN) {
      return 0;
    }
  }
  return 1;
}

static p4est_wrap_t *
p4est_wrap_new_conn (MPI_Comm mpicomm, p4est_connectivity_t *conn,
                     int initial_level)
{
  p4est_wrap_t       *pp;

  pp = SC_ALLOC (p4est_wrap_t, 1);
  pp->p4est_dim = P4EST_DIM;
  pp->p4est_half = P4EST_HALF;
  pp->p4est_faces = P4EST_FACES;
  pp->p4est_children = P4EST_CHILDREN;
  pp->conn = conn;
  pp->p4est = p4est_new_ext (mpicomm, pp->conn,
                             0, initial_level, 1, 0, init_callback, NULL);
  pp->flags = P4EST_ALLOC_ZERO (int8_t, pp->p4est->local_num_quadrants);

  pp->ghost = p4est_ghost_new (pp->p4est, P4EST_CONNECT_FULL);
  pp->mesh = p4est_mesh_new (pp->p4est, pp->ghost, P4EST_CONNECT_FULL);

  pp->ghost_aux = NULL;
  pp->mesh_aux = NULL;
  pp->match_aux = 0;

  return pp;
}

#ifndef P4_TO_P8

p4est_wrap_t       *
p4est_wrap_new_unitsquare (MPI_Comm mpicomm, int initial_level)
{
  return p4est_wrap_new_conn (mpicomm,
                              p4est_connectivity_new_unitsquare (),
                              initial_level);
}

p4est_wrap_t       *
p4est_wrap_new_periodic (MPI_Comm mpicomm, int initial_level)
{
  return p4est_wrap_new_conn (mpicomm,
                              p4est_connectivity_new_periodic (),
                              initial_level);
}

p4est_wrap_t       *
p4est_wrap_new_moebius (MPI_Comm mpicomm, int initial_level)
{
  return p4est_wrap_new_conn (mpicomm,
                              p4est_connectivity_new_moebius (),
                              initial_level);
}

#else

p8est_wrap_t       *
p8est_wrap_new_unitcube (MPI_Comm mpicomm, int initial_level)
{
  return p4est_wrap_new_conn (mpicomm,
                              p8est_connectivity_new_unitcube (),
                              initial_level);
}

#endif

p4est_wrap_t       *
p4est_wrap_new_world (int initial_level)
{
#ifndef P4_TO_P8
  return p4est_wrap_new_unitsquare (MPI_COMM_WORLD, initial_level);
#else
  return p8est_wrap_new_unitcube (MPI_COMM_WORLD, initial_level);
#endif
}

void
p4est_wrap_destroy (p4est_wrap_t * pp)
{
  if (pp->mesh_aux != NULL) {
    p4est_mesh_destroy (pp->mesh_aux);
  }
  if (pp->ghost_aux != NULL) {
    p4est_ghost_destroy (pp->ghost_aux);
  }
 
  p4est_mesh_destroy (pp->mesh);
  p4est_ghost_destroy (pp->ghost);
  
  P4EST_FREE (pp->flags);
  p4est_destroy (pp->p4est);
  p4est_connectivity_destroy (pp->conn);

  SC_FREE (pp);
}

int
p4est_wrap_refine (p4est_wrap_t * pp)
{
  int                 changed;
  size_t              allz, qz;
  p4est_locidx_t      tt;
  p4est_gloidx_t      global_num;
  p4est_t            *p4est = pp->p4est;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *q;
  
  P4EST_ASSERT (pp->mesh != NULL);
  P4EST_ASSERT (pp->ghost != NULL);
  P4EST_ASSERT (pp->mesh_aux == NULL);
  P4EST_ASSERT (pp->ghost_aux == NULL);
  P4EST_ASSERT (pp->match_aux == 0);

  allz = 0;
  for (tt = p4est->first_local_tree; tt <= p4est->last_local_tree; ++tt)
  {
    tree = p4est_tree_array_index (p4est->trees, tt);
    for (qz = 0; qz < tree->quadrants.elem_count; ++qz)
    {
      q = p4est_quadrant_array_index (&tree->quadrants, qz);
      q->p.user_int = (int) pp->flags[allz++];
    }
  }
  P4EST_ASSERT (allz == (size_t) p4est->local_num_quadrants);

  global_num = p4est->global_num_quadrants;
  p4est_refine (p4est, 0, refine_callback, init_callback);
  changed = global_num != p4est->global_num_quadrants;

  global_num = p4est->global_num_quadrants;
  p4est_coarsen (p4est, 0, coarsen_callback, init_callback);
  changed = changed || global_num != p4est->global_num_quadrants;

  if (changed) {
    P4EST_FREE (pp->flags);
    p4est_balance (p4est, P4EST_CONNECT_FULL, init_callback);
    pp->flags = P4EST_ALLOC_ZERO (int8_t, p4est->local_num_quadrants);

    pp->ghost_aux = p4est_ghost_new (p4est, P4EST_CONNECT_FULL);
    pp->mesh_aux = p4est_mesh_new (p4est, pp->ghost_aux, P4EST_CONNECT_FULL);
    pp->match_aux = 1;
  }

  return changed;
}

int
p4est_wrap_partition (p4est_wrap_t * pp)
{
  int                   changed;

  P4EST_ASSERT (pp->ghost != NULL);
  P4EST_ASSERT (pp->mesh != NULL);
  P4EST_ASSERT (pp->ghost_aux != NULL);
  P4EST_ASSERT (pp->mesh_aux != NULL);
  P4EST_ASSERT (pp->match_aux == 1);

  p4est_mesh_destroy (pp->mesh);
  p4est_ghost_destroy (pp->ghost);
  pp->match_aux = 0;
  
  /* In the future the flags could be used to pass partition weights */
  changed = p4est_partition_ext (pp->p4est, 1, NULL) > 0;

  if (changed) {
    P4EST_FREE (pp->flags);
    pp->flags = P4EST_ALLOC_ZERO (int8_t, pp->p4est->local_num_quadrants);

    pp->ghost = p4est_ghost_new (pp->p4est, P4EST_CONNECT_FULL);
    pp->mesh = p4est_mesh_new (pp->p4est, pp->ghost, P4EST_CONNECT_FULL);
  }
  else {
    pp->ghost = pp->ghost_aux;
    pp->mesh = pp->mesh_aux;
    pp->ghost_aux = NULL;
    pp->mesh_aux = NULL;
  }

  return changed;
}

void
p4est_wrap_complete (p4est_wrap_t * pp)
{
  P4EST_ASSERT (pp->ghost != NULL);
  P4EST_ASSERT (pp->mesh != NULL);
  P4EST_ASSERT (pp->ghost_aux != NULL);
  P4EST_ASSERT (pp->mesh_aux != NULL);
  P4EST_ASSERT (pp->match_aux == 0);
  
  p4est_mesh_destroy (pp->mesh_aux);
  p4est_ghost_destroy (pp->ghost_aux);
  pp->ghost_aux = NULL;
  pp->mesh_aux = NULL;
}

static p4est_wrap_leaf_t *
p4est_wrap_leaf_info (p4est_wrap_leaf_t * leaf)
{
#ifdef P4EST_DEBUG
  int                 nface;
  p4est_mesh_t       *mesh =
    leaf->pp->match_aux ? leaf->pp->mesh_aux : leaf->pp->mesh;
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
    printf ("C: Leaf face %d neighbor leaf %d\n", nface,
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
