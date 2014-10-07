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

static int
refine_callback (p4est_t * p4est, p4est_topidx_t which_tree,
                 p4est_quadrant_t * q)
{
  p4est_wrap_t       *pp = (p4est_wrap_t *) p4est->user_pointer;
  const p4est_locidx_t old_counter = pp->inside_counter++;
  const uint8_t       flag = pp->flags[old_counter];

  P4EST_ASSERT (0 <= old_counter);
  P4EST_ASSERT (0 <= pp->num_replaced
                && pp->num_replaced <= pp->num_refine_flags);

  /* copy current flag since we cannot be certain that refinement occurs */
  pp->flags[old_counter] = 0;
  pp->temp_flags[old_counter + (P4EST_CHILDREN - 1) * pp->num_replaced] =
    flag & ~P4EST_WRAP_REFINE;

  return flag & P4EST_WRAP_REFINE;
}

static void
replace_on_refine (p4est_t * p4est, p4est_topidx_t which_tree,
                   int num_outgoing, p4est_quadrant_t * outgoing[],
                   int num_incoming, p4est_quadrant_t * incoming[])
{
  p4est_wrap_t       *pp = (p4est_wrap_t *) p4est->user_pointer;
  const p4est_locidx_t new_counter =
    pp->inside_counter - 1 + (P4EST_CHILDREN - 1) * pp->num_replaced++;
  const uint8_t       flag = pp->temp_flags[new_counter];
  int                 k;

  /* this function is only called when refinement actually happens */
  P4EST_ASSERT (num_outgoing == 1 && num_incoming == P4EST_CHILDREN);
  P4EST_ASSERT (!(flag & (P4EST_WRAP_REFINE | P4EST_WRAP_COARSEN)));

  for (k = 1; k < P4EST_CHILDREN; ++k) {
    pp->temp_flags[new_counter + k] = flag;
  }
}

static int
coarsen_callback (p4est_t * p4est, p4est_topidx_t which_tree,
                  p4est_quadrant_t * q[])
{
  p4est_wrap_t       *pp = (p4est_wrap_t *) p4est->user_pointer;
  const p4est_locidx_t old_counter = pp->inside_counter++;
  int                 k;

  /* are we not coarsening at all, just counting? */
  if (q[1] == NULL) {
    return 0;
  }

  /* now we are possibly coarsening */
  for (k = 0; k < P4EST_CHILDREN; ++k) {
    if (!(pp->temp_flags[old_counter + k] & P4EST_WRAP_COARSEN)) {
      return 0;
    }
  }

  /* we are definitely coarsening */
  pp->inside_counter += P4EST_CHILDREN - 1;
  ++pp->num_replaced;
  return 1;
}

p4est_wrap_t       *
p4est_wrap_new_conn (sc_MPI_Comm mpicomm, p4est_connectivity_t * conn,
                     int initial_level)
{
  p4est_wrap_t       *pp;

  pp = P4EST_ALLOC (p4est_wrap_t, 1);
  pp->user_pointer = NULL;

  pp->p4est_dim = P4EST_DIM;
  pp->p4est_half = P4EST_HALF;
  pp->p4est_faces = P4EST_FACES;
  pp->p4est_children = P4EST_CHILDREN;
  pp->conn = conn;
  pp->p4est = p4est_new_ext (mpicomm, pp->conn,
                             0, initial_level, 1, 0, NULL, NULL);
  pp->weight_exponent = 0;
  pp->flags = P4EST_ALLOC_ZERO (uint8_t, pp->p4est->local_num_quadrants);
  pp->temp_flags = NULL;
  pp->num_refine_flags = pp->inside_counter = pp->num_replaced = 0;

  pp->ghost = p4est_ghost_new (pp->p4est, P4EST_CONNECT_FULL);
  pp->mesh = p4est_mesh_new_ext (pp->p4est, pp->ghost, 1, 1,
                                 P4EST_CONNECT_FULL);

  pp->ghost_aux = NULL;
  pp->mesh_aux = NULL;
  pp->match_aux = 0;

  pp->p4est->user_pointer = pp;

  return pp;
}

#ifndef P4_TO_P8

p4est_wrap_t       *
p4est_wrap_new_unitsquare (sc_MPI_Comm mpicomm, int initial_level)
{
  return p4est_wrap_new_conn (mpicomm,
                              p4est_connectivity_new_unitsquare (),
                              initial_level);
}

p4est_wrap_t       *
p4est_wrap_new_periodic (sc_MPI_Comm mpicomm, int initial_level)
{
  return p4est_wrap_new_conn (mpicomm,
                              p4est_connectivity_new_periodic (),
                              initial_level);
}

p4est_wrap_t       *
p4est_wrap_new_rotwrap (sc_MPI_Comm mpicomm, int initial_level)
{
  return p4est_wrap_new_conn (mpicomm,
                              p4est_connectivity_new_rotwrap (),
                              initial_level);
}

p4est_wrap_t       *
p4est_wrap_new_corner (sc_MPI_Comm mpicomm, int initial_level)
{
  return p4est_wrap_new_conn (mpicomm,
                              p4est_connectivity_new_corner (),
                              initial_level);
}

p4est_wrap_t       *
p4est_wrap_new_pillow (sc_MPI_Comm mpicomm, int initial_level)
{
  return p4est_wrap_new_conn (mpicomm,
                              p4est_connectivity_new_pillow (),
                              initial_level);
}

p4est_wrap_t       *
p4est_wrap_new_moebius (sc_MPI_Comm mpicomm, int initial_level)
{
  return p4est_wrap_new_conn (mpicomm,
                              p4est_connectivity_new_moebius (),
                              initial_level);
}

p4est_wrap_t       *
p4est_wrap_new_cubed (sc_MPI_Comm mpicomm, int initial_level)
{
  return p4est_wrap_new_conn (mpicomm,
                              p4est_connectivity_new_cubed (), initial_level);
}

p4est_wrap_t       *
p4est_wrap_new_disk (sc_MPI_Comm mpicomm, int initial_level)
{
  return p4est_wrap_new_conn (mpicomm,
                              p4est_connectivity_new_disk (), initial_level);
}

#else

p8est_wrap_t       *
p8est_wrap_new_unitcube (sc_MPI_Comm mpicomm, int initial_level)
{
  return p4est_wrap_new_conn (mpicomm,
                              p8est_connectivity_new_unitcube (),
                              initial_level);
}

p8est_wrap_t       *
p8est_wrap_new_rotwrap (sc_MPI_Comm mpicomm, int initial_level)
{
  return p4est_wrap_new_conn (mpicomm,
                              p8est_connectivity_new_rotwrap (),
                              initial_level);
}

#endif

p4est_wrap_t       *
p4est_wrap_new_world (int initial_level)
{
#ifndef P4_TO_P8
  return p4est_wrap_new_unitsquare (sc_MPI_COMM_WORLD, initial_level);
#else
  return p8est_wrap_new_unitcube (sc_MPI_COMM_WORLD, initial_level);
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
  P4EST_FREE (pp->temp_flags);

  p4est_destroy (pp->p4est);
  p4est_connectivity_destroy (pp->conn);

  P4EST_FREE (pp);
}

p4est_ghost_t      *
p4est_wrap_get_ghost (p4est_wrap_t * pp)
{
  return pp->match_aux ? pp->ghost_aux : pp->ghost;
}

p4est_mesh_t       *
p4est_wrap_get_mesh (p4est_wrap_t * pp)
{
  return pp->match_aux ? pp->mesh_aux : pp->mesh;
}

void
p4est_wrap_mark_refine (p4est_wrap_t * pp,
                        p4est_topidx_t which_tree, p4est_locidx_t which_quad)
{
  p4est_t            *p4est = pp->p4est;
  p4est_tree_t       *tree;
  p4est_locidx_t      pos;
  uint8_t             flag;

  P4EST_ASSERT (p4est->first_local_tree <= which_tree);
  P4EST_ASSERT (which_tree <= p4est->last_local_tree);

  tree = p4est_tree_array_index (p4est->trees, which_tree);
  P4EST_ASSERT (0 <= which_quad &&
                which_quad < (p4est_locidx_t) tree->quadrants.elem_count);
  pos = tree->quadrants_offset + which_quad;
  P4EST_ASSERT (0 <= pos && pos < p4est->local_num_quadrants);

  flag = pp->flags[pos];
  if (!(flag & P4EST_WRAP_REFINE)) {
    pp->flags[pos] |= P4EST_WRAP_REFINE;
    ++pp->num_refine_flags;
  }
  pp->flags[pos] &= ~P4EST_WRAP_COARSEN;
}

void
p4est_wrap_mark_coarsen (p4est_wrap_t * pp,
                         p4est_topidx_t which_tree, p4est_locidx_t which_quad)
{
  p4est_t            *p4est = pp->p4est;
  p4est_tree_t       *tree;
  p4est_locidx_t      pos;
  uint8_t             flag;

  P4EST_ASSERT (p4est->first_local_tree <= which_tree);
  P4EST_ASSERT (which_tree <= p4est->last_local_tree);

  tree = p4est_tree_array_index (p4est->trees, which_tree);
  P4EST_ASSERT (0 <= which_quad &&
                which_quad < (p4est_locidx_t) tree->quadrants.elem_count);
  pos = tree->quadrants_offset + which_quad;
  P4EST_ASSERT (0 <= pos && pos < p4est->local_num_quadrants);

  flag = pp->flags[pos];
  if (flag & P4EST_WRAP_REFINE) {
    pp->flags[pos] &= ~P4EST_WRAP_REFINE;
    --pp->num_refine_flags;
  }
  pp->flags[pos] |= P4EST_WRAP_COARSEN;
}

int
p4est_wrap_adapt (p4est_wrap_t * pp)
{
  int                 changed;
#ifdef P4EST_ENABLE_DEBUG
  p4est_locidx_t      jl, local_num;
#endif
  p4est_gloidx_t      global_num;
  p4est_t            *p4est = pp->p4est;

  P4EST_ASSERT (pp->mesh != NULL);
  P4EST_ASSERT (pp->ghost != NULL);
  P4EST_ASSERT (pp->mesh_aux == NULL);
  P4EST_ASSERT (pp->ghost_aux == NULL);
  P4EST_ASSERT (pp->match_aux == 0);

  P4EST_ASSERT (pp->temp_flags == NULL);
  P4EST_ASSERT (pp->num_refine_flags >= 0 &&
                pp->num_refine_flags <= p4est->local_num_quadrants);

  /* This allocation is optimistic when not all refine requests are honored */
  pp->temp_flags = P4EST_ALLOC_ZERO (uint8_t, p4est->local_num_quadrants +
                                     (P4EST_CHILDREN - 1) *
                                     pp->num_refine_flags);

  /* Execute refinement */
  pp->inside_counter = pp->num_replaced = 0;
#ifdef P4EST_ENABLE_DEBUG
  local_num = p4est->local_num_quadrants;
#endif
  global_num = p4est->global_num_quadrants;
  p4est_refine_ext (p4est, 0, -1, refine_callback, NULL, replace_on_refine);
  P4EST_ASSERT (pp->inside_counter == local_num);
  P4EST_ASSERT (p4est->local_num_quadrants - local_num ==
                pp->num_replaced * (P4EST_CHILDREN - 1));
  changed = global_num != p4est->global_num_quadrants;

  /* Execute coarsening */
  pp->inside_counter = pp->num_replaced = 0;
#ifdef P4EST_ENABLE_DEBUG
  local_num = p4est->local_num_quadrants;
#endif
  global_num = p4est->global_num_quadrants;
  p4est_coarsen_ext (p4est, 0, 1, coarsen_callback, NULL, NULL);
  P4EST_ASSERT (pp->inside_counter == local_num);
  P4EST_ASSERT (local_num - p4est->local_num_quadrants ==
                pp->num_replaced * (P4EST_CHILDREN - 1));
  changed = changed || global_num != p4est->global_num_quadrants;

  /* Free temporary flags */
  P4EST_FREE (pp->temp_flags);
  pp->temp_flags = NULL;

  /* Only if refinement and/or coarsening happened do we need to balance */
  if (changed) {
    P4EST_FREE (pp->flags);
    p4est_balance (p4est, P4EST_CONNECT_FULL, NULL);
    pp->flags = P4EST_ALLOC_ZERO (uint8_t, p4est->local_num_quadrants);

    pp->ghost_aux = p4est_ghost_new (p4est, P4EST_CONNECT_FULL);
    pp->mesh_aux = p4est_mesh_new_ext (p4est, pp->ghost_aux, 1, 1,
                                       P4EST_CONNECT_FULL);
    pp->match_aux = 1;
  }
#ifdef P4EST_ENABLE_DEBUG
  else {
    for (jl = 0; jl < p4est->local_num_quadrants; ++jl) {
      P4EST_ASSERT (pp->flags[jl] == 0);
    }
  }
#endif
  pp->num_refine_flags = 0;

  return changed;
}

static int
partition_weight (p4est_t * p4est, p4est_topidx_t which_tree,
                  p4est_quadrant_t * quadrant)
{
  p4est_wrap_t       *pp = (p4est_wrap_t *) p4est->user_pointer;

  return 1 << ((int) quadrant->level * pp->weight_exponent);
}

int
p4est_wrap_partition (p4est_wrap_t * pp, int weight_exponent)
{
  int                 changed;

  P4EST_ASSERT (pp->ghost != NULL);
  P4EST_ASSERT (pp->mesh != NULL);
  P4EST_ASSERT (pp->ghost_aux != NULL);
  P4EST_ASSERT (pp->mesh_aux != NULL);
  P4EST_ASSERT (pp->match_aux == 1);

  p4est_mesh_destroy (pp->mesh);
  p4est_ghost_destroy (pp->ghost);
  pp->match_aux = 0;

  /* In the future the flags could be used to pass partition weights */
  /* We need to lift the restriction on 64 bits for the global weight sum */
  P4EST_ASSERT (weight_exponent == 0 || weight_exponent == 1);
  pp->weight_exponent = weight_exponent;
  changed =
    p4est_partition_ext (pp->p4est, 1,
                         weight_exponent ? partition_weight : NULL) > 0;

  if (changed) {
    P4EST_FREE (pp->flags);
    pp->flags = P4EST_ALLOC_ZERO (uint8_t, pp->p4est->local_num_quadrants);

    pp->ghost = p4est_ghost_new (pp->p4est, P4EST_CONNECT_FULL);
    pp->mesh = p4est_mesh_new_ext (pp->p4est, pp->ghost, 1, 1,
                                   P4EST_CONNECT_FULL);
  }
  else {
    memset (pp->flags, 0, sizeof (uint8_t) * pp->p4est->local_num_quadrants);

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
#if 0
#ifdef P4EST_ENABLE_DEBUG
  int                 nface;
  p4est_mesh_t       *mesh = p4est_wrap_get_mesh (leaf->pp);
#endif
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

#if 0
#ifdef P4EST_ENABLE_DEBUG
  printf ("C: Leaf level %d tree %d tree_leaf %d local_leaf %d\n",
          leaf->level, leaf->which_tree, leaf->which_quad, leaf->total_quad);
  for (nface = 0; nface < P4EST_FACES; ++nface) {
    printf ("C: Leaf face %d neighbor leaf %d\n", nface,
            mesh->quad_to_quad[P4EST_FACES * leaf->total_quad + nface]);
  }
#endif
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

  leaf = P4EST_ALLOC (p4est_wrap_leaf_t, 1);
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
      P4EST_FREE (leaf);
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
