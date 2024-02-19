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

#ifndef P4_TO_P8
#include <p4est_algorithms.h>
#include <p4est_bits.h>
#include <p4est_wrap.h>
#else
#include <p8est_algorithms.h>
#include <p8est_bits.h>
#include <p8est_wrap.h>
#endif

static int
refine_callback (p4est_t * p4est, p4est_topidx_t which_tree,
                 p4est_quadrant_t * q)
{
  p4est_wrap_t       *pp = (p4est_wrap_t *) p4est->user_pointer;
  const p4est_locidx_t old_counter = pp->inside_counter++;
  const uint8_t       flag = pp->flags[old_counter];

  P4EST_ASSERT (pp->params.coarsen_delay >= 0);
  P4EST_ASSERT (0 <= old_counter);
  P4EST_ASSERT (0 <= pp->num_replaced
                && pp->num_replaced <= pp->num_refine_flags);

  /* copy current flag since we cannot be certain that refinement occurs */
  pp->flags[old_counter] = 0;
  pp->temp_flags[old_counter + (P4EST_CHILDREN - 1) * pp->num_replaced] =
    flag & ~P4EST_WRAP_REFINE;

  /* increase quadrant's counter of most recent adaptation */
  /* if refinement actually occurs, it will be reset to zero in all children */
  if (pp->params.coarsen_delay && q->p.user_int >= 0) {
    ++q->p.user_int;
  }

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

  /* we have set the first flag in the refinement callback, do the others */
  for (k = 1; k < P4EST_CHILDREN; ++k) {
    pp->temp_flags[new_counter + k] = flag;
  }

  /* reset the counter for most recent adaptation */
  P4EST_ASSERT (pp->params.coarsen_delay >= 0);
  if (pp->params.coarsen_delay) {
    for (k = 0; k < P4EST_CHILDREN; ++k) {
      incoming[k]->p.user_int = 0;
    }
  }

  /* pass the replaced quadrants to the user-provided function */
  if (pp->params.replace_fn != NULL) {
    pp->params.replace_fn (p4est, which_tree,
                    num_outgoing, outgoing, num_incoming, incoming);
  }
}

static int
coarsen_callback (p4est_t * p4est, p4est_topidx_t which_tree,
                  p4est_quadrant_t * q[])
{
  p4est_wrap_t       *pp = (p4est_wrap_t *) p4est->user_pointer;
  const p4est_locidx_t old_counter = pp->inside_counter++;
  int                 k;

  P4EST_ASSERT (pp->params.coarsen_delay >= 0);

  /* are we not coarsening at all, just counting? */
  if (q[1] == NULL) {
    return 0;
  }

  /* now we are possibly coarsening */
  for (k = 0; k < P4EST_CHILDREN; ++k) {
    if (!(pp->temp_flags[old_counter + k] & P4EST_WRAP_COARSEN)) {
      /* coarsening flag was not set */
      return 0;
    }
    if (pp->params.coarsen_delay && q[k]->p.user_int >= 0 &&
        q[k]->p.user_int <= pp->params.coarsen_delay) {
      /* most recent adaptation has been too recent */
      return 0;
    }
  }

  /* we are definitely coarsening */
  pp->inside_counter += P4EST_CHILDREN - 1;
  ++pp->num_replaced;
  return 1;
}

static void
replace_on_coarsen (p4est_t * p4est, p4est_topidx_t which_tree,
                    int num_outgoing, p4est_quadrant_t * outgoing[],
                    int num_incoming, p4est_quadrant_t * incoming[])
{
  p4est_wrap_t       *pp = (p4est_wrap_t *) p4est->user_pointer;
  P4EST_ASSERT (num_incoming == 1 && num_outgoing == P4EST_CHILDREN);
  P4EST_ASSERT (pp->params.coarsen_delay > 0);

  /* reset most recent adaptation timer */
  incoming[0]->p.user_int = pp->params.coarsen_affect ? 0 : -1;

  /* pass the replaced quadrants to the user-provided function */
  if (pp->params.replace_fn != NULL) {
    pp->params.replace_fn (p4est, which_tree,
                    num_outgoing, outgoing, num_incoming, incoming);
  }
}

static void
replace_on_balance (p4est_t * p4est, p4est_topidx_t which_tree,
                    int num_outgoing, p4est_quadrant_t * outgoing[],
                    int num_incoming, p4est_quadrant_t * incoming[])
{
  p4est_wrap_t       *pp = (p4est_wrap_t *) p4est->user_pointer;
  int                 k;

  /* this function is called when refinement occurs in balance */
  P4EST_ASSERT (num_outgoing == 1 && num_incoming == P4EST_CHILDREN);
  P4EST_ASSERT (pp->params.coarsen_delay > 0);

  /* negative value means coarsening is allowed next time */
  for (k = 0; k < P4EST_CHILDREN; ++k) {
    incoming[k]->p.user_int = -1;
  }

  /* pass the replaced quadrants to the user-provided function */
  if (pp->params.replace_fn != NULL) {
    pp->params.replace_fn (p4est, which_tree,
                           num_outgoing, outgoing, num_incoming, incoming);
  }
}

void
p4est_wrap_params_init (p4est_wrap_params_t *params)
{
  memset (params, 0, sizeof (p4est_wrap_params_t));

  params->hollow = 1;
  p4est_mesh_params_init (&params->mesh_params);
  params->replace_fn = NULL;
  params->coarsen_delay = 0;
  params->coarsen_affect = 0;
  params->partition_for_coarsening = 1;
  params->user_pointer = NULL;
}

p4est_wrap_t       *
p4est_wrap_new_conn (sc_MPI_Comm mpicomm, p4est_connectivity_t * conn,
                     int initial_level)
{
  p4est_wrap_params_t params;

  p4est_wrap_params_init (&params);
  params.hollow = 0;
  params.mesh_params.btype = P4EST_CONNECT_FULL;
  params.mesh_params.compute_level_lists = 1;
  params.mesh_params.compute_tree_index = 1;

  return p4est_wrap_new_params (mpicomm, conn, initial_level, &params);
}

p4est_wrap_t       *
p4est_wrap_new_p4est (p4est_t * p4est, int hollow, p4est_connect_type_t btype,
                      p4est_replace_t replace_fn, void *user_pointer)
{
  p4est_wrap_params_t params;

  p4est_wrap_params_init (&params);
  params.hollow = hollow;
  params.mesh_params.btype = btype;
  params.mesh_params.compute_level_lists = 1;
  params.mesh_params.compute_tree_index = 1;
  params.replace_fn = replace_fn;
  params.user_pointer = user_pointer;

  return p4est_wrap_new_p4est_params (p4est, &params);
}

p4est_wrap_t       *
p4est_wrap_new_p4est_params (p4est_t * p4est, p4est_wrap_params_t * params)
{
  p4est_wrap_t       *pp;

  P4EST_ASSERT (p4est_is_valid (p4est));
  P4EST_ASSERT (p4est->user_pointer == NULL);

  pp = P4EST_ALLOC_ZERO (p4est_wrap_t, 1);

  /* store wrap parameters in wrap */
  if (params != NULL) {
    pp->params = *params;
    params = NULL;
  }
  else {
    p4est_wrap_params_init (&pp->params);
  }

  sc_refcount_init (&pp->conn_rc, p4est_package_id);
  pp->conn = p4est->connectivity;
  pp->conn_owner = NULL;

  pp->p4est_dim = P4EST_DIM;
  pp->p4est_half = P4EST_HALF;
  pp->p4est_faces = P4EST_FACES;
  pp->p4est_children = P4EST_CHILDREN;
  pp->p4est = p4est;
  pp->weight_exponent = 0;      /* keep this even though using ALLOC_ZERO */

  if (!pp->params.hollow) {
    pp->flags = P4EST_ALLOC_ZERO (uint8_t, pp->p4est->local_num_quadrants);
    pp->ghost = p4est_ghost_new (pp->p4est, pp->params.mesh_params.btype);
    pp->mesh =
      p4est_mesh_new_params (pp->p4est, pp->ghost, &pp->params.mesh_params);
  }

  /* reset the data size since changing the p4est_wrap will affect p.user_int */
  p4est_reset_data (pp->p4est, 0, NULL, NULL);

  pp->p4est->user_pointer = pp;

  return pp;
}

p4est_wrap_t       *
p4est_wrap_new_ext (sc_MPI_Comm mpicomm, p4est_connectivity_t * conn,
                    int initial_level, int hollow, p4est_connect_type_t btype,
                    p4est_replace_t replace_fn, void *user_pointer)
{
  p4est_wrap_params_t params;

  p4est_wrap_params_init (&params);
  params.hollow = hollow;
  params.mesh_params.btype = btype;
  params.mesh_params.compute_level_lists = 1;
  params.mesh_params.compute_tree_index = 1;
  params.replace_fn = replace_fn;
  params.user_pointer = user_pointer;

  return p4est_wrap_new_params (mpicomm, conn, initial_level, &params);
}

p4est_wrap_t       *
p4est_wrap_new_params (sc_MPI_Comm mpicomm, p4est_connectivity_t * conn,
                       int initial_level, p4est_wrap_params_t * params)
{
  P4EST_ASSERT (p4est_connectivity_is_valid (conn));
  return p4est_wrap_new_p4est_params (p4est_new_ext (mpicomm, conn,
                                                     0, initial_level, 1,
                                                     0, NULL, NULL), params);
}

p4est_wrap_t       *
p4est_wrap_new_copy (p4est_wrap_t * source, size_t data_size,
                     p4est_replace_t replace_fn, void *user_pointer)
{
  p4est_wrap_t       *pp;

  P4EST_ASSERT (source != NULL);

  pp = P4EST_ALLOC_ZERO (p4est_wrap_t, 1);

  /* copy the sources wrap paramters; however the copy will is hollow */
  pp->params = source->params;
  pp->params.hollow = 1;

  sc_refcount_init_invalid (&pp->conn_rc);
  pp->conn_owner = (source->conn_owner != NULL ? source->conn_owner : source);
  pp->conn = pp->conn_owner->conn;
  sc_refcount_ref (&pp->conn_owner->conn_rc);

  pp->p4est_dim = P4EST_DIM;
  pp->p4est_half = P4EST_HALF;
  pp->p4est_faces = P4EST_FACES;
  pp->p4est_children = P4EST_CHILDREN;
  pp->params.replace_fn = replace_fn;
  pp->p4est = p4est_copy (source->p4est, 0);
  if (data_size > 0) {
    p4est_reset_data (pp->p4est, data_size, NULL, NULL);
  }

  pp->weight_exponent = 0;      /* keep this even though using ALLOC_ZERO */

  pp->p4est->user_pointer = pp;
  pp->params.user_pointer = user_pointer;

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
p4est_wrap_new_disk (sc_MPI_Comm mpicomm, int px, int py, int initial_level)
{
  return p4est_wrap_new_conn
    (mpicomm, p4est_connectivity_new_disk (px, py), initial_level);
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
p4est_wrap_new_brick (sc_MPI_Comm mpicomm, int bx, int by,
#ifdef P4_TO_P8
                      int bz,
#endif
                      int px, int py,
#ifdef P4_TO_P8
                      int pz,
#endif
                      int initial_level)
{
  P4EST_ASSERT (bx > 0 && by > 0);
#ifdef P4_TO_P8
  P4EST_ASSERT (bz > 0);
#endif
  return p4est_wrap_new_conn (mpicomm, p4est_connectivity_new_brick (bx, by,
#ifdef P4_TO_P8
                                                                     bz,
#endif
                                                                     px, py
#ifdef P4_TO_P8
                                                                     , pz
#endif
                              ), initial_level);
}

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

  if (!pp->params.hollow) {
    p4est_mesh_destroy (pp->mesh);
    p4est_ghost_destroy (pp->ghost);
  }

  P4EST_FREE (pp->flags);
  P4EST_FREE (pp->temp_flags);

  p4est_destroy (pp->p4est);

  /* safety checks for connectivity ownership */
  if (pp->conn_owner != NULL) {
    /* we are a copy of a wrap and have borrowed its connectivity */
    P4EST_ASSERT (!sc_refcount_is_active (&pp->conn_rc));
    P4EST_EXECUTE_ASSERT_FALSE (sc_refcount_unref (&pp->conn_owner->conn_rc));
  }
  else {
    /* we are the original wrap that owns the connectivity */
    P4EST_EXECUTE_ASSERT_TRUE (sc_refcount_unref (&pp->conn_rc));
    p4est_connectivity_destroy (pp->conn);
  }

  P4EST_FREE (pp);
}

void
p4est_wrap_set_hollow (p4est_wrap_t * pp, int hollow)
{
  /* Verify consistency */
  if (!pp->params.hollow) {
    P4EST_ASSERT (pp->flags != NULL);
    P4EST_ASSERT (pp->ghost != NULL);
    P4EST_ASSERT (pp->mesh != NULL);
  }
  else {
    P4EST_ASSERT (pp->flags == NULL);
    P4EST_ASSERT (pp->ghost == NULL);
    P4EST_ASSERT (pp->mesh == NULL);
  }

  /* Make sure a full wrap is only set to hollow outside of adaptation cycle */
  P4EST_ASSERT (!pp->match_aux);
  P4EST_ASSERT (pp->temp_flags == NULL);
  P4EST_ASSERT (pp->ghost_aux == NULL);
  P4EST_ASSERT (pp->mesh_aux == NULL);

  /* Do nothing if the status is right */
  if (hollow == pp->params.hollow) {
    return;
  }

  if (pp->params.hollow) {
    /* Allocate the ghost, mesh, and flag members */
    pp->flags = P4EST_ALLOC_ZERO (uint8_t, pp->p4est->local_num_quadrants);
    pp->ghost = p4est_ghost_new (pp->p4est, pp->params.mesh_params.btype);
    pp->mesh =
      p4est_mesh_new_params (pp->p4est, pp->ghost, &pp->params.mesh_params);
  }
  else {
    /* Free and nullify the ghost, mesh, and flag members */
    p4est_mesh_destroy (pp->mesh);
    p4est_ghost_destroy (pp->ghost);
    P4EST_FREE (pp->flags);
    pp->ghost = NULL;
    pp->mesh = NULL;
    pp->flags = NULL;
  }
  pp->num_refine_flags = pp->inside_counter = pp->num_replaced = 0;
  pp->params.hollow = hollow;
}

void
p4est_wrap_set_coarsen_delay (p4est_wrap_t * pp,
                              int coarsen_delay, int coarsen_affect)
{
  size_t              zz;
  p4est_topidx_t      tt;
  p4est_t            *p4est;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *quadrant;
  sc_array_t         *tquadrants;

  P4EST_ASSERT (pp != NULL);
  P4EST_ASSERT (coarsen_delay >= 0);

  pp->params.coarsen_delay = coarsen_delay;
  pp->params.coarsen_affect = coarsen_affect;
  p4est = pp->p4est;
  P4EST_ASSERT (p4est->data_size == 0);

  /* initialize delay memory in the quadrants' user field */
  for (tt = p4est->first_local_tree; tt <= p4est->last_local_tree; ++tt) {
    tree = p4est_tree_array_index (p4est->trees, tt);
    tquadrants = &tree->quadrants;
    for (zz = 0; zz < tquadrants->elem_count; ++zz) {
      quadrant = p4est_quadrant_array_index (tquadrants, zz);
      quadrant->p.user_int = 0;
    }
  }
}

void
p4est_wrap_set_partitioning (p4est_wrap_t *pp, int partition_for_coarsening)
{
  P4EST_ASSERT (pp != NULL);
  pp->params.partition_for_coarsening = partition_for_coarsening;
}

p4est_ghost_t      *
p4est_wrap_get_ghost (p4est_wrap_t * pp)
{
  P4EST_ASSERT (!pp->params.hollow);

  return pp->match_aux ? pp->ghost_aux : pp->ghost;
}

p4est_mesh_t       *
p4est_wrap_get_mesh (p4est_wrap_t * pp)
{
  P4EST_ASSERT (!pp->params.hollow);

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

  P4EST_ASSERT (!pp->params.hollow);
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

  P4EST_ASSERT (!pp->params.hollow);
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

  P4EST_ASSERT (!pp->params.hollow);
  P4EST_ASSERT (pp->params.coarsen_delay >= 0);

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
  p4est_coarsen_ext (p4est, 0, 1, coarsen_callback, NULL,
                     pp->params.coarsen_delay ? replace_on_coarsen :
                     pp->params.replace_fn);
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
    p4est_balance_ext (p4est, pp->params.mesh_params.btype, NULL,
                       pp->params.coarsen_delay ? replace_on_balance :
                       pp->params.replace_fn);
    pp->flags = P4EST_ALLOC_ZERO (uint8_t, p4est->local_num_quadrants);

    pp->ghost_aux = p4est_ghost_new (p4est, pp->params.mesh_params.btype);
    pp->mesh_aux =
      p4est_mesh_new_params (p4est, pp->ghost_aux, &pp->params.mesh_params);
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

static void
p4est_wrap_partition_unchanged (p4est_gloidx_t pre_me,
                                p4est_gloidx_t pre_next,
                                p4est_gloidx_t post_me,
                                p4est_gloidx_t post_next,
                                p4est_locidx_t * unchanged_first,
                                p4est_locidx_t * unchanged_length,
                                p4est_locidx_t * unchanged_old_first)
{
  p4est_locidx_t      uf, ul, uof;
  p4est_gloidx_t      unext;

  /* consistency checks */
  P4EST_ASSERT (0 <= pre_me && pre_me <= pre_next);
  P4EST_ASSERT (0 <= post_me && post_me <= post_next);

  /* initialize the case that no quadrants stay on this processor */
  uf = ul = uof = 0;

  /* check whether any quadrants stay at all, and which ones */
  if (pre_me < post_next && post_me < pre_next) {
    unext = SC_MIN (pre_next, post_next);
    if (pre_me <= post_me) {
      uof = (p4est_locidx_t) (post_me - pre_me);
      ul = (p4est_locidx_t) (unext - post_me);
    }
    else {
      uf = (p4est_locidx_t) (pre_me - post_me);
      ul = (p4est_locidx_t) (unext - pre_me);
    }
  }

  /* consistency checks */
  P4EST_ASSERT (uf >= 0 && ul >= 0 && uof >= 0);
  P4EST_ASSERT ((p4est_gloidx_t) uf + ul <= post_next - post_me);
  P4EST_ASSERT ((p4est_gloidx_t) uof + ul <= pre_next - pre_me);

  /* assign to output variables */
  if (unchanged_first != NULL) {
    *unchanged_first = uf;
  }
  if (unchanged_length != NULL) {
    *unchanged_length = ul;
  }
  if (unchanged_old_first != NULL) {
    *unchanged_old_first = uof;
  }
}

int
p4est_wrap_partition (p4est_wrap_t * pp, int weight_exponent,
                      p4est_locidx_t * unchanged_first,
                      p4est_locidx_t * unchanged_length,
                      p4est_locidx_t * unchanged_old_first)
{
  int                 changed;
  p4est_gloidx_t      pre_me, pre_next;
  p4est_gloidx_t      post_me, post_next;

  P4EST_ASSERT (!pp->params.hollow);

  P4EST_ASSERT (pp->ghost != NULL);
  P4EST_ASSERT (pp->mesh != NULL);
  P4EST_ASSERT (pp->ghost_aux != NULL);
  P4EST_ASSERT (pp->mesh_aux != NULL);
  P4EST_ASSERT (pp->match_aux == 1);

  p4est_mesh_destroy (pp->mesh);
  p4est_ghost_destroy (pp->ghost);
  pp->match_aux = 0;

  /* Remember the window onto global quadrant sequence before partition */
  pre_me = pp->p4est->global_first_quadrant[pp->p4est->mpirank];
  pre_next = pp->p4est->global_first_quadrant[pp->p4est->mpirank + 1];

  /* Initialize output for the case that the partition does not change */
  if (unchanged_first != NULL) {
    *unchanged_first = 0;
  }
  if (unchanged_length != NULL) {
    *unchanged_length = pp->p4est->local_num_quadrants;
  }
  if (unchanged_old_first != NULL) {
    *unchanged_old_first = 0;
  }

  /* In the future the flags could be used to pass partition weights */
  /* We need to lift the restriction on 64 bits for the global weight sum */
  P4EST_ASSERT (weight_exponent == 0 || weight_exponent == 1);
  pp->weight_exponent = weight_exponent;
  changed =
    p4est_partition_ext (pp->p4est, pp->params.partition_for_coarsening,
                         weight_exponent ? partition_weight : NULL) > 0;

  if (changed) {
    P4EST_FREE (pp->flags);
    pp->flags = P4EST_ALLOC_ZERO (uint8_t, pp->p4est->local_num_quadrants);

    pp->ghost = p4est_ghost_new (pp->p4est, pp->params.mesh_params.btype);
    pp->mesh =
      p4est_mesh_new_params (pp->p4est, pp->ghost, &pp->params.mesh_params);

    /* Query the window onto global quadrant sequence after partition */
    if (unchanged_first != NULL || unchanged_length != NULL ||
        unchanged_old_first != NULL) {

      /* compute new windof of local quadrants */
      post_me = pp->p4est->global_first_quadrant[pp->p4est->mpirank];
      post_next = pp->p4est->global_first_quadrant[pp->p4est->mpirank + 1];

      /* compute the range of quadrants that have stayed on this processor */
      p4est_wrap_partition_unchanged (pre_me, pre_next, post_me, post_next,
                                      unchanged_first, unchanged_length,
                                      unchanged_old_first);
    }
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
  P4EST_ASSERT (!pp->params.hollow);

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
#ifdef P4EST_ENABLE_DEBUG
  p4est_t            *p4est = leaf->pp->p4est;
#endif
#if 0
  p4est_quadrant_t    corner;
#endif
  p4est_quadrant_t   *mirror;

  /* complete information on current quadrant */
  leaf->local_quad = leaf->tree->quadrants_offset + leaf->which_quad;
  leaf->quad =
    p4est_quadrant_array_index (leaf->tquadrants, leaf->which_quad);

#if 0
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
#endif

#if 0
#ifdef P4EST_ENABLE_DEBUG
  printf ("C: Leaf level %d tree %d tree_leaf %d local_leaf %d\n",
          (int) leaf->quad->level, leaf->which_tree, leaf->which_quad,
          leaf->local_quad);
#endif
#endif

  /* track parallel mirror quadrants */
  if (leaf->mirrors != NULL) {
    if (leaf->local_quad == leaf->next_mirror_quadrant) {
      if (++leaf->nm + 1 < (p4est_locidx_t) leaf->mirrors->elem_count) {
        mirror = p4est_quadrant_array_index (leaf->mirrors, leaf->nm + 1);
        leaf->next_mirror_quadrant = mirror->p.piggy3.local_num;
        P4EST_ASSERT (leaf->next_mirror_quadrant > leaf->local_quad);
        P4EST_ASSERT (leaf->next_mirror_quadrant <
                      p4est->local_num_quadrants);
      }
      else {
        leaf->next_mirror_quadrant = -1;
      }
      leaf->is_mirror = 1;
    }
    else {
      leaf->is_mirror = 0;
    }
  }

  return leaf;
}

p4est_wrap_leaf_t  *
p4est_wrap_leaf_first (p4est_wrap_t * pp, int track_mirrors)
{
  p4est_t            *p4est = pp->p4est;
  p4est_wrap_leaf_t  *leaf;
  p4est_quadrant_t   *mirror;

  if (p4est->local_num_quadrants == 0) {
    P4EST_ASSERT (p4est->first_local_tree == -1);
    P4EST_ASSERT (p4est->last_local_tree == -2);
    return NULL;
  }

  /* prepare internal state of the leaf iterator */
  leaf = P4EST_ALLOC (p4est_wrap_leaf_t, 1);
  leaf->pp = pp;
  leaf->which_tree = p4est->first_local_tree;
  P4EST_ASSERT (leaf->which_tree >= 0);
  leaf->tree = p4est_tree_array_index (p4est->trees, leaf->which_tree);
  leaf->tquadrants = &leaf->tree->quadrants;
  P4EST_ASSERT (leaf->tquadrants->elem_size > 0);
  leaf->which_quad = 0;

  /* initialize mirror tracking if desired */
  leaf->nm = leaf->next_mirror_quadrant = -1;
  if (track_mirrors) {
    leaf->mirrors = &(p4est_wrap_get_ghost (pp))->mirrors;
    if (leaf->mirrors->elem_count > 0) {
      mirror = p4est_quadrant_array_index (leaf->mirrors, 0);
      leaf->next_mirror_quadrant = (int) mirror->p.piggy3.local_num;
      P4EST_ASSERT (leaf->next_mirror_quadrant >= 0);
      P4EST_ASSERT (leaf->next_mirror_quadrant < p4est->local_num_quadrants);
    }
  }
  else {
    leaf->mirrors = NULL;
    leaf->is_mirror = 0;
  }

  /* complete leaf and mirror information */
  return p4est_wrap_leaf_info (leaf);
}

p4est_wrap_leaf_t  *
p4est_wrap_leaf_next (p4est_wrap_leaf_t * leaf)
{
  p4est_t            *p4est = leaf->pp->p4est;

  P4EST_ASSERT (leaf != NULL);

  if ((size_t) leaf->which_quad + 1 == leaf->tquadrants->elem_count) {
    ++leaf->which_tree;
    if (leaf->which_tree > p4est->last_local_tree) {
#ifdef P4EST_ENABLE_DEBUG
      if (leaf->mirrors != NULL) {
        P4EST_ASSERT (leaf->nm + 1 ==
                      (p4est_locidx_t) leaf->mirrors->elem_count);
        P4EST_ASSERT (leaf->next_mirror_quadrant == -1);
      }
#endif
      P4EST_FREE (leaf);
      return NULL;
    }
    leaf->tree = p4est_tree_array_index (p4est->trees, leaf->which_tree);
    leaf->tquadrants = &leaf->tree->quadrants;
    P4EST_ASSERT (leaf->tquadrants->elem_size > 0);
    leaf->which_quad = 0;
  }
  else {
    ++leaf->which_quad;
  }

  return p4est_wrap_leaf_info (leaf);
}
