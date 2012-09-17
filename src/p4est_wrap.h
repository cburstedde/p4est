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

#include <p4est_mesh.h>

SC_EXTERN_C_BEGIN;

/*** COMPLETE INTERNAL STATE OF P4EST ***/

typedef enum p4est_wrap_flags {
  P4EST_WRAP_REFINE = 1,
  P4EST_WRAP_COARSEN = 2
}
p4est_wrap_flags_t;

typedef struct p4est_wrap
{
  int                 p4est_dim;
  int                 p4est_half;
  int                 p4est_faces;
  int                 p4est_children;
  p4est_connectivity_t *conn;
  p4est_t            *p4est;
  int8_t             *flags;
  p4est_ghost_t      *ghost;
  p4est_mesh_t       *mesh;
  p4est_ghost_t      *ghost_aux;
  p4est_mesh_t       *mesh_aux;
  int                 match_aux;
}
p4est_wrap_t;

/** Create p4est and auxiliary data structures.
 * Expects MPI_Init to be called beforehand.
 * The pp->flags array is initialized to all zeros.
 */
p4est_wrap_t       *p4est_wrap_new_unitsquare (MPI_Comm mpicomm,
                                               int initial_level);
p4est_wrap_t       *p4est_wrap_new_periodic (MPI_Comm mpicomm,
                                             int initial_level);
p4est_wrap_t       *p4est_wrap_new_pillow (MPI_Comm mpicomm,
                                           int initial_level);
p4est_wrap_t       *p4est_wrap_new_moebius (MPI_Comm mpicomm,
                                            int initial_level);

/** Passes MPI_COMM_WORLD to p4est_wrap_new_unitsquare. */
p4est_wrap_t       *p4est_wrap_new_world (int initial_level);
void                p4est_wrap_destroy (p4est_wrap_t * pp);

/** Call p4est_refine, coarsen, balance to update pp->p4est.
 * Checks pp->flags as per-quadrant input against p4est_wrap_flags_t.
 * The pp->flags array is updated along with p4est and initialized to zeros.
 * Creates ghost_aux and mesh_aux to represent the intermediate mesh.
 * \return          boolean whether p4est has changed.
 *                  If true, partition must be called.
 *                  If false, partition must not be called.
 */
int                 p4est_wrap_refine (p4est_wrap_t * pp);

/** Call p4est_partition for equal leaf distribution.
 * Frees the old ghost and mesh first and updates pp->flags along with p4est.
 * Creates ghost and mesh to represent the new mesh.
 * \return          boolean whether p4est has changed.
 *                  If true, complete must be called.
 *                  If false, complete must not be called.
 */
int                 p4est_wrap_partition (p4est_wrap_t * pp);

/** Free memory for the intermediate mesh.
 * Sets mesh_aux and ghost_aux to NULL.
 */
void                p4est_wrap_complete (p4est_wrap_t * pp);

/*** ITERATOR OVER THE FOREST LEAVES ***/

typedef struct p4est_wrap_leaf
{
  p4est_wrap_t       *pp;
  int                 level;
  p4est_topidx_t      which_tree;
  p4est_locidx_t      which_quad;
  p4est_locidx_t      total_quad;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *quad;
  double              lowerleft[3];
  double              upperright[3];
}
p4est_wrap_leaf_t;

/* Create an iterator over the leaves in the forest.
 * Returns a newly allocated state containing the first leaf,
 * or NULL if the local partition of the tree is empty.
 */
p4est_wrap_leaf_t  *p4est_wrap_leaf_first (p4est_wrap_t * pp);

/* Move the forest leaf iterator forward.
 * Returns the state that was input with information for the next leaf,
 * or NULL and deallocates the input if called with the last leaf.
 */
p4est_wrap_leaf_t  *p4est_wrap_leaf_next (p4est_wrap_leaf_t * leaf);

SC_EXTERN_C_END;
