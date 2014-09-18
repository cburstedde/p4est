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

#ifndef P4EST_WRAP_H
#define P4EST_WRAP_H

#include <p4est_mesh.h>

SC_EXTERN_C_BEGIN;

/*** COMPLETE INTERNAL STATE OF P4EST ***/

typedef enum p4est_wrap_flags
{
  P4EST_WRAP_NONE = 0,
  P4EST_WRAP_REFINE = 0x01,
  P4EST_WRAP_COARSEN = 0x02
}
p4est_wrap_flags_t;

typedef struct p4est_wrap
{
  /* this member is never used or changed by p4est_wrap */
  void               *user_pointer;     /**< Convenience member for users */

  /* these members are considered public and read-only */
  int                 p4est_dim;
  int                 p4est_half;
  int                 p4est_faces;
  int                 p4est_children;
  p4est_connectivity_t *conn;
  p4est_t            *p4est;    /**< p4est->user_pointer is used internally */

  /* anything below here is considered private und should not be touched */
  int                 weight_exponent;
  uint8_t            *flags, *temp_flags;
  p4est_locidx_t      num_refine_flags, inside_counter, num_replaced;

  /* for ghost and mesh use p4est_wrap_get_ghost, _mesh declared below */
  p4est_ghost_t      *ghost;
  p4est_mesh_t       *mesh;
  p4est_ghost_t      *ghost_aux;
  p4est_mesh_t       *mesh_aux;
  int                 match_aux;
}
p4est_wrap_t;

/** Create a p4est wrapper from a given connectivity structure.
 * \param [in] mpicomm        We expect sc_MPI_Init to be called already.
 * \param [in] conn           Connectivity structure.  Wrap takes ownership.
 * \param [in] initial_level  Initial level of uniform refinement.
 * \return                    A fully initialized p4est_wrap structure.
 */
p4est_wrap_t       *p4est_wrap_new_conn (sc_MPI_Comm mpicomm,
                                         p4est_connectivity_t * conn,
                                         int initial_level);

/** Create p4est and auxiliary data structures.
 * Expects sc_MPI_Init to be called beforehand.
 */
p4est_wrap_t       *p4est_wrap_new_unitsquare (sc_MPI_Comm mpicomm,
                                               int initial_level);
p4est_wrap_t       *p4est_wrap_new_periodic (sc_MPI_Comm mpicomm,
                                             int initial_level);
p4est_wrap_t       *p4est_wrap_new_rotwrap (sc_MPI_Comm mpicomm,
                                            int initial_level);
p4est_wrap_t       *p4est_wrap_new_corner (sc_MPI_Comm mpicomm,
                                           int initial_level);
p4est_wrap_t       *p4est_wrap_new_pillow (sc_MPI_Comm mpicomm,
                                           int initial_level);
p4est_wrap_t       *p4est_wrap_new_moebius (sc_MPI_Comm mpicomm,
                                            int initial_level);
p4est_wrap_t       *p4est_wrap_new_cubed (sc_MPI_Comm mpicomm,
                                          int initial_level);
p4est_wrap_t       *p4est_wrap_new_disk (sc_MPI_Comm mpicomm,
                                         int initial_level);

/** Passes sc_MPI_COMM_WORLD to p4est_wrap_new_unitsquare. */
p4est_wrap_t       *p4est_wrap_new_world (int initial_level);
void                p4est_wrap_destroy (p4est_wrap_t * pp);

/** Return the appropriate ghost layer.
 * This function is necessary since two versions may exist simultaneously
 * after refinement and before partition/complete.
 * */
p4est_ghost_t      *p4est_wrap_get_ghost (p4est_wrap_t * pp);

/** Return the appropriate mesh structure.
 * This function is necessary since two versions may exist simultaneously
 * after refinement and before partition/complete.
 * */
p4est_mesh_t       *p4est_wrap_get_mesh (p4est_wrap_t * pp);

/** Mark a local element for refinement.
 * This will cancel any coarsening mark set previously for this element.
 * \param [in,out] wrap The p4est wrapper to work with.
 * \param [in] which_tree The number of the tree this element lives in.
 * \param [in] which_quad The number of this element relative to its tree.
 */
void                p4est_wrap_mark_refine (p4est_wrap_t * pp,
                                            p4est_topidx_t which_tree,
                                            p4est_locidx_t which_quad);

/** Mark a local element for coarsening.
 * This will cancel any refinement mark set previously for this element.
 * \param [in,out] wrap The p4est wrapper to work with.
 * \param [in] which_tree The number of the tree this element lives in.
 * \param [in] which_quad The number of this element relative to its tree.
 */
void                p4est_wrap_mark_coarsen (p4est_wrap_t * pp,
                                             p4est_topidx_t which_tree,
                                             p4est_locidx_t which_quad);

/** Call p4est_refine, coarsen, and balance to update pp->p4est.
 * Checks pp->flags as per-quadrant input against p4est_wrap_flags_t.
 * The pp->flags array is updated along with p4est and reset to zeros.
 * Creates ghost_aux and mesh_aux to represent the intermediate mesh.
 * \return          boolean whether p4est has changed.
 *                  If true, partition must be called.
 *                  If false, partition must not be called, and
 *                  complete must not be called either.
 */
int                 p4est_wrap_adapt (p4est_wrap_t * pp);

/** Call p4est_partition for equal leaf distribution.
 * Frees the old ghost and mesh first and updates pp->flags along with p4est.
 * The pp->flags array is reset to zeros.
 * Creates ghost and mesh to represent the new mesh.
 * \param [in] weight_exponent      Integer weight assigned to each leaf
 *                  according to 2 ** (level * exponent).  Passing 0 assigns
 *                  equal weight to all leaves.  Passing 1 increases the
 *                  leaf weight by a factor of two for each level increase.
 *                  CURRENTLY ONLY 0 AND 1 ARE LEGAL VALUES.
 * \return          boolean whether p4est has changed.
 *                  If true, complete must be called.
 *                  If false, complete must not be called.
 */
int                 p4est_wrap_partition (p4est_wrap_t * pp,
                                          int weight_exponent);

/** Free memory for the intermediate mesh.
 * Sets mesh_aux and ghost_aux to NULL.
 * This function must be used if both refinement and partition effect changes.
 * After this call, we are ready for another mark-refine-partition cycle.
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

#endif /* !P4EST_WRAP_H */
