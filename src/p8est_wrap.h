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

#ifndef P8EST_WRAP_H
#define P8EST_WRAP_H

/** \file p8est_wrap.h
 * The logic in p8est_wrap encapsulates core p4est data structures and provides
 * functions that clarify the mark-adapt-partition cycle.  There is also an
 * element iterator that can replace the nested loops over trees and tree
 * quadrants, respectively, which can help make application code cleaner.
 */

#include <p8est_mesh.h>
#include <p8est_extended.h>
#include <sc_refcount.h>

SC_EXTERN_C_BEGIN;

/*** COMPLETE INTERNAL STATE OF P8EST ***/

typedef enum p8est_wrap_flags
{
  P8EST_WRAP_NONE = 0,
  P8EST_WRAP_REFINE = 0x01,
  P8EST_WRAP_COARSEN = 0x02
}
p8est_wrap_flags_t;

typedef struct p8est_wrap
{
  /* this member is never used or changed by p8est_wrap */
  void               *user_pointer;     /**< Convenience member for users */

  /** If true, this wrap has NULL for ghost, mesh, and flag members.
   * If false, they are properly allocated and kept current internally. */
  int                 hollow;

  /** Non-negative integer tells us how many adaptations to wait
   * before any given quadrant may be coarsened again. */
  int                 coarsen_delay;

  /** Boolean: If true, we delay coarsening not only after refinement,
   * but also between subsequent coarsenings of the same quadrant. */
  int                 coarsen_affect;

  /** This reference counter is a workaround for internal use only.
   * Until we have refcounting/copy-on-write for the connectivity,
   * we count the references to conn by copies of this wrap structure.
   * There must be no external references left when this wrap is destroyed.
   */
  sc_refcount_t       conn_rc;
  p8est_connectivity_t *conn;
  struct p8est_wrap  *conn_owner;

  /* these members are considered public and read-only */
  int                 p4est_dim;
  int                 p4est_half;
  int                 p4est_faces;
  int                 p4est_children;
  p8est_connect_type_t btype;
  p8est_replace_t     replace_fn;
  p8est_t            *p4est;    /**< p4est->user_pointer is used internally */

  /* anything below here is considered private und should not be touched */
  int                 weight_exponent;
  uint8_t            *flags, *temp_flags;
  p4est_locidx_t      num_refine_flags, inside_counter, num_replaced;

  /* for ghost and mesh use p8est_wrap_get_ghost, _mesh declared below */
  p8est_ghost_t      *ghost;
  p8est_mesh_t       *mesh;
  p8est_ghost_t      *ghost_aux;
  p8est_mesh_t       *mesh_aux;
  int                 match_aux;
}
p8est_wrap_t;

/** Create a p8est wrapper from a given connectivity structure.
 * The ghost and mesh members are initialized as well as the flags.
 * The btype is set to P8EST_CONNECT_FULL.
 * \param [in] mpicomm        We expect sc_MPI_Init to be called already.
 * \param [in] conn           Connectivity structure.  Wrap takes ownership.
 * \param [in] initial_level  Initial level of uniform refinement.
 * \return                    A fully initialized p8est_wrap structure.
 */
p8est_wrap_t       *p8est_wrap_new_conn (sc_MPI_Comm mpicomm,
                                         p8est_connectivity_t * conn,
                                         int initial_level);

/** Create a wrapper for a given p8est structure.
 * \param [in,out] p8est      Valid p8est object that we will own.
 *                            We take ownership of its connectivity too.
 *                            Its user pointer must be NULL and will be changed.
 * \param [in] hollow         Do not allocate flags, ghost, and mesh members.
 * \param [in] btype          The neighborhood used for balance, ghost, mesh.
 * \param [in] replace_fn     Callback to replace quadrants during refinement,
 *                            coarsening or balancing in \ref p8est_wrap_adapt.
 *                            May be NULL.
 * \param [in] user_pointer   Set the user pointer in \ref p8est_wrap_t.
 *                            Subsequently, we will never access it.
 * \return                    A fully initialized p8est_wrap structure.
 */
p8est_wrap_t       *p8est_wrap_new_p8est (p8est_t * p8est, int hollow,
                                          p8est_connect_type_t btype,
                                          p8est_replace_t replace_fn,
                                          void *user_pointer);

/** Create a p8est wrapper from a given connectivity structure.
 * Like p8est_wrap_new_conn, but with extra parameters \a hollow and \a btype.
 * \param [in] mpicomm        We expect sc_MPI_Init to be called already.
 * \param [in] conn           Connectivity structure.  Wrap takes ownership.
 * \param [in] initial_level  Initial level of uniform refinement.
 *                            No effect if less/equal to zero.
 * \param [in] hollow         Do not allocate flags, ghost, and mesh members.
 * \param [in] btype          The neighborhood used for balance, ghost, mesh.
 * \param [in] replace_fn     Callback to replace quadrants during refinement,
 *                            coarsening or balancing in \ref p8est_wrap_adapt.
 *                            May be NULL.
 * \param [in] user_pointer   Set the user pointer in \ref p8est_wrap_t.
 *                            Subsequently, we will never access it.
 * \return                    A fully initialized p8est_wrap structure.
 */
p8est_wrap_t       *p8est_wrap_new_ext (sc_MPI_Comm mpicomm,
                                        p8est_connectivity_t * conn,
                                        int initial_level, int hollow,
                                        p8est_connect_type_t btype,
                                        p8est_replace_t replace_fn,
                                        void *user_pointer);

/** Create a p8est wrapper from an existing one.
 * \note This wrapper must be destroyed before the original one.
 * We set it to hollow and copy the original p8est data structure.
 * \param [in,out] source   We access the source for debugging purposes.
 * \param [in] data_size    The data size installed in the copied forest.
 * \param [in] replace_fn     Callback to replace quadrants during refinement,
 *                            coarsening or balancing in \ref p8est_wrap_adapt.
 *                            May be NULL.
 * \param [in] user_pointer   Set the user pointer in \ref p8est_wrap_t.
 *                            Subsequently, we will never access it.
 * \return                    A fully initialized p8est_wrap structure.
 */
p8est_wrap_t       *p8est_wrap_new_copy (p8est_wrap_t * source,
                                         size_t data_size,
                                         p8est_replace_t replace_fn,
                                         void *user_pointer);

/** Create p8est and auxiliary data structures.
 * Expects sc_MPI_Init to be called beforehand.
 */
p8est_wrap_t       *p8est_wrap_new_unitcube (sc_MPI_Comm mpicomm,
                                             int initial_level);
p8est_wrap_t       *p8est_wrap_new_rotwrap (sc_MPI_Comm mpicomm,
                                            int initial_level);
p8est_wrap_t       *p8est_wrap_new_brick (sc_MPI_Comm mpicomm,
                                          int bx, int by, int bz,
                                          int px, int py, int pz,
                                          int initial_level);

/** Passes sc_MPI_COMM_WORLD to p8est_wrap_new_unitcube. */
p8est_wrap_t       *p8est_wrap_new_world (int initial_level);
void                p8est_wrap_destroy (p8est_wrap_t * pp);

/** Change hollow status of the wrap.
 * It is legal to set to the current hollow status.
 * \param [in,out] pp   The present wrap structure, hollow or not.
 * \param [in] hollow   The desired hollow status.
 */
void                p8est_wrap_set_hollow (p8est_wrap_t * pp, int hollow);

/** Set a parameter that delays coarsening after adaptation.
 * If positive each quadrant counts the number of adaptations it has survived.
 * Calling this function initializes all quadrant counters to zero.
 * On adaptation we only coarsen a quadrant if it is old enough.
 * Optionally, we can also delay the time between subsequent coarsenings.
 * \param [in,out] pp           A valid p8est_wrap structure.
 * \param [in] coarsen_delay    Set how many adaptation cycles a quadrant has
 *                              to wait to be allowed to coarsen.
 *                              Non-negative number; 0 disables the feature.
 *                              Suggested default value: not larger than 2.
 * \param [in] coarsen_affect   Boolean; If true, we not only count from the
 *                              most recent refinement but also between
 *                              subsequent coarsenings.
 *                              Suggested default: 0.
 */
void                p8est_wrap_set_coarsen_delay (p8est_wrap_t * pp,
                                                  int coarsen_delay,
                                                  int coarsen_affect);

/** Return the appropriate ghost layer.
 * This function is necessary since two versions may exist simultaneously
 * after refinement and before partition/complete.
 * \param [in,out] pp The p8est wrapper to work with, must not be hollow.
 * */
p8est_ghost_t      *p8est_wrap_get_ghost (p8est_wrap_t * pp);

/** Return the appropriate mesh structure.
 * This function is necessary since two versions may exist simultaneously
 * after refinement and before partition/complete.
 * \param [in,out] pp The p8est wrapper to work with, must not be hollow.
 * */
p8est_mesh_t       *p8est_wrap_get_mesh (p8est_wrap_t * pp);

/** Mark a local element for refinement.
 * This will cancel any coarsening mark set previously for this element.
 * \param [in,out] pp The p8est wrapper to work with, must not be hollow.
 * \param [in] which_tree The number of the tree this element lives in.
 * \param [in] which_quad The number of this element relative to its tree.
 */
void                p8est_wrap_mark_refine (p8est_wrap_t * pp,
                                            p4est_topidx_t which_tree,
                                            p4est_locidx_t which_quad);

/** Mark a local element for coarsening.
 * This will cancel any refinement mark set previously for this element.
 * \param [in,out] pp The p8est wrapper to work with, must not be hollow.
 * \param [in] which_tree The number of the tree this element lives in.
 * \param [in] which_quad The number of this element relative to its tree.
 */
void                p8est_wrap_mark_coarsen (p8est_wrap_t * pp,
                                             p4est_topidx_t which_tree,
                                             p4est_locidx_t which_quad);

/** Call p8est_refine, coarsen, and balance to update pp->p8est.
 * Checks pp->flags as per-quadrant input against p8est_wrap_flags_t.
 * The pp->flags array is updated along with p8est and reset to zeros.
 * Creates ghost_aux and mesh_aux to represent the intermediate mesh.
 * \param [in,out] pp The p8est wrapper to work with, must not be hollow.
 * \return          boolean whether p8est has changed.
 *                  If true, partition must be called.
 *                  If false, partition must not be called,
 *                  and complete must not be called either.
 */
int                 p8est_wrap_adapt (p8est_wrap_t * pp);

/** Call p8est_partition for equal leaf distribution.
 * Frees the old ghost and mesh first and updates pp->flags along with p8est.
 * Creates ghost and mesh to represent the new mesh.
 * \param [in,out] pp The p8est wrapper to work with, must not be hollow.
 * \param [in] weight_exponent      Integer weight assigned to each leaf
 *                  according to 2 ** (level * exponent).  Passing 0 assigns
 *                  equal weight to all leaves.  Passing 1 increases the
 *                  leaf weight by a factor of two for each level increase.
 *                  CURRENTLY ONLY 0 AND 1 ARE LEGAL VALUES.
 * \param [out] unchanged_first
 *                  If not NULL, is assigned the processor-local index of the
 *                  first local quadrant that has stayed on this processor.  If
 *                  no quadrant has stayed, the value is set to zero.
 *                  This number is in reference to the new (output) partition.
 * \param [out] unchanged_length
 *                  If not NULL, is assigned the number of quadrants that have
 *                  stayed on this processor.  If no quadrant has stayed, the
 *                  value is set to zero.
 * \param [out] unchanged_old_first
 *                  If not NULL, is assigned the processor-local index of the
 *                  first local quadrant that has stayed with reference to
 *                  the old (input) partition.  If no quadrant has stayed,
 *                  the value is set to zero.
 * \return          boolean whether p8est has changed.
 *                  If true, complete must be called.
 *                  If false, complete must not be called.
 */
int                 p8est_wrap_partition (p8est_wrap_t * pp,
                                          int weight_exponent,
                                          p4est_locidx_t * unchanged_first,
                                          p4est_locidx_t * unchanged_length,
                                          p4est_locidx_t *
                                          unchanged_old_first);

/** Free memory for the intermediate mesh.
 * Sets mesh_aux and ghost_aux to NULL.
 * This function must be used if both refinement and partition effect changes.
 * After this call, we are ready for another mark-refine-partition cycle.
 * \param [in,out] pp The p8est wrapper to work with, must not be hollow.
 */
void                p8est_wrap_complete (p8est_wrap_t * pp);

/*** ITERATOR OVER THE FOREST LEAVES ***/

typedef struct p8est_wrap_leaf
{
  p8est_wrap_t       *pp;             /**< Must contain a valid ghost */

  /* Information about the current quadrant */
  p4est_topidx_t      which_tree;     /**< Current tree number */
  p4est_locidx_t      which_quad;     /**< Quadrant number relative to tree */
  p4est_locidx_t      local_quad;     /**< Quadrant number relative to proc */
  p8est_tree_t       *tree;           /**< Current tree */
  sc_array_t         *tquadrants;     /**< Current tree's quadrants */
  p8est_quadrant_t   *quad;           /**< Current quadrant */
#if 0                           /* DEPRECATED -- anyone using them? */
  int                 level;
  double              lowerleft[3];
  double              upperright[3];
#endif

  /* Information about parallel neighbors */
  int                 is_mirror;      /**< Quadrant at parallel boundary? */
  sc_array_t         *mirrors;        /**< If not NULL, from pp's ghost */
  p4est_locidx_t      nm;             /**< Internal: mirror counter */
  p4est_locidx_t      next_mirror_quadrant;     /**< Internal: next */
}
p8est_wrap_leaf_t;

/** Determine whether we have just entered a different tree */
#define P8EST_LEAF_IS_FIRST_IN_TREE(wleaf) ((wleaf)->which_quad == 0)

/* Create an iterator over the local leaves in the forest.
 * Returns a newly allocated state containing the first leaf,
 * or NULL if the local partition of the tree is empty.
 * \param [in] pp   Legal p8est_wrap structure, hollow or not.
 * \param [in] track_mirrors    If true, \a pp must not be hollow and mirror
 *                              information from the ghost layer is stored.
 * \return          NULL if processor is empty, otherwise a leaf iterator for
 *                  subsequent use with \a p8est_wrap_leaf_next.
 */
p8est_wrap_leaf_t  *p8est_wrap_leaf_first (p8est_wrap_t * pp,
                                           int track_mirrors);

/* Move the forest leaf iterator forward.
 * \param [in,out] leaf     A non-NULL leaf iterator created by
 *                          \ref p8est_wrap_leaf_first.
 * \return          The state that was input with updated information for the
 *                  next leaf, or NULL and deallocates the input if called with
 *                  the last leaf on this processor.
 */
p8est_wrap_leaf_t  *p8est_wrap_leaf_next (p8est_wrap_leaf_t * leaf);

SC_EXTERN_C_END;

#endif /* !P8EST_WRAP_H */
