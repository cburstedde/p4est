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

/** \file p4est_ghost_ext.h
 *
 * Interface to create ghost layer for 2D forest using recursive top-down tree traversal.
 *
 * \ingroup p4est
 */

#ifndef P4EST_GHOST_EXT_H
#define P4EST_GHOST_EXT_H

#include <p4est_ghost.h>

SC_EXTERN_C_BEGIN;

/** Flags for operating with unbalanced forest. */
typedef enum
{
  P4EST_GHOST_UNBALANCED_ABORT = 0,
  P4EST_GHOST_UNBALANCED_FAIL,
  P4EST_GHOST_UNBALANCED_ALLOW
}
p4est_ghost_tolerance_t;

#ifdef P4EST_ENABLE_MPI

/** Data structure that contains temporary mirror information */
typedef struct p4est_ghost_mirror
{
  int                 mpisize, mpirank;
  int                 known;    /* Flag for avoiding duplicate entries. */
  p4est_locidx_t      sum_all_procs;    /* sum of mirrors by processor */
  sc_array_t         *send_bufs;        /* Array of quadrants which are ghosts to other partitions. */
  sc_array_t         *mirrors;  /* mirror quadrants. */
  sc_array_t         *offsets_by_proc;  /* a p4est_locidx_t array per proc */
}
p4est_ghost_mirror_t;

/** Initializes the ghost mirror structure from an allocated/static ghost structure.
 * \param [in]      ghost      The ghost structure.
 * \param [in]      mpirank    Partition rank.
 * \param [in]      send_bufs  Array of ghost quadrants that need to be communicated.
 * \param [in,out]  m          The ghost mirror structure that needs to be initialized.
 */
void                p4est_ghost_mirror_init (p4est_ghost_t * ghost,
                                             int mpirank,
                                             sc_array_t * send_bufs,
                                             p4est_ghost_mirror_t * m);

/** Potentially record a quadrant that is to be sent as a mirror.
 * \param [in] m      The temporary ghost mirror data structure to work on.
 * \param [in] treeid The tree number looped through by the current rank.
 * \param [in] q      The quadrant currently looked at by current rank.
 * \param [in] p      The rank that \a q should be sent to.
 */
void                p4est_ghost_mirror_add (p4est_ghost_mirror_t * m,
                                            p4est_topidx_t treeid,
                                            p4est_locidx_t number,
                                            p4est_quadrant_t * q, int p);

/** Populate the mirror fields in the ghost layer with final data.
 * The elements in the temporary p4est_ghost_mirror_t structure are freed.
 * \param [in,out]  ghost     The recepient ghost structure.
 * \param [in]      m         The ghost mirror structure which holds the final data.
 * \param [in]      populate  Flag to enable copying data from \a m to \a ghost.
 */
void                p4est_ghost_mirror_reset (p4est_ghost_t * ghost,
                                              p4est_ghost_mirror_t * m,
                                              int populate);

/** Test if a quadrant can be a potential ghost irrespective of unbalanced/balanced forest.
 * If true then make a call to p4est_ghost_mirror_add to update the ghost mirror.
 * \param [in]      p4est       The forest structure.
 * \param [in,out]  m           The ghost mirror structure which needs to be updated.
 * \param [in]      q           The quadrant that needs to be tested.
 * \param [in]      t           The tree index of the quadrant.
 * \param [in]      nq          The neighbor quadrant to \a q.
 * \param [in]      nt          The tree index of the neighbor quadrant.
 * \param [in]      touch       Flag to check if this quadrant touches forest boundary.
 * \param [in]      rank        The rank of the owner partition of \a q
 * \param [in]      local_num   The linear morton index of \a q.
*/
void                p4est_ghost_test_add (p4est_t * p4est,
                                          p4est_ghost_mirror_t * m,
                                          p4est_quadrant_t * q,
                                          p4est_topidx_t t,
                                          p4est_quadrant_t * nq,
                                          p4est_topidx_t nt, int32_t touch,
                                          int rank,
#if 0
                                          sc_array_t * send_bufs,
#endif
                                          p4est_locidx_t local_num);

/** This adds a quadrant to the end of a buffer.
 *
 * It crams the tree id into the user_data field of the quadrant in
 * the buffer and only adds the quadrant to the end of the buffer if
 * it is unique.
 *
 * \param [in,out] buf    \a q is added to the end if it is not already there.
 * \param [in,out] q      the quadrant to be added.  The \c user_data field
 *                        is filled with \a treeid.
 * \param [in]            treeid the tree id of \a q.
 * \return                true if the ghost was added, false if duplicate.
 */
int                 p4est_ghost_add_to_buf (sc_array_t * buf,
                                            p4est_topidx_t treeid,
                                            p4est_locidx_t number,
                                            const p4est_quadrant_t * q);

/** Obtain the 1D array at a specified row index in a 2D array.
 * \param [in]  array   The 2D array.
 * \param [in]  i       The row index.
 * \return              The pointer to 1D array at row \a i.
 */
sc_array_t         *p4est_ghost_array_index (sc_array_t * array, int i);

#endif

/** Get the owner tree index of a quadrant from an quadrant array.
 * \param [in]  array   The array which stores quadrants.
 * \param [in]  zindex  The index of quadrant.
 * \return              The owner tree index for quadrant at index \a zindex.
*/
size_t              p4est_ghost_tree_type (sc_array_t * array, size_t zindex,
                                           void *data);

/** Checks if a quadrant's face is on the boundary of the forest.
 *
 * \param [in] p4est  The forest in which to search for \a q
 * \param [in] treeid The tree to which \a q belongs.
 * \param [in] q      The quadrant that is in question.
 * \param [in] face   The face of the quadrant that is in question.
 *
 * \return true if the quadrant's face is on the boundary of the forest and
 *         false otherwise.
 */
int                 p4est_quadrant_on_face_boundary (p4est_t * p4est,
                                                     p4est_topidx_t treeid,
                                                     int face,
                                                     const p4est_quadrant_t *
                                                     q);

SC_EXTERN_C_END;

#endif
