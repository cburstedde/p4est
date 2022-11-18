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
 * passing quadrants and data to neighboring processes
 *
 * \ingroup p4est
 */

#ifndef P4EST_GHOST_EXT_H
#define P4EST_GHOST_EXT_H

#include <p4est_ghost.h>

SC_EXTERN_C_BEGIN;

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
  int                 known;    /* was this mirror added before? */
  p4est_locidx_t      sum_all_procs;    /* sum of mirrors by processor */
  sc_array_t         *send_bufs;        /* lives in p4est_ghost_new_check */
  sc_array_t         *mirrors;  /* lives in p4est_ghost_t */
  sc_array_t         *offsets_by_proc;  /* a p4est_locidx_t array per proc */
}
p4est_ghost_mirror_t;

void                p4est_ghost_mirror_init (p4est_ghost_t * ghost,
                                             int mpirank,
                                             sc_array_t * send_bufs,
                                             p4est_ghost_mirror_t * m);

/** Potentially record a quadrant that is to be sent as a mirror
 * \param [in] m      The temporary data structure to work on.
 * \param [in] treeid The tree number looped through by the current rank.
 * \param [in] q      The quadrant currently looked at by current rank.
 * \param [in] p      The rank that \a q should be sent to.
 */
void                p4est_ghost_mirror_add (p4est_ghost_mirror_t * m,
                                            p4est_topidx_t treeid,
                                            p4est_locidx_t number,
                                            p4est_quadrant_t * q, int p);

/** Populate the mirror fields in the ghost layer with final data.
 * The elements in the temporary p4est_ghost_mirror_t structure are freed. */
void                p4est_ghost_mirror_reset (p4est_ghost_t * ghost,
                                              p4est_ghost_mirror_t * m,
                                              int populate);

/** */
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
int                 p4est_add_ghost_to_buf (sc_array_t * buf,
                                            p4est_topidx_t treeid,
                                            p4est_locidx_t number,
                                            const p4est_quadrant_t * q);

/* Helper routines */
sc_array_t         *p4est_ghost_array_index (sc_array_t * array, int i);

#endif

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
